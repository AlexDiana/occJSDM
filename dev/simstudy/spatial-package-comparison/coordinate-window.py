#!/usr/bin/env python3
"""Run an authorised window or completion queue through the frozen scientific CLI."""
import fcntl, hashlib, json, math, os, pathlib, signal, subprocess, sys, time
from launch import atomic_json


def run_window(launcher, script, run, policy_path, record_path):
    root=pathlib.Path(run).resolve(); policy_path=pathlib.Path(policy_path).resolve()
    policy=json.loads(policy_path.read_text()); record=pathlib.Path(record_path).resolve()
    if pathlib.Path(policy['run_root']).resolve()!=root: raise RuntimeError('Window run identity mismatch')
    mode=policy.get('mode','bounded_window')
    if mode not in ('bounded_window','run_to_completion'): raise RuntimeError('Unknown scheduling mode')
    bounded=mode=='bounded_window';workers=policy['workers']
    if type(workers) is not int or not 1<=workers<=(6 if bounded else 8): raise RuntimeError('Invalid worker count')
    if bounded:
        cutoff=policy['dispatch_until'];deadline=policy['stop_by']
        if not math.isfinite(deadline) or not math.isfinite(cutoff) or cutoff>deadline:
            raise RuntimeError('Invalid deadline')
    else:
        if 'dispatch_until' in policy or 'stop_by' in policy:
            raise RuntimeError('Completion mode must not carry window deadlines')
        cutoff=deadline=float('inf')
    recorded_deadline=deadline if bounded else None
    prefix='window' if bounded else 'fit'
    grace=float(policy.get('grace_seconds',10)); poll=float(policy.get('poll_seconds',1))
    if not 0<grace<=10 or not 0<poll<=1: raise RuntimeError('Invalid cleanup timing')
    def verify_control():
        for name,digest in policy['files_sha256'].items():
            if hashlib.sha256(pathlib.Path(name).read_bytes()).hexdigest()!=digest:
                raise RuntimeError('Window source identity mismatch: '+name)
    verify_control()
    control=root/'control';control.mkdir(exist_ok=True)
    lock=(control/'dispatch.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    record.mkdir(parents=True,exist_ok=False)
    atomic_json(record/'owner.json',dict(pid=os.getpid(),started=time.time(),policy=str(policy_path),policy_sha256=hashlib.sha256(policy_path.read_bytes()).hexdigest()))
    active={}; queue_proc=None; queue_logs=None; completed=[]; interrupted=[]; deferred=[]
    outcome='paused_window'; queue_interrupted=False; error=None
    def update(stage,**extra):
        atomic_json(root/'continuation-status.json',dict(stage=stage,updated=time.time(),amendment=policy.get('amendment','E14'),mode=mode,window_record=str(record),deadline=recorded_deadline,**extra))
    def stop_signal(signum,_): raise SystemExit(128+signum)
    old_handlers={s:signal.signal(s,stop_signal) for s in (signal.SIGTERM,signal.SIGINT)}
    def group_members(proc):
        rows=subprocess.check_output(['/bin/ps','-axo','pid=,pgid=,stat='],text=True)
        members=[]
        for line in rows.splitlines():
            fields=line.split()
            if len(fields)==3 and int(fields[1])==proc.pid and 'Z' not in fields[2]:
                members.append(int(fields[0]))
        if proc.returncode is not None and proc.pid in members:
            raise RuntimeError('Reaped leader PID was reused; refuse to signal it')
        return members
    def terminate(processes):
        # Each child starts a new session. Retain its group until all live
        # members disappear, even when the leader exits before descendants.
        for p in processes:
            if group_members(p):
                try: os.killpg(p.pid,signal.SIGTERM)
                except ProcessLookupError: pass
        until=time.monotonic()+grace
        while any(group_members(p) for p in processes) and time.monotonic()<until:
            time.sleep(min(.05,poll))
        for p in processes:
            if group_members(p):
                try: os.killpg(p.pid,signal.SIGKILL)
                except ProcessLookupError: pass
            p.wait()
        until=time.monotonic()+2
        while any(group_members(p) for p in processes) and time.monotonic()<until:
            time.sleep(min(.05,poll))
        if any(group_members(p) for p in processes):
            raise RuntimeError('Owned process group has not exited; do not declare drive quiet')
    def finish(pid,status='exited',reason=None):
        proc,task,logs,started=active[pid]
        terminate([proc])
        del active[pid]
        for log in logs: log.close()
        e=dict(exit_status=proc.returncode,elapsed_seconds=time.time()-started,status=status)
        if reason: e['reason']=reason
        atomic_json(pathlib.Path(task['directory'])/'exit.json',e)
        (interrupted if status=='resource_interrupted' else completed).append(dict(directory=task['directory'],**e))
    def budget(task):
        directory=pathlib.Path(task['directory'])
        prediction=directory.name=='prediction'
        key=directory.parent.name if prediction else directory.parent.parent.name
        package=task['package']; arm='oracle' if package=='oracle' else key.split('-')[1]
        category='predict' if prediction else ('fit_initial' if directory.name=='01' or package in ('sjSDM','gllvm') else 'fit_extended')
        value=float(policy['budgets'].get(':'.join((package,arm,category)),float('inf')))
        return value if value>0 and math.isfinite(value) else float('inf')
    try:
        while time.time()<cutoff:
            verify_control();update(prefix+'_queue')
            queue=record/'dispatch.json'
            command=[launcher,'-f',script,'--args','--mode=queue','--run='+str(root),'--phase=full','--kind=all','--output='+str(queue)]
            queue_logs=[(record/'queue-stdout.log').open('ab',buffering=0),(record/'queue-stderr.log').open('ab',buffering=0)]
            queue_proc=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=queue_logs[0],stderr=queue_logs[1],start_new_session=True)
            atomic_json(record/'queue-owner.json',dict(pid=queue_proc.pid,command=command,coordinator_pid=os.getpid(),started=time.time()))
            while queue_proc.poll() is None and time.time()<deadline: time.sleep(poll)
            if queue_proc.poll() is None:
                terminate([queue_proc]);queue_interrupted=True
            queue_exit=queue_proc.returncode;terminate([queue_proc]);queue_proc=None
            for log in queue_logs: log.close()
            queue_logs=None
            if queue_interrupted: break
            if queue_exit: raise RuntimeError('Scientific queue failed; inspect window queue logs')
            state=json.loads(queue.read_text());atomic_json(control/'dispatch.json',state)
            if state['status'] not in ('running','complete'): raise RuntimeError('Scientific queue requires attention: '+state['status'])
            if not state['tasks']:
                if state['status']!='complete': raise RuntimeError('Scientific queue requires attention: '+state['status'])
                outcome='complete';break
            pending=list(state['tasks']);ran=False
            update(prefix+'_running',queued=len(pending),workers=workers)
            while pending or active:
                now=time.time()
                if now>=deadline:
                    terminate([v[0] for v in active.values()])
                    for pid in list(active): finish(pid,'resource_interrupted','authorised_window_deadline')
                    break
                for task in list(pending):
                    if len(active)>=workers or time.time()>=cutoff: break
                    if bounded and budget(task)>deadline-time.time(): continue
                    directory=pathlib.Path(task['directory']);request=pathlib.Path(task['request'])
                    if root.joinpath('jobs').resolve() not in directory.resolve().parents or request.parent.resolve()!=directory.resolve():
                        raise RuntimeError('Task path outside owned job archive')
                    if (directory/'owner.json').exists() or (directory/'exit.json').exists(): raise RuntimeError('Task already has process records: '+str(directory))
                    verify_control()
                    logs=[(directory/'stdout.log').open('ab',buffering=0),(directory/'stderr.log').open('ab',buffering=0)]
                    command=[launcher,'-f',script,'--args','--mode=worker','--task='+str(request)]
                    proc=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=logs[0],stderr=logs[1],start_new_session=True)
                    started=time.time();active[proc.pid]=(proc,task,logs,started)
                    atomic_json(directory/'owner.json',dict(pid=proc.pid,command=command,started=started,coordinator_pid=os.getpid(),window_record=str(record)))
                    pending.remove(task);ran=True
                for pid,(proc,task,logs,started) in list(active.items()):
                    if proc.poll() is not None: finish(pid)
                    elif time.time()-started>task['timeout_seconds']:
                        terminate([proc]);finish(pid,'resource_interrupted','per_attempt_limit')
                if interrupted: break
                if not active: break
                time.sleep(poll)
            deferred=[t['directory'] for t in pending]
            if interrupted:
                if not bounded: outcome='needs_attention';error='Per-attempt limit reached; recovery decision required'
                break
            if time.time()>=cutoff or not ran: break
            # Revisit the queue after useful work: shorter predictions/other tasks
            # can fit the remaining window even when a long extension cannot.
    except BaseException as exc:
        error=str(exc);outcome='needs_attention';raise
    finally:
        for s in old_handlers: signal.signal(s,signal.SIG_IGN)
        if queue_proc is not None:
            terminate([queue_proc]);queue_interrupted=True
        if queue_logs:
            for log in queue_logs: log.close()
        if active:
            terminate([v[0] for v in active.values()])
            for pid in list(active): finish(pid,'resource_interrupted','window_controller_stopped')
        result=dict(status=outcome,finished=time.time(),mode=mode,deadline=recorded_deadline,completed=completed,interrupted=interrupted,deferred=deferred,queue_interrupted=queue_interrupted,error=error,owned_processes_remaining=0,requires_authorisation_to_resume=(outcome!='complete'))
        atomic_json(record/'complete.json',result)
        update('paused_by_user' if outcome=='paused_window' else (prefix+'_complete' if outcome=='complete' else prefix+'_needs_attention'),receipt=str(record/'complete.json'),interrupted=len(interrupted))
        print('SPATIAL WINDOW COMPLETE '+outcome,flush=True)
        for s,h in old_handlers.items(): signal.signal(s,h)
        lock.close()
    return 0

if __name__=='__main__':
    raise SystemExit(run_window(*sys.argv[1:]))
