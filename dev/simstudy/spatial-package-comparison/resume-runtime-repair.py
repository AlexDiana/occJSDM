#!/usr/bin/env python3
"""Drain the held pilot, recover validation failures, then resume its gated run."""
import hashlib, json, os, pathlib, signal, subprocess, sys, time
from launch import atomic_json

def process(pid):
    result = subprocess.run(['/bin/ps','-p',str(pid),'-o','pid=,ppid=,stat=,command='],capture_output=True,text=True)
    if result.returncode: return None
    parts = result.stdout.strip().split(None,3)
    if len(parts)!=4: raise RuntimeError('Cannot inspect owned process '+str(pid))
    return dict(pid=int(parts[0]),parent=int(parts[1]),state=parts[2],command=parts[3])

def terminate_owned(pid, coordinator, request, grace_seconds=10):
    def running():
        state=process(pid)
        if state is None or 'Z' in state['state']: return False
        if state['parent']!=coordinator or str(request) not in state['command']:
            raise RuntimeError('Worker identity changed during timeout handling')
        return True
    if running(): os.killpg(pid,signal.SIGTERM)
    deadline=time.time()+grace_seconds
    while running() and time.time()<deadline: time.sleep(.1)
    if running(): os.killpg(pid,signal.SIGKILL)

def drain_held(root, poll_seconds=5):
    root=pathlib.Path(root)
    hold=json.loads((root/'runtime-repair-hold.json').read_text())
    coordinator=int(hold['coordinator_pid'])
    phase_path=root/'runtime-repair-draining.json'
    current=process(coordinator)
    if not phase_path.exists() and (current is None or current['command'].strip()!=hold['command'].strip() or 'T' not in current['state']):
        raise RuntimeError('Held coordinator identity/state changed')
    owners={}
    for path in (root/'jobs').glob('*/attempts/*/owner.json'):
        owner=json.loads(path.read_text())
        if owner['pid'] in hold['workers']:
            if owner['pid'] in owners: raise RuntimeError('Duplicate worker owner')
            owners[owner['pid']]=(path.parent,owner)
    if set(owners)!=set(hold['workers']): raise RuntimeError('Missing worker owner')
    owner_hashes={str(directory/'owner.json'):hashlib.md5((directory/'owner.json').read_bytes()).hexdigest() for directory,_ in owners.values()}
    hold_hash=hashlib.md5((root/'runtime-repair-hold.json').read_bytes()).hexdigest()
    timeout_path=root/'runtime-repair-timeouts.json'
    timed_out=set(json.loads(timeout_path.read_text())['workers']) if timeout_path.exists() else set()
    if phase_path.exists():
        phase=json.loads(phase_path.read_text())
        if phase['owner_hashes']!=owner_hashes or phase['hold_hash']!=hold_hash:
            raise RuntimeError('Drain receipt identity changed')
        timed_out=set(phase['timed_out']); pending=set()
    else: pending=set(owners)
    while pending:
        for pid in list(pending):
            directory,owner=owners[pid]
            live=process(pid)
            if live is not None and live['parent']!=coordinator:
                raise RuntimeError('Worker ownership changed')
            finished=live is None or 'Z' in live['state']
            if not finished:
                if str(directory/'request.rds') not in live['command']: raise RuntimeError('Worker command changed')
                if time.time()-owner['started']>24*3600:
                    timed_out.add(pid)
                    atomic_json(timeout_path,dict(workers=sorted(timed_out)))
                    terminate_owned(pid,coordinator,directory/'request.rds')
                continue
            observation=directory/'held-worker-finished.json'
            if not observation.exists():
                atomic_json(observation,dict(pid=pid,observed_finished=time.time(),
                    elapsed_seconds_upper_bound=time.time()-owner['started'],poll_seconds=poll_seconds,
                    note='Conservative upper bound; first inspection may occur after completion'))
            pending.remove(pid)
        if pending: time.sleep(poll_seconds)
    # Every worker is finished. The original coordinator can now reap its children
    # and write the real exit codes without interrupting any fit.
    if not phase_path.exists():
        atomic_json(phase_path,dict(owner_hashes=owner_hashes,hold_hash=hold_hash,timed_out=sorted(timed_out)))
    current=process(coordinator)
    if current is not None and 'Z' not in current['state']:
        if current['command'].strip()!=hold['command'].strip(): raise RuntimeError('Coordinator changed before drain')
        for pid in owners:
            worker=process(pid)
            if worker is not None and 'Z' not in worker['state']: raise RuntimeError('Recorded finished worker is unexpectedly active')
        os.kill(coordinator,signal.SIGTERM)
        os.kill(coordinator,signal.SIGCONT)
    deadline=time.time()+30
    while time.time()<deadline:
        current=process(coordinator)
        if current is None or 'Z' in current['state']: break
        time.sleep(.2)
    else: raise RuntimeError('Original coordinator did not finish draining')
    for pid,(directory,_) in owners.items():
        if not (directory/'exit.json').exists(): raise RuntimeError('Missing drained worker exit record')
        # The stopped coordinator reaps a finished child late. Preserve its raw
        # record, then use the independently observed finish-time upper bound so
        # idle time waiting for the other fits is not treated as computing time.
        exit_path=directory/'exit.json'
        original=directory/'held-coordinator-exit.json'
        if not original.exists():
            pending_original=original.with_suffix('.pending')
            pending_original.write_bytes(exit_path.read_bytes()); os.replace(pending_original,original)
        record=json.loads(original.read_text())
        observation=json.loads((directory/'held-worker-finished.json').read_text())
        record['elapsed_seconds']=observation['elapsed_seconds_upper_bound']
        record['timing_method']='observed_finish_upper_bound_during_coordinator_hold'
        record['raw_coordinator_elapsed_seconds']=json.loads(original.read_text())['elapsed_seconds']
        if pid in timed_out: record['status']='resource_interrupted'
        if exit_path.read_bytes()!=original.read_bytes() and json.loads(exit_path.read_text())!=record:
            raise RuntimeError('Drain timing record changed')
        atomic_json(exit_path,record)
    atomic_json(root/'runtime-repair-drained.json',dict(coordinator_pid=coordinator,workers=hold['workers'],finished=time.time(),timed_out=sorted(timed_out)))
    if timed_out: raise RuntimeError('Resource-interrupted attempts require a budget amendment; all held workers safely drained')

def continue_repaired(root, control):
    root=pathlib.Path(root).resolve(); control=pathlib.Path(control).resolve()
    def verify_control():
        record=json.loads((control/'manifest.json').read_text())
        if any(hashlib.md5((control/name).read_bytes()).hexdigest()!=digest for name,digest in record['files'].items()):
            raise RuntimeError('Runtime repair source hash mismatch')
    def status(stage,**extra):
        atomic_json(root/'continuation-status.json',dict(stage=stage,updated=time.time(),amendment='E12',**extra))
        print(stage,flush=True)
    verify_control()
    status('waiting_for_active_pilot_fits')
    if not (root/'runtime-repair-drained.json').exists(): drain_held(root)
    if json.loads((root/'runtime-repair-drained.json').read_text()).get('timed_out'):
        raise RuntimeError('Resource-interrupted attempts require a budget amendment')
    launcher=str(root/'snapshot/study/run-r'); script=str(control/'run.R')
    def cli(mode,*args): return [launcher,'-f',script,'--args','--mode='+mode,'--run='+str(root),*args]
    def pool(phase,workers): return [sys.executable,str(control/'coordinate.py'),launcher,script,str(root),phase,str(workers),'all','']
    checkpoint_code='x <- readRDS('+json.dumps(str(root/'runtime-checkpoint.rds'))+'); stopifnot(identical(x$snapshot_hash,unname(tools::md5sum('+json.dumps(str(root/'snapshot-record.rds'))+'))), identical(x$limit_hours,readRDS('+json.dumps(str(root/'config.rds'))+')$remaining_wall_hours)); quit(status=if(isTRUE(x$passed)) 0L else 21L)'
    commands=[('recover_hmsc',cli('recover-runtime')),
        ('finish_pilot',pool('pilot',4)),('checkpoint',cli('checkpoint')),
        ('checking_checkpoint',[launcher,'-e',checkpoint_code]),
        ('full_fit',pool('full',6)),('verify_full',cli('verify','--phase=full')),
        ('score_full',cli('score','--input-mode=reconstructed','--phase=full'))]
    for stage,command in commands:
        if stage=='checkpoint' and (root/'runtime-checkpoint.rds').exists(): continue
        verify_control(); status(stage)
        result=subprocess.run(command,stdin=subprocess.DEVNULL)
        if result.returncode:
            status(stage+'_needs_attention',exit_status=result.returncode)
            return result.returncode if result.returncode>0 else 2
    status('scientific_results_ready',lesson_export_complete=False)
    return 0

if __name__=='__main__':
    try:
        raise SystemExit(continue_repaired(*sys.argv[1:]))
    except Exception as exc:
        root=pathlib.Path(sys.argv[1])
        atomic_json(root/'continuation-status.json',dict(stage='runtime_repair_needs_attention',error=str(exc),updated=time.time()))
        raise
