#!/usr/bin/env python3
"""Retire one verified controller between batches; never launch or retry jobs."""
import fcntl, hashlib, json, os, pathlib, signal, subprocess, sys, time
from launch import atomic_json

def read(path):
    return json.loads(path.read_text())

def processes():
    rows=subprocess.check_output(['/bin/ps','-axo','pid=,ppid=,pgid=,stat=,command='],text=True)
    result={}
    for row in rows.splitlines():
        f=row.split(None,4)
        if len(f)==5 and 'Z' not in f[3]:
            result[int(f[0])]=dict(pid=int(f[0]),ppid=int(f[1]),pgid=int(f[2]),stat=f[3],command=f[4])
    return result

def identity_matches(expected):
    result=subprocess.run(['/bin/ps','-p',str(expected['pid']),'-o','lstart=,command='],capture_output=True,text=True)
    fields=result.stdout.strip().split(None,5)
    return result.returncode==0 and len(fields)==6 and ' '.join(fields[:5])==expected['lstart'] and fields[5]==expected['command']

def retire(config_path,record_path):
    config_path=pathlib.Path(config_path).resolve();cfg=read(config_path)
    root=pathlib.Path(cfg['run_root']).resolve();old=pathlib.Path(cfg['coordinator_record']).resolve()
    if root not in old.parents: raise RuntimeError('Controller record is outside run root')
    owner_path=old/'owner.json';owner=read(owner_path);expected=cfg['identity'];pid=expected['pid']
    if owner['pid']!=pid or hashlib.sha256(owner_path.read_bytes()).hexdigest()!=cfg['owner_sha256']:
        raise RuntimeError('Controller ownership changed')
    poll=float(cfg.get('poll_seconds',1));cleanup=float(cfg.get('cleanup_seconds',30))
    queue_wait=float(cfg.get('queue_wait_seconds',1800))
    if not 0<poll<=1 or not 1<=cleanup<=40 or not 1<=queue_wait<=3600: raise RuntimeError('Invalid transition timings')
    lock=(root/'control/retire.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    record=pathlib.Path(record_path).resolve();record.mkdir(parents=True,exist_ok=False)
    atomic_json(record/'owner.json',dict(pid=os.getpid(),started=time.time(),config=str(config_path),config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest()))
    frozen=False;signalled=False
    def checked_signal(sig):
        if not identity_matches(expected): raise RuntimeError('Controller live identity mismatch; refusing signal')
        os.kill(pid,sig)
    def interrupted(signum,_): raise SystemExit(128+signum)
    handlers={s:signal.signal(s,interrupted) for s in (signal.SIGINT,signal.SIGTERM)}
    try:
        while True:
            if (old/'complete.json').exists():
                atomic_json(record/'complete.json',dict(status='predecessor_finished_without_transition',finished=time.time(),predecessor=str(old)))
                return 0
            if not identity_matches(expected): raise RuntimeError('Controller live identity mismatch')
            status=read(root/'continuation-status.json')
            if status.get('window_record')!=str(old): raise RuntimeError('Active controller record changed')
            if status['stage']=='fit_running': time.sleep(poll);continue
            if status['stage']!='fit_queue': raise RuntimeError('Controller is not in a transferable state')
            # A fresh read after SIGSTOP closes the queue-to-worker dispatch race.
            frozen=True;checked_signal(signal.SIGSTOP)
            until=time.monotonic()+2
            while True:
                table=processes()
                if pid in table and 'T' in table[pid]['stat']: break
                if time.monotonic()>=until: raise RuntimeError('Controller stop was not confirmed')
                time.sleep(.02)
            status=read(root/'continuation-status.json')
            if status.get('window_record')!=str(old): raise RuntimeError('Active controller record changed while stopping')
            if status['stage']!='fit_queue':
                checked_signal(signal.SIGCONT);frozen=False;time.sleep(poll);continue
            descendants=set();changed=True
            while changed:
                added={p for p,row in table.items() if row['ppid'] in descendants|{pid}}-descendants
                changed=bool(added);descendants.update(added)
            archive_names={str(root/'jobs'),str(pathlib.Path(cfg['run_root'])/'jobs'),str((root/'jobs').resolve())}
            worker_rows=[row for row in table.values() if '--mode=worker' in row['command'] and any(name in row['command'] for name in archive_names)]
            if worker_rows: raise RuntimeError('A study worker is still live; refusing retirement')
            queue_owner_path=old/'queue-owner.json'
            queue_owner=read(queue_owner_path) if queue_owner_path.exists() else {}
            if queue_owner.get('started',0)<status['updated']:
                checked_signal(signal.SIGCONT);frozen=False;time.sleep(poll);continue
            outputs={'--output='+str(old/'dispatch.json'),'--output='+str(pathlib.Path(cfg['coordinator_record'])/'dispatch.json')}
            if queue_owner.get('coordinator_pid')!=pid or '--mode=queue' not in queue_owner.get('command',[]) or not outputs.intersection(queue_owner.get('command',[])):
                raise RuntimeError('Queue ownership does not match the controller')
            direct={child for child in descendants if table[child]['ppid']==pid}
            if direct-{queue_owner['pid']}:
                raise RuntimeError('Unexpected live descendant; refusing retirement')
            if queue_owner['pid'] in table:
                cmd=table[queue_owner['pid']]['command']
                if '--mode=queue' not in cmd or not any(value in cmd for value in outputs):
                    raise RuntimeError('Queue live command does not match its ownership')
            groups={table[p]['pgid'] for p in descendants|{pid}}
            atomic_json(record/'boundary.json',dict(checked_at=time.time(),status=status,controller=expected,descendants=[table[p] for p in sorted(descendants)],owned_groups=sorted(groups),live_workers=0))
            # Request creation is not a single atomic transaction. Let the
            # administrative queue finish before retiring its stopped parent.
            queue_groups=groups-{table[pid]['pgid']}
            until=time.monotonic()+queue_wait
            while any(row['pgid'] in queue_groups for row in processes().values()):
                if time.monotonic()>=until: raise RuntimeError('Queue drain timed out; resuming original controller')
                time.sleep(poll)
            table=processes()
            if pid not in table or 'T' not in table[pid]['stat'] or read(root/'continuation-status.json')!=status:
                raise RuntimeError('Controller resumed or queue generation changed during drain; refusing retirement')
            if any(row['ppid']==pid or ('--mode=worker' in row['command'] and any(name in row['command'] for name in archive_names)) for row in table.values()):
                raise RuntimeError('A worker or unexpected descendant appeared during drain; refusing retirement')
            queue=old/'dispatch.json'
            if not queue.exists() or queue.stat().st_mtime<queue_owner['started']:
                raise RuntimeError('No fresh completed queue result; resuming original controller')
            state=read(queue)
            if state.get('status') not in ('running','complete'):
                raise RuntimeError('Completed queue requires attention; resuming original controller')
            if not state['tasks']:
                checked_signal(signal.SIGCONT);frozen=False
                atomic_json(record/'complete.json',dict(status='no_successor_needed',finished=time.time(),predecessor=str(old),successor_launched=False))
                return 0
            atomic_json(record/'queue-drained.json',dict(finished=time.time(),queue_owner=queue_owner,queue_sha256=hashlib.sha256(queue.read_bytes()).hexdigest(),tasks=len(state['tasks'])))
            # TERM invokes the old controller's tested cleanup after every
            # worker and the administrative queue have already finished.
            checked_signal(signal.SIGTERM);signalled=True
            checked_signal(signal.SIGCONT);frozen=False
            until=time.monotonic()+cleanup
            while True:
                table=processes();remaining=[row for row in table.values() if row['pgid'] in groups]
                if (old/'complete.json').exists() and not remaining: break
                if time.monotonic()>=until: raise RuntimeError('Controller cleanup not confirmed; successor must not start')
                time.sleep(poll)
            complete=read(old/'complete.json')
            if complete['interrupted'] or complete['owned_processes_remaining']:
                raise RuntimeError('Unexpected interrupted work; successor must not start')
            atomic_json(record/'complete.json',dict(status='retired_at_queue',finished=time.time(),predecessor=str(old),predecessor_complete_sha256=hashlib.sha256((old/'complete.json').read_bytes()).hexdigest(),owned_processes_remaining=0,scientific_attempts_interrupted=0,queue_finished_naturally=True,successor_launched=False))
            return 0
    except BaseException as exc:
        atomic_json(record/'error.json',dict(error=str(exc),finished=time.time(),term_sent=signalled,successor_launched=False))
        raise
    finally:
        for s in handlers: signal.signal(s,signal.SIG_IGN)
        if frozen and identity_matches(expected): os.kill(pid,signal.SIGCONT)
        for s,handler in handlers.items(): signal.signal(s,handler)
        lock.close()

if __name__=='__main__': raise SystemExit(retire(*sys.argv[1:]))
