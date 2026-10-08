#!/usr/bin/env python3
"""Bounded fresh-R worker pool. Only this pool's tracked process groups are signalled."""
import fcntl, json, os, pathlib, signal, subprocess, sys, time
from launch import atomic_json

def coordinate(launcher, script, run, phase, workers, kind, jobs):
    root = pathlib.Path(run)
    control = root / 'control'
    control.mkdir(exist_ok=True)
    lock = (control/'dispatch.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    active = {}
    def stop_owned(signum=None, frame=None):
        for proc, task, logs, started in active.values():
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
        if signum is not None:
            raise SystemExit(128+signum)
    signal.signal(signal.SIGTERM, stop_owned)
    signal.signal(signal.SIGINT, stop_owned)
    try:
        while True:
            queue = control/'dispatch.json'
            subprocess.run([launcher,'-f',script,'--args','--mode=queue','--run='+run,
                            '--phase='+phase,'--kind='+kind,'--output='+str(queue)]+(['--jobs='+jobs] if jobs else []),check=True)
            state = json.loads(queue.read_text())
            tasks = state['tasks']
            if not tasks:
                print('SPATIAL DISPATCH COMPLETE '+state['status'],flush=True)
                return 0 if state['status']=='complete' else 2
            pending = list(tasks)
            while pending or active:
                while pending and len(active)<workers:
                    task = pending.pop(0); directory = pathlib.Path(task['directory'])
                    directory.mkdir(parents=True,exist_ok=True)
                    out = (directory/'stdout.log').open('ab',buffering=0)
                    err = (directory/'stderr.log').open('ab',buffering=0)
                    command = [launcher,'-f',script,'--args','--mode=worker','--task='+task['request']]
                    proc = subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=out,stderr=err,start_new_session=True)
                    started = time.time()
                    atomic_json(directory/'owner.json',dict(pid=proc.pid,command=command,started=started,coordinator_pid=os.getpid()))
                    active[proc.pid]=(proc,task,(out,err),started)
                for pid,(proc,task,logs,started) in list(active.items()):
                    elapsed=time.time()-started
                    timeout=elapsed>task['timeout_seconds']
                    if timeout and proc.poll() is None:
                        os.killpg(pid,signal.SIGTERM)
                        try: proc.wait(timeout=10)
                        except subprocess.TimeoutExpired:
                            os.killpg(pid,signal.SIGKILL);proc.wait()
                    if proc.poll() is not None:
                        for log in logs: log.close()
                        atomic_json(pathlib.Path(task['directory'])/'exit.json',dict(exit_status=proc.returncode,
                            elapsed_seconds=elapsed,status='resource_interrupted' if timeout else 'exited'))
                        del active[pid]
                if active: time.sleep(1)
    finally:
        stop_owned()
        for proc,task,logs,started in active.values():
            try: proc.wait(timeout=10)
            except subprocess.TimeoutExpired: os.killpg(proc.pid,signal.SIGKILL);proc.wait()
            for log in logs: log.close()
            atomic_json(pathlib.Path(task['directory'])/'exit.json',dict(exit_status=proc.returncode,
                elapsed_seconds=time.time()-started,status='interrupted'))

if __name__=='__main__':
    launcher,script,run,phase,workers,kind,jobs=sys.argv[1:]
    workers=int(workers)
    if not 1<=workers<=6: raise SystemExit('Workers must be 1 through 6')
    raise SystemExit(coordinate(launcher,script,run,phase,workers,kind,jobs))
