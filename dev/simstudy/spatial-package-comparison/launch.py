#!/usr/bin/env python3
"""Detached, logged command supervisor. Completion is job-level, not an MCMC checkpoint."""
import datetime, fcntl, hashlib, json, os, pathlib, signal, subprocess, sys, time

def atomic_json(path, value):
    path = pathlib.Path(path)
    temporary = path.with_name(path.name + '.pending-' + str(os.getpid()))
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    os.replace(temporary, path)

def supervise(root, command):
    root = pathlib.Path(root)
    signal.signal(signal.SIGHUP, signal.SIG_IGN)
    lock = (root / 'supervisor.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    started = time.time()
    atomic_json(root / 'owner.json', dict(pid=os.getpid(), command=command, started=started))
    child = None
    def interrupted(signum, _):
        if child is not None and child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
        raise SystemExit(128 + signum)
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    status = 1
    try:
        if sys.platform == 'darwin' and pathlib.Path('/usr/bin/caffeinate').exists():
            command = ['/usr/bin/caffeinate', '-i', '-s'] + command
        with (root / 'stdout.log').open('ab', buffering=0) as stdout, (root / 'stderr.log').open('ab', buffering=0) as stderr:
            child = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=stdout, stderr=stderr, start_new_session=True)
            atomic_json(root / 'child.json', dict(pid=child.pid, command=command, started=time.time()))
            status = child.wait()
    except SystemExit as exc:
        status = int(exc.code)
    finally:
        if child is not None and child.poll() is None:
            try:
                child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
        record = dict(exit_status=status, elapsed_seconds=time.time()-started,
                      finished=datetime.datetime.now(datetime.timezone.utc).isoformat())
        atomic_json(root / 'exit.json', record)
        print('SPATIAL COMMAND COMPLETE exit_status=' + str(status), flush=True)
        record['logs'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in [root/'stdout.log', root/'stderr.log'] if p.exists()}
        atomic_json(root / 'complete.json', record)
    return status

if __name__ == '__main__':
    if len(sys.argv) < 4:
        raise SystemExit('Usage: launch.sh NEW_LOG_DIRECTORY COMMAND [ARG ...]')
    if sys.argv[1] == '--supervise':
        raise SystemExit(supervise(sys.argv[2], sys.argv[3:]))
    root = pathlib.Path(sys.argv[1]).resolve()
    root.mkdir(parents=True, exist_ok=False)
    with (root / 'launcher.log').open('ab', buffering=0) as log:
        process = subprocess.Popen([sys.executable, str(pathlib.Path(__file__).resolve()), '--supervise', str(root)] + sys.argv[2:],
                                   stdin=subprocess.DEVNULL, stdout=log, stderr=log, start_new_session=True, close_fds=True)
    for _ in range(100):
        if (root / 'owner.json').exists():
            print(json.dumps(dict(supervisor_pid=process.pid, logs=str(root))))
            break
        if process.poll() is not None:
            raise SystemExit('Supervisor failed; inspect ' + str(root/'launcher.log'))
        time.sleep(.05)
    else:
        raise SystemExit('Supervisor startup was not confirmed; inspect logs before retrying')
