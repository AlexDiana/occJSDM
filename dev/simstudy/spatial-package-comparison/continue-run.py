#!/usr/bin/env python3
"""Continue the authorised study using only its frozen, gated CLI commands."""
import datetime, hashlib, json, os, pathlib, subprocess, sys, time

def continue_run(root):
    root=pathlib.Path(root).resolve()
    def status(stage, **extra):
        value=dict(stage=stage,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**extra)
        pending=root/'continuation-status.pending.json'
        pending.write_text(json.dumps(value,indent=2)+'\n')
        os.replace(pending,root/'continuation-status.json')
        print(stage,flush=True)
    record=root/(pathlib.Path(__file__).stem+'-record.json')
    if not record.exists(): record=root/'continuation-driver.json'
    if record.exists():
        expected=json.loads(record.read_text())['sha256']
        if hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()!=expected:
            status('driver_identity_mismatch'); return 2
    policy_path=root/'execution-policy.json'
    policy=json.loads(policy_path.read_text()) if policy_path.exists() else None
    def dispatcher_verified():
        files=dict(policy.get('support_files',{}))
        files[str(dispatcher)]=policy['coordinator_sha256']
        return all(pathlib.Path(p).is_file() and hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()==h for p,h in files.items())
    if policy is not None:
        dispatcher=pathlib.Path(policy['coordinator_path'])
        if (not isinstance(policy['workers'],int) or not 1<=policy['workers']<=6 or
            not dispatcher_verified()):
            status('execution_policy_identity_mismatch'); return 2
    status('waiting_for_pilot')
    finished=root/'launch-pilot-01/complete.json'
    while not finished.exists():
        owner=root/'launch-pilot-01/owner.json'
        if owner.exists():
            pid=json.loads(owner.read_text()).get('pid')
            if pid is not None:
                try: os.kill(int(pid),0)
                except ProcessLookupError:
                    status('pilot_needs_attention',reason='Pilot supervisor disappeared without a completion record'); return 2
        time.sleep(30)
    pilot=json.loads(finished.read_text())
    if pilot['exit_status']!=0:
        status('pilot_needs_attention',exit_status=pilot['exit_status']); return 2
    launcher=str(root/'snapshot/study/run-r'); script=str(root/'snapshot/study/run.R')
    def command(mode, *args):
        return [launcher,'-f',script,'--args','--mode='+mode,'--run='+str(root),*args]
    full_command=command('fit','--input-mode=reconstructed','--phase=full','--workers=4')
    if policy is not None:
        full_command=[sys.executable,str(dispatcher),launcher,script,str(root),'full',str(policy['workers']),'all','']
    commands=[('checkpoint',command('checkpoint')),
      ('checkpoint_needs_attention',[launcher,'-e','x <- readRDS('+json.dumps(str(root/'runtime-checkpoint.rds'))+'); quit(status=if(isTRUE(x$passed)) 0L else 21L)']),
      ('full_fit',full_command),
      ('verify_full',command('verify','--phase=full')),
      ('score_full',command('score','--input-mode=reconstructed','--phase=full'))]
    for stage,args in commands:
        if stage=='full_fit' and policy is not None and not dispatcher_verified():
            status('execution_policy_identity_mismatch'); return 2
        status('checking_checkpoint' if stage=='checkpoint_needs_attention' else stage)
        result=subprocess.run(args,stdin=subprocess.DEVNULL)
        if result.returncode:
            status(stage if stage=='checkpoint_needs_attention' else stage+'_needs_attention',exit_status=result.returncode)
            return result.returncode if result.returncode>0 else 2
    status('scientific_results_ready',lesson_export_complete=False)
    return 0

if __name__=='__main__':
    raise SystemExit(continue_run(sys.argv[1]))
