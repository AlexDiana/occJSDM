import importlib.util, json, os, pathlib, signal, subprocess, sys, tempfile, time, unittest
from unittest import mock
SOURCE=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(SOURCE))
spec=importlib.util.spec_from_file_location('repair',SOURCE/'resume-runtime-repair.py')
repair=importlib.util.module_from_spec(spec); spec.loader.exec_module(repair)

class RepairTests(unittest.TestCase):
    def test_draining_preserves_running_workers_and_exit_codes(self):
        with tempfile.TemporaryDirectory() as temp:
            root=pathlib.Path(temp)
            script=root/'coordinator.py'
            script.write_text('''import pathlib, subprocess, sys, json, time, signal
root=pathlib.Path(sys.argv[1]); active=[]
for i in range(2):
 d=root/'jobs'/str(i)/'attempts'/'01'; d.mkdir(parents=True)
 command=[sys.executable,'-c',"import sys,time,pathlib; time.sleep(1); pathlib.Path(sys.argv[1]).with_name('finished').write_text('finished'); sys.exit(int(sys.argv[2]))",str(d/'request.rds'),str(i)]
 p=subprocess.Popen(command,start_new_session=True); active.append((p,d))
 (d/'owner.json').write_text(json.dumps(dict(pid=p.pid,started=time.time())))
def stop(*args):
 for p,d in active:
  result=p.wait(); (d/'exit.json').write_text(json.dumps(dict(exit_status=result,elapsed_seconds=2)))
 sys.exit(0)
signal.signal(signal.SIGTERM,stop)
(root/'ready').write_text(json.dumps([p.pid for p,d in active]))
while True: time.sleep(.1)
''')
            parent=subprocess.Popen([sys.executable,str(script),str(root)])
            try:
                deadline=time.time()+5
                while not (root/'ready').exists() and time.time()<deadline: time.sleep(.02)
                children=json.loads((root/'ready').read_text())
                command=repair.process(parent.pid)['command']
                os.kill(parent.pid,signal.SIGSTOP)
                (root/'runtime-repair-hold.json').write_text(json.dumps(dict(coordinator_pid=parent.pid,workers=children,command=command)))
                repair.drain_held(root,poll_seconds=.05)
                self.assertEqual(parent.wait(timeout=5),0)
                for i in range(2):
                    directory=root/'jobs'/str(i)/'attempts'/'01'
                    self.assertTrue((directory/'finished').exists())
                    self.assertEqual(json.loads((directory/'exit.json').read_text())['exit_status'],i)
                    self.assertTrue((directory/'held-worker-finished.json').exists())
                # Simulate a crash after exits were finalised but before the final
                # drained marker: parent is gone and raw timing archives exist.
                (root/'runtime-repair-drained.json').unlink()
                repair.drain_held(root,poll_seconds=.05)
                self.assertTrue((root/'runtime-repair-drained.json').exists())
            finally:
                if parent.poll() is None:
                    os.kill(parent.pid,signal.SIGCONT); parent.terminate(); parent.wait(timeout=5)

    def pipeline(self,failing_stage):
        with tempfile.TemporaryDirectory() as temp:
            root=pathlib.Path(temp); control=root/'patch'; control.mkdir()
            (control/'manifest.json').write_text('{"files":{}}')
            (root/'runtime-repair-drained.json').write_text('{}')
            calls=[]
            def run(command,**kwargs):
                calls.append(command)
                fail=('--mode=recover-runtime' in command if failing_stage=='recover' else '-e' in command)
                return subprocess.CompletedProcess(command,21 if fail else 0)
            with mock.patch.object(repair.subprocess,'run',side_effect=run):
                self.assertEqual(repair.continue_repaired(root,control),21)
            self.assertFalse(any('full' in call for call in calls))
            self.assertFalse(any('--mode=score' in call for call in calls))
    def test_failed_recovery_blocks_dispatch(self): self.pipeline('recover')
    def test_failed_checkpoint_blocks_full_run(self): self.pipeline('checkpoint')

    def test_timeout_escalates_only_the_owned_worker(self):
        state=dict(parent=12,state='R',command='R request.rds')
        with mock.patch.object(repair,'process',return_value=state), mock.patch.object(repair.os,'killpg') as kill:
            repair.terminate_owned(15,12,'request.rds',grace_seconds=0)
            self.assertEqual(kill.call_args_list,[mock.call(15,signal.SIGTERM),mock.call(15,signal.SIGKILL)])
        with mock.patch.object(repair,'process',return_value=state), mock.patch.object(repair.os,'killpg') as kill:
            with self.assertRaisesRegex(RuntimeError,'identity changed'): repair.terminate_owned(15,13,'request.rds',grace_seconds=0)
            kill.assert_not_called()

    def test_prior_resource_interruption_blocks_resume(self):
        with tempfile.TemporaryDirectory() as temp:
            root=pathlib.Path(temp); control=root/'patch'; control.mkdir()
            (control/'manifest.json').write_text('{"files":{}}')
            (root/'runtime-repair-drained.json').write_text('{"timed_out":[123]}')
            with mock.patch.object(repair.subprocess,'run') as run:
                with self.assertRaisesRegex(RuntimeError,'budget amendment'): repair.continue_repaired(root,control)
                run.assert_not_called()

    def test_changed_patch_blocks_all_execution(self):
        with tempfile.TemporaryDirectory() as temp:
            root=pathlib.Path(temp); control=root/'patch'; control.mkdir()
            (control/'source').write_text('changed')
            (control/'manifest.json').write_text('{"files":{"source":"wrong"}}')
            with mock.patch.object(repair.subprocess,'run') as run:
                with self.assertRaisesRegex(RuntimeError,'hash mismatch'): repair.continue_repaired(root,control)
                run.assert_not_called()

if __name__=='__main__': unittest.main()
