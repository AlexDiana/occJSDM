import importlib.util, json, pathlib, subprocess, sys, tempfile, unittest
source=pathlib.Path(__file__).resolve().parents[1]/'continue-run.py'
class ContinuationTests(unittest.TestCase):
    def fixture(self, gate=True, pilot=0):
        tmp=tempfile.TemporaryDirectory(); self.addCleanup(tmp.cleanup); root=pathlib.Path(tmp.name)
        (root/'launch-pilot-01').mkdir(); (root/'snapshot/study').mkdir(parents=True)
        (root/'launch-pilot-01/complete.json').write_text(json.dumps({'exit_status':pilot}))
        launcher=root/'snapshot/study/run-r'
        launcher.write_text('#!/usr/bin/env python3\nimport json,pathlib,sys\nr=pathlib.Path(__file__).resolve().parents[2]\nwith (r/"calls.jsonl").open("a") as f: f.write(json.dumps(sys.argv[1:])+"\\n")\nraise SystemExit('+str(0 if gate else 21)+' if "-e" in sys.argv else 0)\n')
        launcher.chmod(0o755)
        return root
    def run_driver(self,root):
        return subprocess.run([sys.executable,str(source),str(root)],capture_output=True,text=True)
    def test_passed_gate_continues_in_order(self):
        root=self.fixture(); p=self.run_driver(root); self.assertEqual(p.returncode,0,p.stderr)
        calls=[json.loads(x) for x in (root/'calls.jsonl').read_text().splitlines()]
        modes=[next((a for a in x if a.startswith('--mode=')),'gate') for x in calls]
        self.assertEqual(modes,['--mode=checkpoint','gate','--mode=fit','--mode=verify','--mode=score'])
        self.assertIn('--workers=4',calls[2]); self.assertIn('--phase=full',calls[2])
        self.assertEqual(json.loads((root/'continuation-status.json').read_text())['stage'],'scientific_results_ready')
    def test_failed_gate_never_dispatches_full_fits(self):
        root=self.fixture(gate=False); p=self.run_driver(root); self.assertEqual(p.returncode,21,p.stderr)
        self.assertEqual(len((root/'calls.jsonl').read_text().splitlines()),2)
        self.assertEqual(json.loads((root/'continuation-status.json').read_text())['stage'],'checkpoint_needs_attention')
    def test_failed_pilot_never_runs_checkpoint_or_more_fits(self):
        root=self.fixture(pilot=2); p=self.run_driver(root); self.assertEqual(p.returncode,2,p.stderr)
        self.assertFalse((root/'calls.jsonl').exists())
        self.assertEqual(json.loads((root/'continuation-status.json').read_text())['stage'],'pilot_needs_attention')
    def test_missing_supervisor_stops_without_dispatch(self):
        root=self.fixture(); (root/'launch-pilot-01/complete.json').unlink()
        (root/'launch-pilot-01/owner.json').write_text(json.dumps({'pid':2147483647}))
        p=subprocess.run([sys.executable,str(source),str(root)],capture_output=True,text=True,timeout=1)
        self.assertEqual(p.returncode,2,p.stderr)
        self.assertFalse((root/'calls.jsonl').exists())
if __name__=='__main__': unittest.main()
