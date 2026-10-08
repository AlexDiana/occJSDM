import json,pathlib,subprocess,sys,tempfile,unittest
coordinator=pathlib.Path(__file__).resolve().parents[1]/'coordinate.py'
class CoordinatorTests(unittest.TestCase):
    def test_six_owned_workers_are_the_upper_bound(self):
        with tempfile.TemporaryDirectory() as td:
            root=pathlib.Path(td); launcher=root/'dummy-r'; script=root/'unused.R'
            launcher.write_text('''#!/usr/bin/env python3
import json,pathlib,sys,time
args={x[2:].split('=',1)[0]:x.split('=',1)[1] for x in sys.argv if x.startswith('--') and '=' in x}
if args['mode']=='queue':
    root=pathlib.Path(args['run']); done=list((root/'jobs').glob('*/done'))
    tasks=[] if len(done)==8 else [dict(request=str(root/'jobs'/str(i)/'request.json'),directory=str(root/'jobs'/str(i)),timeout_seconds=20) for i in range(8)]
    pathlib.Path(args['output']).write_text(json.dumps(dict(tasks=tasks,status='running' if tasks else 'complete')))
else:
    root=pathlib.Path(args['task']).parent
    (root/'started').write_text(str(time.time()))
    time.sleep(1)
    (root/'done').write_text(str(time.time()))
''')
            launcher.chmod(0o755)
            p=subprocess.run([sys.executable,str(coordinator),str(launcher),str(script),str(root),'full','6','all',''],capture_output=True,text=True,timeout=30)
            self.assertEqual(p.returncode,0,p.stderr)
            jobs=list((root/'jobs').iterdir()); self.assertEqual(len(jobs),8)
            events=[]
            for job in jobs:
                events += [(float((job/'started').read_text()),1),(float((job/'done').read_text()),-1)]
                self.assertEqual(json.loads((job/'exit.json').read_text())['exit_status'],0)
            active=maximum=0
            for _,delta in sorted(events): active+=delta; maximum=max(maximum,active)
            self.assertEqual(maximum,6)
    def test_seven_workers_are_rejected(self):
        p=subprocess.run([sys.executable,str(coordinator),'x','y','z','full','7','all',''],capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0); self.assertIn('1 through 6',p.stderr)
if __name__=='__main__': unittest.main()
