"""Real subprocess checks for bounded dispatch, deadline cleanup and preserved outputs."""
import importlib.util, json, os, pathlib, signal, subprocess, sys, tempfile, time, unittest

SOURCE=pathlib.Path(__file__).resolve().parents[1]/'coordinate-window.py'

class WindowTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root=pathlib.Path(self.tmp.name); (self.root/'control').mkdir(); self.record=self.root/'window'
        self.launcher=self.root/'launcher'; self.launcher.write_text('''#!/usr/bin/env python3
import json,os,pathlib,signal,subprocess,sys,time
args=dict(x[2:].split('=',1) for x in sys.argv if x.startswith('--') and '=' in x)
r=pathlib.Path(__file__).parent
if args['mode']=='queue':
 time.sleep(float((r/'queue-sleep').read_text()) if (r/'queue-sleep').exists() else 0)
 tasks=json.loads((r/'tasks.json').read_text())
 tasks=[x for x in tasks if not (pathlib.Path(x['directory'])/'exit.json').exists()]
 status=(r/'queue-status').read_text() if (r/'queue-status').exists() else ('running' if tasks else 'complete')
 pathlib.Path(args['output']).write_text(json.dumps(dict(status=status,tasks=tasks)))
else:
 p=pathlib.Path(args['task']); x=json.loads(p.read_text())
 if x.get('ignore_term'): signal.signal(signal.SIGTERM,signal.SIG_IGN)
 if (r/'spawn-child').exists():
  child=subprocess.Popen([sys.executable,'-c','import signal,time;signal.signal(signal.SIGTERM,signal.SIG_IGN);time.sleep(20)'])
  (r/'descendant-pid').write_text(str(child.pid))
 (p.parent/'started').write_text(str(time.time()))
 time.sleep(x['sleep'])
 (p.parent/'saved-result').write_text('finished')
'''); self.launcher.chmod(0o755)
        self.policy=dict(run_root=str(self.root),workers=2,dispatch_until=time.time()+3,stop_by=time.time()+4,budgets={'occJSDM:env:fit_initial':.2},grace_seconds=.1,poll_seconds=.02,files_sha256={})
    def task(self,name,sleep=.05,ignore_term=False):
        d=self.root/'jobs'/('occJSDM-env-'+name)/'attempts'/'01';d.mkdir(parents=True)
        p=d/'request.rds';p.write_text(json.dumps(dict(sleep=sleep,ignore_term=ignore_term)))
        return dict(directory=str(d),request=str(p),package='occJSDM',timeout_seconds=24*3600)
    def run_pool(self,tasks):
        self.assertTrue(SOURCE.exists(),'Deadline-aware coordinator is not implemented')
        (self.root/'tasks.json').write_text(json.dumps(tasks));p=self.root/'policy.json';p.write_text(json.dumps(self.policy))
        result=subprocess.run([sys.executable,str(SOURCE),str(self.launcher),'unused.R',str(self.root),str(p),str(self.record)],capture_output=True,text=True,timeout=8)
        self.assertEqual(result.returncode,0,result.stderr)
        return json.loads((self.record/'complete.json').read_text())
    def test_completed_outputs_survive_and_workers_overlap(self):
        a,b=self.task('a',.15),self.task('b',.15)
        out=self.run_pool([a,b]);self.assertEqual(out['status'],'complete')
        self.assertEqual(len(out['completed']),2)
        starts=[float((pathlib.Path(t['directory'])/'started').read_text()) for t in [a,b]]
        self.assertLess(abs(starts[0]-starts[1]),.13)
        for t in [a,b]: self.assertEqual((pathlib.Path(t['directory'])/'saved-result').read_text(),'finished')
    def test_too_long_tasks_are_deferred_without_worker_or_exit(self):
        t=self.task('long');self.policy['budgets']['occJSDM:env:fit_initial']=100
        out=self.run_pool([t]);self.assertEqual(out['status'],'paused_window')
        self.assertFalse((pathlib.Path(t['directory'])/'started').exists())
        self.assertFalse((pathlib.Path(t['directory'])/'exit.json').exists())
    def test_hard_deadline_kills_only_owned_unfinished_worker(self):
        other=subprocess.Popen([sys.executable,'-c','import time;time.sleep(10)']);self.addCleanup(lambda: (other.terminate(),other.wait()))
        t=self.task('hung',10,True);self.policy.update(dispatch_until=time.time()+.5,stop_by=time.time()+.8)
        out=self.run_pool([t]);self.assertEqual(out['status'],'paused_window')
        self.assertEqual(len(out['interrupted']),1);self.assertIsNone(other.poll())
        e=json.loads((pathlib.Path(t['directory'])/'exit.json').read_text());self.assertEqual(e['status'],'resource_interrupted')
        self.assertFalse((pathlib.Path(t['directory'])/'saved-result').exists())
        owner=json.loads((pathlib.Path(t['directory'])/'owner.json').read_text())
        with self.assertRaises(ProcessLookupError): os.kill(owner['pid'],0)
    def test_queue_read_is_also_bounded(self):
        (self.root/'queue-sleep').write_text('10');self.policy.update(dispatch_until=time.time()+.3,stop_by=time.time()+.5)
        out=self.run_pool([]);self.assertEqual(out['status'],'paused_window');self.assertTrue(out['queue_interrupted'])
    def test_no_dispatch_after_cutoff(self):
        t=self.task('never');self.policy['dispatch_until']=time.time()-1
        out=self.run_pool([t]);self.assertEqual(out['status'],'paused_window');self.assertFalse((pathlib.Path(t['directory'])/'started').exists())
    def test_missing_budget_fails_closed(self):
        t=self.task('unknown');self.policy['budgets']={}
        out=self.run_pool([t]);self.assertEqual(out['status'],'paused_window');self.assertFalse((pathlib.Path(t['directory'])/'started').exists())
    def test_descendant_is_killed_when_its_leader_exits_on_term(self):
        (self.root/'spawn-child').touch();t=self.task('tree',10)
        self.policy.update(dispatch_until=time.time()+.5,stop_by=time.time()+.8)
        def cleanup():
            p=self.root/'descendant-pid'
            if p.exists():
                try: os.kill(int(p.read_text()),signal.SIGKILL)
                except ProcessLookupError: pass
        self.addCleanup(cleanup)
        self.run_pool([t]);pid=int((self.root/'descendant-pid').read_text())
        stat=subprocess.run(['/bin/ps','-p',str(pid),'-o','stat='],capture_output=True,text=True)
        self.assertTrue(stat.returncode!=0 or 'Z' in stat.stdout,'Owned descendant escaped cleanup')
    def test_descendant_is_killed_after_natural_leader_exit(self):
        (self.root/'spawn-child').touch();t=self.task('tree',.1)
        def cleanup():
            p=self.root/'descendant-pid'
            if p.exists():
                try: os.kill(int(p.read_text()),signal.SIGKILL)
                except ProcessLookupError: pass
        self.addCleanup(cleanup)
        self.run_pool([t]);pid=int((self.root/'descendant-pid').read_text())
        stat=subprocess.run(['/bin/ps','-p',str(pid),'-o','stat='],capture_output=True,text=True)
        self.assertTrue(stat.returncode!=0 or 'Z' in stat.stdout,'Owned descendant escaped normal completion')
    def test_attention_queue_stops_before_any_other_worker(self):
        self.assertTrue(SOURCE.exists());t=self.task('blocked')
        (self.root/'queue-status').write_text('needs_attention')
        (self.root/'tasks.json').write_text(json.dumps([t]));p=self.root/'policy.json';p.write_text(json.dumps(self.policy))
        result=subprocess.run([sys.executable,str(SOURCE),str(self.launcher),'unused.R',str(self.root),str(p),str(self.record)],capture_output=True,text=True,timeout=8)
        self.assertNotEqual(result.returncode,0)
        receipt=json.loads((self.record/'complete.json').read_text())
        self.assertIn('Scientific queue requires attention',receipt['error'])
        self.assertFalse((pathlib.Path(t['directory'])/'started').exists(),'Runnable task launched despite attention gate')
    def completion_policy(self):
        self.policy.update(mode='run_to_completion',amendment='E15')
        for key in ('dispatch_until','stop_by','budgets'): del self.policy[key]
    def test_completion_mode_runs_without_deadline_or_admission_budget(self):
        self.completion_policy();t=self.task('complete',.2)
        out=self.run_pool([t]);self.assertEqual(out['status'],'complete')
        self.assertIsNone(out['deadline'])
        self.assertFalse(out['requires_authorisation_to_resume'])
        self.assertEqual((pathlib.Path(t['directory'])/'saved-result').read_text(),'finished')
        status=json.loads((self.root/'continuation-status.json').read_text())
        self.assertEqual(status['amendment'],'E15');self.assertEqual(status['stage'],'fit_complete')
    def test_completion_mode_retains_attempt_timeout_and_does_not_retry(self):
        self.completion_policy();t=self.task('limit',10,True);t['timeout_seconds']=.2
        out=self.run_pool([t]);self.assertEqual(out['status'],'needs_attention')
        self.assertEqual(len(out['interrupted']),1);self.assertEqual(out['completed'],[])
        self.assertEqual(out['interrupted'][0]['reason'],'per_attempt_limit')
        self.assertTrue(out['requires_authorisation_to_resume'])
        owner=json.loads((pathlib.Path(t['directory'])/'owner.json').read_text())
        with self.assertRaises(ProcessLookupError): os.kill(owner['pid'],0)
    def test_completion_mode_retains_attention_gate(self):
        self.completion_policy();self.test_attention_queue_stops_before_any_other_worker()
    def test_completion_mode_rejects_conflicting_deadlines(self):
        self.policy.update(mode='run_to_completion')
        (self.root/'tasks.json').write_text('[]');p=self.root/'policy.json';p.write_text(json.dumps(self.policy))
        result=subprocess.run([sys.executable,str(SOURCE),str(self.launcher),'unused.R',str(self.root),str(p),str(self.record)],capture_output=True,text=True,timeout=8)
        self.assertNotEqual(result.returncode,0)
        self.assertFalse((self.record/'owner.json').exists())
    def test_completion_mode_eight_workers_overlap_without_ninth(self):
        self.completion_policy();self.policy['workers']=8
        tasks=[self.task(str(i),.5) for i in range(9)]
        out=self.run_pool(tasks);self.assertEqual(len(out['completed']),9)
        starts=[float((pathlib.Path(t['directory'])/'started').read_text()) for t in tasks]
        self.assertLess(max(starts[:8])-min(starts[:8]),.45)
        self.assertGreater(starts[8]-min(starts[:8]),.45)
    def test_worker_limits_reject_before_starting_any_job(self):
        for mode,workers in [('bounded_window',7),('run_to_completion',9)]:
            with self.subTest(mode=mode,workers=workers):
                self.policy.update(mode=mode,workers=workers)
                if mode=='run_to_completion':
                    for key in ('dispatch_until','stop_by','budgets'): self.policy.pop(key,None)
                p=self.root/'policy.json';p.write_text(json.dumps(self.policy))
                result=subprocess.run([sys.executable,str(SOURCE),str(self.launcher),'unused.R',str(self.root),str(p),str(self.record)],capture_output=True,text=True,timeout=8)
                self.assertNotEqual(result.returncode,0)
                self.assertFalse((self.record/'owner.json').exists())
    def test_changed_control_is_rejected_before_dispatch(self):
        self.assertTrue(SOURCE.exists(),'Deadline-aware coordinator is not implemented')
        t=self.task('tamper');self.policy['files_sha256']={str(self.launcher):'0'*64}
        (self.root/'tasks.json').write_text(json.dumps([t]));p=self.root/'policy.json';p.write_text(json.dumps(self.policy))
        result=subprocess.run([sys.executable,str(SOURCE),str(self.launcher),'unused.R',str(self.root),str(p),str(self.record)],capture_output=True,text=True,timeout=8)
        self.assertNotEqual(result.returncode,0);self.assertFalse((pathlib.Path(t['directory'])/'started').exists())

if __name__=='__main__': unittest.main()
