"""Real process checks for replacing a fixed-concurrency controller between batches."""
import hashlib, importlib.util, json, os, pathlib, signal, subprocess, sys, time, unittest
from unittest import mock

HERE=pathlib.Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('window_fixture',HERE/'test-window.py')
fixture=importlib.util.module_from_spec(spec);spec.loader.exec_module(fixture)
SOURCE=HERE.parent/'retire-at-queue.py'

class RetirementTests(unittest.TestCase):
    setUp=fixture.WindowTests.setUp
    task=fixture.WindowTests.task
    completion_policy=fixture.WindowTests.completion_policy
    def wait_for(self,condition):
        until=time.monotonic()+8
        while not condition():
            self.assertLess(time.monotonic(),until,'Timed out waiting for process evidence')
            time.sleep(.02)
    def start_controller(self,wrapped=False):
        self.completion_policy();self.task_record=self.task('preserved',.6)
        (self.root/'tasks.json').write_text(json.dumps([self.task_record]))
        p=self.root/'policy.json';p.write_text(json.dumps(self.policy))
        self.logs=(self.root/'controller.log').open('w');self.addCleanup(self.logs.close)
        command=[sys.executable,str(fixture.SOURCE),str(self.launcher),'unused.R',str(self.root),str(p),str(self.record)]
        if wrapped:
            self.supervisor=self.root/'supervisor';self.supervisor.mkdir()
            command=[sys.executable,str(SOURCE.parent/'launch.py'),'--supervise',str(self.supervisor)]+command
        self.proc=subprocess.Popen(command,stdout=self.logs,stderr=self.logs,start_new_session=True)
        def cleanup():
            if self.proc.poll() is None: self.proc.terminate();self.proc.wait(timeout=4)
        self.addCleanup(cleanup)
        self.wait_for(lambda:(pathlib.Path(self.task_record['directory'])/'started').exists())
        # Keep the next scientific queue alive long enough to test retirement.
        (self.root/'queue-sleep').write_text('1.2')
        self.next_task=self.task('next-batch',.2)
        (self.root/'tasks.json').write_text(json.dumps([self.task_record,self.next_task]))
    def config(self,stale=False):
        controller_pid=json.loads((self.record/'owner.json').read_text())['pid']
        line=subprocess.check_output(['/bin/ps','-p',str(controller_pid),'-o','lstart=,command='],text=True).strip()
        fields=line.split(None,5)
        value=dict(run_root=str(self.root),coordinator_record=str(self.record),owner_sha256=hashlib.sha256((self.record/'owner.json').read_bytes()).hexdigest(),identity=dict(pid=controller_pid,lstart=' '.join(fields[:5]),command=fields[5]),poll_seconds=.02,cleanup_seconds=4)
        if hasattr(self,'supervisor'):
            child_path=self.supervisor/'child.json';child=json.loads(child_path.read_text())
            rows=subprocess.check_output(['/bin/ps','-axo','pid=,ppid=,pgid=,lstart=,command='],text=True)
            helpers=[row.split(None,8) for row in rows.splitlines() if row.split(None,2)[1]==str(controller_pid)]
            helpers=[f for f in helpers if f[8]==' '.join(child['command'])]
            self.assertEqual(len(helpers),1,helpers);f=helpers[0]
            value['caffeinate']=dict(pid=int(f[0]),ppid=int(f[1]),pgid=int(f[2]),lstart=' '.join(f[3:8]),command=f[8],supervisor_child=str(child_path),supervisor_child_sha256=hashlib.sha256(child_path.read_bytes()).hexdigest())
        if stale: value['identity']['lstart']='stale process identity'
        path=self.root/'transition.json';path.write_text(json.dumps(value));return path
    def retire(self,path):
        self.assertTrue(SOURCE.exists(),'Queue-boundary retirement has not been implemented')
        return subprocess.run([sys.executable,str(SOURCE),str(path),str(self.root/'retirement')],capture_output=True,text=True,timeout=10)
    def test_waits_for_worker_then_retires_queue_without_interrupting_fit(self):
        self.start_controller();result=self.retire(self.config())
        self.assertEqual(result.returncode,0,result.stderr)
        self.proc.wait(timeout=3)
        task=pathlib.Path(self.task_record['directory'])
        self.assertEqual((task/'saved-result').read_text(),'finished')
        self.assertEqual(json.loads((task/'exit.json').read_text())['exit_status'],0)
        old=json.loads((self.record/'complete.json').read_text())
        self.assertEqual(old['interrupted'],[])
        receipt=json.loads((self.root/'retirement/complete.json').read_text())
        self.assertEqual(receipt['status'],'retired_at_queue')
        self.assertEqual(receipt['owned_processes_remaining'],0)
        queue=json.loads((self.record/'dispatch.json').read_text())
        self.assertEqual(queue['tasks'],[self.next_task],'Administrative queue was interrupted before its result was saved')
        self.assertFalse((pathlib.Path(self.next_task['directory'])/'started').exists())
        self.assertFalse((self.root/'successor').exists(),'Retirement must never launch a successor')
    @unittest.skipUnless(sys.platform=='darwin','macOS caffeinate topology')
    def test_real_supervisor_caffeinate_exits_after_natural_queue_retirement(self):
        self.start_controller(wrapped=True);path=self.config();cfg=json.loads(path.read_text())
        result=self.retire(path);self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(self.proc.wait(timeout=3),143)
        receipt=json.loads((self.root/'retirement/complete.json').read_text())
        self.assertTrue(receipt['queue_finished_naturally'])
        self.assertEqual(json.loads((self.record/'complete.json').read_text())['interrupted'],[])
        self.assertFalse((pathlib.Path(self.next_task['directory'])/'started').exists())
        remaining=subprocess.run(['/bin/ps','-p',str(cfg['caffeinate']['pid']),'-o','stat='],capture_output=True,text=True)
        self.assertTrue(remaining.returncode!=0 or 'Z' in remaining.stdout)
    @unittest.skipUnless(sys.platform=='darwin','macOS caffeinate topology')
    def test_caffeinate_stale_identity_refuses_retirement_and_resumes_controller(self):
        self.start_controller(wrapped=True);path=self.config();cfg=json.loads(path.read_text())
        cfg['caffeinate']['lstart']='stale helper identity';path.write_text(json.dumps(cfg))
        result=self.retire(path);self.assertNotEqual(result.returncode,0)
        self.assertIn('Caffeinate live identity mismatch',result.stderr)
        state=subprocess.check_output(['/bin/ps','-p',str(cfg['identity']['pid']),'-o','stat='],text=True)
        self.assertNotIn('T',state)
        self.assertFalse((self.record/'complete.json').exists())
    @unittest.skipUnless(sys.platform=='darwin','macOS caffeinate topology')
    def test_caffeinate_requires_bound_supervisor_child_record(self):
        self.start_controller(wrapped=True);path=self.config();cfg=json.loads(path.read_text())
        cfg['caffeinate']['supervisor_child_sha256']='0'*64;path.write_text(json.dumps(cfg))
        result=self.retire(path);self.assertNotEqual(result.returncode,0)
        self.assertIn('Caffeinate ownership mismatch',result.stderr)
        self.assertIsNone(self.proc.poll())
    def test_live_study_worker_guard_resumes_controller_without_signalling_worker(self):
        self.start_controller()
        self.wait_for(lambda:json.loads((self.root/'continuation-status.json').read_text())['stage']=='fit_queue')
        stranger=subprocess.Popen([sys.executable,'-c','import time;time.sleep(10)','--mode=worker',str(self.root/'jobs/stray')])
        self.addCleanup(lambda:(stranger.terminate(),stranger.wait()))
        result=self.retire(self.config())
        self.assertNotEqual(result.returncode,0)
        self.assertIn('worker is still live',result.stderr)
        self.assertIsNone(stranger.poll());self.assertIsNone(self.proc.poll())
        status=subprocess.check_output(['/bin/ps','-p',str(self.proc.pid),'-o','stat='],text=True)
        self.assertNotIn('T',status,'Refused transition left old controller frozen')
    def test_wrong_live_identity_refuses_to_signal_controller(self):
        self.start_controller();result=self.retire(self.config(stale=True))
        self.assertNotEqual(result.returncode,0)
        self.assertIsNone(self.proc.poll())
        self.wait_for(lambda:(pathlib.Path(self.task_record['directory'])/'exit.json').exists())
        self.assertEqual((pathlib.Path(self.task_record['directory'])/'saved-result').read_text(),'finished')
    def test_interruption_immediately_after_stop_resumes_old_controller(self):
        self.start_controller()
        self.wait_for(lambda:json.loads((self.root/'continuation-status.json').read_text())['stage']=='fit_queue')
        sys.path.insert(0,str(SOURCE.parent))
        spec=importlib.util.spec_from_file_location('retire_fixture',SOURCE)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        original_kill=os.kill
        def interrupt_after_stop(pid,sig):
            original_kill(pid,sig)
            if sig==signal.SIGSTOP: raise SystemExit(143)
        try:
            with mock.patch.object(module.os,'kill',side_effect=interrupt_after_stop):
                with self.assertRaises(SystemExit): module.retire(self.config(),self.root/'interrupted-retirement')
            status=subprocess.check_output(['/bin/ps','-p',str(self.proc.pid),'-o','stat='],text=True)
            self.assertNotIn('T',status,'Interrupted helper left the controller frozen')
        finally:
            original_kill(self.proc.pid,signal.SIGCONT)
    def test_queue_owner_must_be_published_before_retiring_controller(self):
        self.start_controller()
        self.wait_for(lambda:json.loads((self.root/'continuation-status.json').read_text())['stage']=='fit_queue')
        worker_end=json.loads((pathlib.Path(self.task_record['directory'])/'owner.json').read_text())['started']
        owner_path=self.record/'queue-owner.json'
        self.wait_for(lambda:json.loads(owner_path.read_text())['started']>worker_end)
        saved=owner_path.read_text();owner_path.unlink()
        helper=subprocess.Popen([sys.executable,str(SOURCE),str(self.config()),str(self.root/'retirement')],stdout=self.logs,stderr=self.logs)
        try:
            self.wait_for(lambda:(self.root/'retirement/owner.json').exists())
            time.sleep(.2)
            self.assertIsNone(helper.poll(),'Helper did not wait for queue ownership publication')
            self.assertIsNone(self.proc.poll())
            owner_path.write_text(saved)
            self.assertEqual(helper.wait(timeout=5),0)
        finally:
            if helper.poll() is None: helper.terminate();helper.wait(timeout=3)
    def test_queue_drain_timeout_resumes_original_controller(self):
        self.start_controller();(self.root/'queue-sleep').write_text('3')
        path=self.config();cfg=json.loads(path.read_text());cfg['queue_wait_seconds']=1;path.write_text(json.dumps(cfg))
        result=self.retire(path)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('Queue drain timed out',result.stderr)
        self.assertIsNone(self.proc.poll())
        status=subprocess.check_output(['/bin/ps','-p',str(self.proc.pid),'-o','stat='],text=True)
        self.assertNotIn('T',status)
        self.assertFalse((self.record/'complete.json').exists())
    def test_finished_study_resumes_without_replacement(self):
        self.start_controller();(self.root/'tasks.json').write_text('[]')
        result=self.retire(self.config());self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(self.proc.wait(timeout=4),0)
        self.assertEqual(json.loads((self.record/'complete.json').read_text())['status'],'complete')
        self.assertEqual(json.loads((self.root/'retirement/complete.json').read_text())['status'],'no_successor_needed')
    def test_external_resume_during_queue_drain_cannot_interrupt_new_workers(self):
        self.start_controller()
        # Make the second batch last long enough to observe whether TERM hurt it.
        pathlib.Path(self.next_task['request']).write_text(json.dumps(dict(sleep=1.5,ignore_term=False)))
        helper=subprocess.Popen([sys.executable,str(SOURCE),str(self.config()),str(self.root/'retirement')],stdout=self.logs,stderr=self.logs)
        try:
            self.wait_for(lambda:(self.root/'retirement/boundary.json').exists())
            os.kill(self.proc.pid,signal.SIGCONT)
            self.assertNotEqual(helper.wait(timeout=6),0)
            self.wait_for(lambda:(pathlib.Path(self.next_task['directory'])/'exit.json').exists())
            exit_record=json.loads((pathlib.Path(self.next_task['directory'])/'exit.json').read_text())
            self.assertEqual(exit_record['exit_status'],0,'External resume caused interruption of next-batch worker')
            self.assertEqual((pathlib.Path(self.next_task['directory'])/'saved-result').read_text(),'finished')
        finally:
            if helper.poll() is None: helper.terminate();helper.wait(timeout=3)

if __name__=='__main__':unittest.main()
