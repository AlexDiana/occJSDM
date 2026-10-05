"""Breaks caught: ignored options, native drift, incomplete publication, unsafe paths and unbounded output."""
import asyncio
import csv
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
import pytest
from fastmcp import Client

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
NATIVE = json.loads((ROOT/'notebooks/quickstart/reference/native_execution_results.json').read_text())

def module():
    try:
        return importlib.import_module('tools.quickstart')
    except ModuleNotFoundError:
        pytest.fail('Four-tool wrapper has not been implemented')

def vector(v):
    return v if isinstance(v,list) else [v]

def options(case):
    return dict(occCovariates=vector(case['occCovariates']),collCovariates=vector(case['collCovariates']),threshold=1,listParams={'n_factors':2},MCMCparams=dict(case['settings']))

@pytest.fixture(autouse=True)
def root_env(monkeypatch):
    monkeypatch.setenv('P2A_PROJECT_ROOT',str(ROOT))
    monkeypatch.setenv('P2A_RSCRIPT','/usr/local/bin/Rscript')
    monkeypatch.setenv('P2A_R_PROJECT',str(ROOT/'r-runtime'))

async def call(name, args):
    async with Client(module().quickstart_mcp) as client:
        return (await client.call_tool(name,args)).data

@pytest.mark.parametrize('mutation',[{'spatCovariates':['x','y']},{'threshold':0},{'threshold':1.5},{'listParams':{'n_factors':-1}},{'listParams':{'n_factors':2,'prior':0}},{'MCMCparams':{'nchain':1,'nburn':0,'niter':10,'nthin':1}}])
def test_unknown_and_unsupported_options_fail_before_outputs(mutation):
    o=options(NATIVE['cases']['two_stage']);o.update(mutation)
    with pytest.raises(Exception) as failure:
        asyncio.run(call('fit_model',dict(data_path='notebooks/quickstart/data/two_stage_data.rds',options=o,seed=1702)))
    assert 'validation' in str(failure.value).lower() or 'input' in str(failure.value).lower()

@pytest.mark.parametrize('name',['binary','occupancy','two_stage','no_traits','renamed_single_covariates','permuted_trait_rows','categorical_traits','changed_seed_thinning2','changed_seed_same_schedule','alignment_global_ordered','alignment_global_shuffled','alignment_local_ordered','alignment_local_shuffled','alignment_character_ordered','alignment_character_shuffled','alignment_identifier_covariates','alignment_occupancy_shuffled','boundary_two_stage_NA_OTU'])
def test_native_fit_and_postprocessing_agreement(name):
    # Wrong seed/options, Python replacement math, changed row order or wrong draw orientation fail this comparison.
    c=NATIVE['cases'][name];o=options(c)
    valid=asyncio.run(call('validate_data',dict(data_path=f'notebooks/quickstart/data/{name}_data.rds',occCovariates=o['occCovariates'],collCovariates=o['collCovariates'])))
    assert valid['eligible'] is True
    fit=asyncio.run(call('fit_model',dict(data_path=f'notebooks/quickstart/data/{name}_data.rds',options=o,seed=c['seed'])))
    assert fit['model']==c['inferred_model']
    assert fit['kept_draws_per_chain']==c['settings']['niter']
    fp=next(a['path'] for a in fit['artifacts'] if a['kind']=='fit')
    compare=Path(__file__).with_name('compare_native.R')
    run=subprocess.run(['/usr/local/bin/Rscript','--vanilla',str(compare),str(ROOT),fp,str(ROOT/f'notebooks/quickstart/reference/{name}_fit.rds')],capture_output=True,text=True)
    assert run.returncode==0,run.stdout+run.stderr
    diag=asyncio.run(call('diagnostics',dict(fit_id=fit['fit_id'],max_rows=2)))
    summary=asyncio.run(call('summarise_fit',dict(fit_id=fit['fit_id'],max_rows=2)))
    assert len(diag['rows'])<=2 and len(summary['rows'])<=2
    assert 'baseline' in summary['estimand'].lower()
    for result,kind,ref in [(diag,'diagnostics',f'{name}_diagnostics.csv'),(summary,'summary',f'{name}_occupancy_summary.csv')]:
        actual=Path(next(a['path'] for a in result['artifacts'] if a['kind']==kind)).read_text()
        expected=(ROOT/'notebooks/quickstart/reference'/ref).read_text()
        assert actual==expected
        json.dumps(result,allow_nan=False)
    manifest=json.loads(Path(next(a['path'] for a in fit['artifacts'] if a['kind']=='manifest')).read_text())
    assert manifest['seed']==c['seed'] and manifest['effective_settings']['spatCovariates']==[]
    assert manifest['source_revision']=='b7b7001e56cea8e3931b09c917ea0e29e6cef3c6'
    assert manifest['runtime']['occJSDM_path'].startswith(str(ROOT/'r-runtime'))

@pytest.mark.parametrize('name',['wrong_dimensions','missing_covariate','missing_trait_species','missing_Primer','NA_site_ID','duplicated_species','within_site_covariate_conflict','within_sample_covariate_conflict','infinite_OTU','binary_NA_OTU','occupancy_NA_OTU','binary_unsorted_uniqueSite','occupancy_localSampleIDs','NA_occupancy_covariate','Inf_occupancy_covariate','NaN_occupancy_covariate','no_stage_count_response'])
def test_validation_rejects_native_failures_and_misleading_successes(name):
    c=NATIVE['cases']['boundary_'+name];o=options(c)
    result=asyncio.run(call('validate_data',dict(data_path=f'notebooks/quickstart/data/boundary_{name}_data.rds',occCovariates=o['occCovariates'],collCovariates=o['collCovariates'])))
    assert result['eligible'] is False and result['issues']
    with pytest.raises(Exception):
        asyncio.run(call('fit_model',dict(data_path=f'notebooks/quickstart/data/boundary_{name}_data.rds',options=o,seed=1702)))

@pytest.mark.parametrize('path',['../README.md','/etc/passwd','notebooks/../quickstart/data/two_stage_data.rds'])
def test_input_paths_confined(path):
    with pytest.raises(Exception):
        asyncio.run(call('validate_data',dict(data_path=path,occCovariates=[],collCovariates=[])))

def test_symlink_escape(tmp_path):
    link=ROOT/'tests/code/implementation/escape.rds';link.symlink_to(tmp_path/'outside.rds')
    try:
        with pytest.raises(Exception): asyncio.run(call('validate_data',dict(data_path=str(link),occCovariates=[],collCovariates=[])))
    finally: link.unlink()

@pytest.mark.parametrize('max_rows',[0,21])
def test_rows_cap_and_arbitrary_fit_path(max_rows):
    with pytest.raises(Exception): asyncio.run(call('diagnostics',dict(fit_id='../../reference/two_stage_fit.rds',max_rows=max_rows)))

def test_repeated_calls_fresh_artifacts_and_failed_fit_unpublished():
    c=NATIVE['cases']['binary'];o=options(c)
    a=asyncio.run(call('fit_model',dict(data_path='notebooks/quickstart/data/binary_data.rds',options=o,seed=1702)))
    b=asyncio.run(call('fit_model',dict(data_path='notebooks/quickstart/data/binary_data.rds',options=o,seed=1702)))
    assert a['fit_id']!=b['fit_id']
    assert {x['path'] for x in a['artifacts']}.isdisjoint(x['path'] for x in b['artifacts'])
    for x in a['artifacts']:
        import hashlib
        assert hashlib.sha256(Path(x['path']).read_bytes()).hexdigest()==x['sha256']
    o=options(NATIVE['cases']['two_stage']);o['MCMCparams']['niter']=1
    with pytest.raises(Exception,match='order.max') as failure:asyncio.run(call('fit_model',dict(data_path='notebooks/quickstart/data/two_stage_data.rds',options=o,seed=1702)))
    failed_dir=Path(str(failure.value).split('Local logs: ')[-1])
    assert failed_dir.is_dir() and not (failed_dir/'manifest.json').exists() and not (failed_dir/'fit.rds').exists()

def test_effective_factor_cap_and_pinned_runtime_identity():
    o=options(NATIVE['cases']['binary']);o['listParams']['n_factors']=99
    fit=asyncio.run(call('fit_model',dict(data_path='notebooks/quickstart/data/binary_data.rds',options=o,seed=1702)))
    manifest=json.loads(Path(next(a['path'] for a in fit['artifacts'] if a['kind']=='manifest')).read_text())
    assert manifest['effective_settings']['listParams']['n_factors']==4
    assert manifest['requested_settings']['listParams']['n_factors']==99
    assert manifest['runtime']['source_revision']=='b7b7001e56cea8e3931b09c917ea0e29e6cef3c6'


def register_reference_for_test(path):
    # Test-only seed of a server-owned registry permits native postprocessing fixtures without a new import tool.
    from tools import runtime
    import uuid
    fit_id='fit_'+uuid.uuid4().hex
    out=runtime.fresh_output('fixture')
    fp=out/'fit.rds';import shutil;shutil.copyfile(path,fp)
    mp=out/'manifest.json'
    runtime.write_json_atomic(mp,dict(status='complete',fit_id=fit_id,fit_path=str(fp),fit_sha256=runtime.digest(fp),model='binary',effective_settings={}))
    registry=ROOT/'artifacts/registry';registry.mkdir(parents=True,exist_ok=True)
    runtime.write_json_atomic(registry/(fit_id+'.json'),dict(manifest_path=str(mp),manifest_sha256=runtime.digest(mp)))
    return fit_id,fp


def test_constant_chain_null_and_tampered_fit_rejected():
    module()
    fit_id,fp=register_reference_for_test(ROOT/'notebooks/quickstart/data/diagnostic_constant_species_fit.rds')
    r=asyncio.run(call('diagnostics',dict(fit_id=fit_id,max_rows=20)))
    row=next(x for x in r['rows'] if x['param']=='beta0_psi' and x['idx1']==1)
    assert row['rhat'] is None and r['screening']['unavailable_rhat']>=1
    json.dumps(r,allow_nan=False)
    fp.write_bytes(b'broken fit')
    with pytest.raises(Exception):asyncio.run(call('diagnostics',dict(fit_id=fit_id)))


def fake_rscript(tmp_path):
    script=tmp_path/'fake Rscript'
    script.write_text('#!'+sys.executable+'\nimport json,os,sys,time\nfrom pathlib import Path\nr=json.loads(Path(sys.argv[4]).read_text())\no=Path(r["output_dir"])\n(o/"child.pid").write_text(str(os.getpid()))\n(o/"fit.rds.pending").write_bytes(b"partial")\ntime.sleep(60)\n')
    script.chmod(0o755)
    return script


def test_timeout_reaps_process_and_keeps_pending_fit_unpublished(tmp_path,monkeypatch):
    module();from tools import runtime
    monkeypatch.setenv('P2A_RSCRIPT',str(fake_rscript(tmp_path)))
    monkeypatch.setenv('P2A_TIMEOUT_SECONDS','1')
    out=runtime.fresh_output('fit_model')
    with pytest.raises(RuntimeError,match='exceeded'):
        asyncio.run(runtime.run_r('fit_model',{},out))
    pid=int((out/'child.pid').read_text())
    with pytest.raises(ProcessLookupError):os.kill(pid,0)
    assert not (out/'fit.rds').exists() and not (out/'manifest.json').exists()


def test_cancellation_reaps_native_child(tmp_path,monkeypatch):
    module();from tools import runtime
    monkeypatch.setenv('P2A_RSCRIPT',str(fake_rscript(tmp_path)))
    out=runtime.fresh_output('fit_model')
    async def cancel():
        task=asyncio.create_task(runtime.run_r('fit_model',{},out))
        for _ in range(100):
            if (out/'child.pid').exists():break
            await asyncio.sleep(.02)
        assert (out/'child.pid').exists()
        task.cancel()
        with pytest.raises(asyncio.CancelledError):await task
    asyncio.run(cancel())
    with pytest.raises(ProcessLookupError):os.kill(int((out/'child.pid').read_text()),0)
    assert not (out/'manifest.json').exists()

@pytest.fixture(scope='module')
def extra_inputs():
    import uuid
    folder=ROOT/'artifacts'/('test-inputs-'+uuid.uuid4().hex);folder.mkdir(parents=True,exist_ok=True)
    script=Path(__file__).with_name('extra_inputs.R')
    done=subprocess.run(['/usr/local/bin/Rscript','--vanilla',str(script),str(ROOT),str(folder)],capture_output=True,text=True)
    assert done.returncode==0,done.stdout+done.stderr
    return folder

@pytest.mark.parametrize('name,eligible',[('aligned_rows',True),('default_info_otu_rows',True),('misaligned_rows',False),('malformed',False),('missing_species',False),('nonnumeric_otu',False),('negative_reads',False),('nan_reads',False),('continuous',False),('constant_covariate',False),('traits_alias',False)])
def test_actual_rds_contents_and_provided_row_ids(extra_inputs,name,eligible):
    r=asyncio.run(call('validate_data',dict(data_path=str(extra_inputs/(name+'.rds')),occCovariates=['X_psi.EnvCov.1'],collCovariates=[])))
    assert r['eligible'] is eligible
    if not eligible:assert r['issues']

@pytest.mark.parametrize('n_factors,threshold,nchain,nburn',[(0,2,2,20),(1,3,3,0)])
def test_changed_fit_options_match_direct_native(n_factors,threshold,nchain,nburn):
    o=options(NATIVE['cases']['two_stage']);o['listParams']['n_factors']=n_factors;o['threshold']=threshold
    o['MCMCparams'].update(nchain=nchain,nburn=nburn,niter=10,nthin=2)
    fit=asyncio.run(call('fit_model',dict(data_path='notebooks/quickstart/data/two_stage_data.rds',options=o,seed=1727)))
    fp=Path(next(a['path'] for a in fit['artifacts'] if a['kind']=='fit'))
    settings=fp.with_name('comparison-request.json');settings.write_text(json.dumps(dict(options=o,seed=1727)))
    result=subprocess.run(['/usr/local/bin/Rscript','--vanilla',str(Path(__file__).with_name('compare_changed_options.R')),str(ROOT),str(fp),str(settings)],capture_output=True,text=True)
    assert result.returncode==0,result.stdout+result.stderr
    assert fit['settings']['threshold']==threshold and fit['settings']['listParams']['n_factors']==n_factors
    assert fit['kept_draws_per_chain']==10 and fit['total_iterations_per_chain']==nburn+20

def test_short_fit_qualification_persists_in_downstream_results():
    o=options(NATIVE['cases']['binary'])
    fit=asyncio.run(call('fit_model',dict(data_path='notebooks/quickstart/data/binary_data.rds',options=o,seed=1702)))
    manifest=json.loads(Path(next(a['path'] for a in fit['artifacts'] if a['kind']=='manifest')).read_text())
    assert manifest.get('qualifications')==fit['qualifications']
    for tool in ('diagnostics','summarise_fit'):
        r=asyncio.run(call(tool,dict(fit_id=fit['fit_id'],max_rows=1)))
        assert any('unsuitable for scientific interpretation' in q for q in r.get('fit_qualifications',[]))
