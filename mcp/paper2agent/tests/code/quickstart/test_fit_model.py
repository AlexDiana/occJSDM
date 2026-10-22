import copy
import json
import pytest
from fastmcp.exceptions import ToolError
from conftest import ROOT,OUT,invoke,artifact,sha

def test_full_native_fit_and_installed_source(numerical_audit):
    assert numerical_audit['all_passed'] is True
    assert len(numerical_audit['saved_reference_cases'])==18
    assert len(numerical_audit['changed_cases'])==3
    assert numerical_audit['installed_source_bodies_identical'] is True

@pytest.mark.parametrize('name',['binary_single','occupancy_changed','two_stage_changed'])
def test_provenance_effective_settings_qualifications(specs,exercised,name):
    c=specs[name];r=exercised[name]['fit_model'];m=json.loads(artifact(r,'manifest').read_text())
    assert m['requested_settings']==c['options'] and m['seed']==c['seed']
    assert m['input_sha256']==sha(c['data_path'])==r['input_sha256']
    assert m['effective_settings']['listParams']==c['options']['listParams']
    assert m['runtime']['source_revision']==m['source_revision']=='b7b7001e56cea8e3931b09c917ea0e29e6cef3c6'
    assert m['runtime']['occJSDM_path'].startswith(str(ROOT/'r-runtime'))
    assert r['kept_draws_per_chain']==20 and r['total_iterations_per_chain']==20+20*c['options']['MCMCparams']['nthin']
    assert any('unsuitable for scientific interpretation' in q for q in m['qualifications'])
    for a in r['artifacts']:assert sha(a['path'])==a['sha256']

@pytest.mark.parametrize('update',[{'spatCovariates':['x','y']},{'listPriors':{}},{'threshold':0},{'threshold':1.5},{'listParams':{'n_factors':0,'splineVars':True}}])
def test_closed_options_are_mcp_errors(specs,update):
    c=specs['binary_single'];o=copy.deepcopy(c['options']);o.update(update)
    with pytest.raises(ToolError,match='validation'):
        invoke('fit_model',dict(data_path=c['data_path'],options=o,seed=c['seed']))

def test_native_single_species_two_stage_failure_retained(exercised):
    assert exercised['two_stage_single']['native_error']=='incorrect number of dimensions'

def test_composed_validation_failure_prevents_publication(specs):
    before=set((ROOT/'artifacts/registry').glob('*.json'))
    c=specs['occupancy_changed']
    with pytest.raises(ToolError,match='row identifiers'):
        invoke('fit_model',dict(data_path=str(OUT/'row_mismatch.rds'),options=c['options'],seed=c['seed']))
    assert set((ROOT/'artifacts/registry').glob('*.json'))==before

def test_repeated_fit_new_id_and_exact_native_results(specs,exercised):
    c=specs['binary_single'];r=invoke('fit_model',dict(data_path=c['data_path'],options=c['options'],seed=c['seed']))
    initial=exercised['binary_single']['fit_model']
    assert r['fit_id']!=initial['fit_id']
    assert {a['path'] for a in r['artifacts']}.isdisjoint(a['path'] for a in initial['artifacts'])
    assert sha(artifact(r,'fit'))==sha(artifact(initial,'fit'))
    (OUT/'repeat_fit_result.json').write_text(json.dumps(r,indent=2)+'\n')

def test_unavailable_configured_rscript_is_mcp_error(specs,monkeypatch):
    monkeypatch.setenv('P2A_RSCRIPT','/missing/verifier/Rscript')
    c=specs['binary_single']
    with pytest.raises(ToolError,match='Rscript unavailable'):
        invoke('fit_model',dict(data_path=c['data_path'],options=c['options'],seed=c['seed']))
