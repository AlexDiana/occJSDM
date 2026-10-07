import json
import pytest
from fastmcp.exceptions import ToolError
from conftest import ROOT,OUT,invoke,artifact

@pytest.mark.parametrize('name',['binary_single','occupancy_changed','two_stage_changed'])
def test_full_direct_diagnostics(exercised,name,numerical_audit):
    r=exercised[name]['diagnostics']
    assert artifact(r,'diagnostics').read_bytes()==(OUT/(name+'_oracle_diagnostics.csv')).read_bytes()
    assert len(r['rows'])==2 and r['screening']['total_rows']==r['total_rows']
    assert any('unsuitable for scientific interpretation' in q for q in r['fit_qualifications'])
    json.dumps(r,allow_nan=False)

@pytest.mark.parametrize('limit',[0,21,1.5])
def test_preview_bounds_valid_fit(exercised,limit):
    with pytest.raises(ToolError,match='validation'):
        invoke('diagnostics',dict(fit_id=exercised['binary_single']['fit_model']['fit_id'],max_rows=limit))

def test_repeat_diagnostic_artifacts_separate(exercised):
    old=exercised['binary_single']['diagnostics'];r=invoke('diagnostics',dict(fit_id=old['fit_id'],max_rows=20))
    assert {a['path'] for a in old['artifacts']}.isdisjoint(a['path'] for a in r['artifacts'])
    assert artifact(old,'diagnostics').read_bytes()==artifact(r,'diagnostics').read_bytes()
    assert len(r['rows'])<=20
    (OUT/'repeat_diagnostics_result.json').write_text(json.dumps(r,indent=2)+'\n')

def test_constant_chain_null_ess_zero_full_native_table(constant_fit):
    fid,fp,mp=constant_fit;r=invoke('diagnostics',dict(fit_id=fid,max_rows=20))
    row=next(x for x in r['rows'] if x['param']=='beta0_psi' and x['idx1']==1)
    assert row['rhat'] is None and row['ess']==0
    assert r['screening']['unavailable_rhat']==1 and r['screening']['unavailable_ess']==0
    assert r['screening']['unavailable_rows']==1
    assert artifact(r,'diagnostics').read_bytes()==(ROOT/'notebooks/quickstart/reference/constant_species_diagnostics.csv').read_bytes()
    (OUT/'constant_diagnostics_result.json').write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')

def test_tampered_native_fit_rejected_before_r(constant_fit):
    fid,fp,mp=constant_fit;data=fp.read_bytes()
    try:
        fp.write_bytes(data+b'corruption')
        with pytest.raises(ToolError,match='Fit artifact integrity'):invoke('diagnostics',dict(fit_id=fid))
    finally:fp.write_bytes(data)
