import json
import pytest
from fastmcp.exceptions import ToolError
from conftest import OUT,invoke,artifact

@pytest.mark.parametrize('name',['binary_single','occupancy_changed','two_stage_changed'])
def test_full_baseline_summary_labels_and_single_species(exercised,name,numerical_audit):
    r=exercised[name]['summarise_fit']
    assert artifact(r,'summary').read_bytes()==(OUT/(name+'_oracle_summary.csv')).read_bytes()
    assert len(r['rows'])==1 and r['rows'][0]['species']=='verifier_species_1'
    assert r['total_rows']==(1 if name=='binary_single' else 3)
    assert 'inverse-logit species intercept' in r['estimand']
    assert any('unsuitable for scientific interpretation' in q for q in r['fit_qualifications'])
    assert any('not occupancy at a particular site' in q for q in r['qualifications'])

def test_repeated_summaries_and_unknown_identifier(exercised):
    old=exercised['binary_single']['summarise_fit'];r=invoke('summarise_fit',dict(fit_id=old['fit_id'],max_rows=20))
    assert {a['path'] for a in old['artifacts']}.isdisjoint(a['path'] for a in r['artifacts'])
    assert artifact(old,'summary').read_bytes()==artifact(r,'summary').read_bytes()
    with pytest.raises(ToolError,match='Unknown fit_id'):invoke('summarise_fit',dict(fit_id='../fit.rds'))
    (OUT/'repeat_summary_result.json').write_text(json.dumps(r,indent=2)+'\n')
