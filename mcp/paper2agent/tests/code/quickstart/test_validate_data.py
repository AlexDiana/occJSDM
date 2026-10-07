import asyncio
import json
import pytest
from fastmcp.exceptions import ToolError
from conftest import ROOT,OUT,invoke,sha

def test_installed_fastmcp_schema_contract():
    from fastmcp import Client
    from tools.quickstart import quickstart_mcp,response
    from typing import get_type_hints
    assert get_type_hints(response)['out'].__name__=='Path'
    async def read():
        async with Client(quickstart_mcp) as c:return await c.list_tools()
    tools={t.name:t.input_schema for t in asyncio.run(read())}
    expected={'validate_data':{'data_path','occCovariates','collCovariates'},'fit_model':{'data_path','options','seed'},'diagnostics':{'fit_id'},'summarise_fit':{'fit_id'}}
    assert set(tools)==set(expected)
    for name,required in expected.items():assert set(tools[name]['required'])==required
    for name in ('diagnostics','summarise_fit'):assert tools[name]['properties']['max_rows']['maximum']==20
    (OUT/'independent_input_schemas.json').write_text(json.dumps(tools,indent=2)+'\n')

@pytest.mark.parametrize('name',['binary_single','occupancy_changed','two_stage_changed'])
def test_fresh_validation_native_design(specs,exercised,name):
    c=specs[name];r=exercised[name]['validate_data']
    assert r['model']==c['model'] and r['dimensions']==c['dimensions']
    assert r['traits_present'] is False and r['design']['species_count']==len(c['species'])
    assert r['issues']==[] and 'OTU' not in r
    assert r['design']['occupancy_columns']==c['options']['occCovariates'] or name!='binary_single'

@pytest.mark.parametrize('fixture,fragment',[('missing_info','info must'),('traits_alias','exact component'),('row_mismatch','row identifiers')])
def test_actual_native_rds_invalid(fixture,fragment):
    r=invoke('validate_data',dict(data_path=str(OUT/(fixture+'.rds')),occCovariates=['habitat_changed'],collCovariates=[]))
    assert r['eligible'] is False and fragment in ' '.join(r['issues'])

def test_source_bindings_immutable():
    bindings=json.loads((ROOT/'reports/tool-bindings.json').read_text())
    assert all(sha(ROOT/'repo/occJSDM'/p)==h for p,h in bindings['source_hashes'].items())

@pytest.mark.parametrize('path',['/etc/passwd','../README.md','notebooks/../README.md'])
def test_path_confinement_is_mcp_error(path):
    with pytest.raises(ToolError):invoke('validate_data',dict(data_path=path,occCovariates=[],collCovariates=[]))
