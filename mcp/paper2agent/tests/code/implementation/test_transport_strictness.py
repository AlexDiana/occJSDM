"""Catch FastMCP transport overriding strict fields and silently coercing numeric inputs."""
import asyncio
import copy
import json
import os
from pathlib import Path
import shutil
import sys
import uuid

from fastmcp import Client
from fastmcp.client.transports import StdioTransport

ROOT = Path(__file__).resolve().parents[3]
OPTIONS = dict(occCovariates=['X_psi.EnvCov.1','X_psi.EnvCov.2'],collCovariates=[],threshold=1,
               listParams={'n_factors':2},MCMCparams={'nchain':2,'nburn':20,'niter':20,'nthin':1})


def stdio_client(folder):
    env = dict(os.environ,P2A_PROJECT_ROOT=str(folder),P2A_R_PROJECT=str(ROOT/'r-runtime'),P2A_RSCRIPT='/usr/local/bin/Rscript')
    return Client(StdioTransport(command=sys.executable,args=[str(ROOT/'src/occJSDM_mcp.py')],env=env,cwd=str(ROOT),keep_alive=False,log_file=folder/'server-stderr.log'))


def project(label):
    folder=ROOT/'reports/cloud-transfer'/f'{label}-{uuid.uuid4().hex}'
    folder.mkdir(parents=True)
    shutil.copyfile(ROOT/'notebooks/quickstart/data/binary_data.rds',folder/'binary.rds')
    return folder


def published(folder):
    return set(folder.glob('artifacts/registry/*.json')),set(folder.glob('artifacts/calls/*/manifest.json'))


def test_stdio_rejects_boolean_seed_without_publishing():
    # Without server-level strictness, True is accepted as1 and produces a scientifically different fit.
    folder=project('boolean-seed')
    async def scenario():
        async with stdio_client(folder) as client:
            result=await client.call_tool('fit_model',dict(data_path='binary.rds',options=OPTIONS,seed=True),raise_on_error=False)
            assert result.is_error, f'Boolean seed was accepted and published: {result.data}'
            text=' '.join(x.text for x in result.content if hasattr(x,'text')).lower()
            assert 'validation' in text and 'int' in text
            assert published(folder)==(set(),set())
    asyncio.run(scenario())


def test_stdio_strict_scalar_nested_types_and_valid_composition():
    folder=project('strict-types')
    async def scenario():
        async with stdio_client(folder) as client:
            inventory=await client.list_tools()
            assert {t.name for t in inventory}=={'validate_data','fit_model','diagnostics','summarise_fit'}
            schema={t.name:t.input_schema for t in inventory}
            assert schema['diagnostics']['properties']['max_rows']['minimum']==1
            assert schema['diagnostics']['properties']['max_rows']['maximum']==20
            fit=(await client.call_tool('fit_model',dict(data_path='binary.rds',options=OPTIONS,seed=1702))).data
            assert fit['model']=='binary' and fit['kept_draws_per_chain']==20
            manifest=json.loads(Path(next(a['path'] for a in fit['artifacts'] if a['kind']=='manifest')).read_text())
            assert manifest['seed']==1702 and type(manifest['seed']) is int
            before=published(folder)
            invalid=[]
            for seed in (True,1.0,'1702'):
                invalid.append(('fit_model',dict(data_path='binary.rds',options=OPTIONS,seed=seed)))
            for location,key,value in [('root','threshold',True),('root','threshold',1.0),('root','threshold','1'),('listParams','n_factors',False),('listParams','n_factors',2.0),('MCMCparams','nchain',2.0),('MCMCparams','nburn',False),('MCMCparams','niter','20'),('MCMCparams','nthin',1.0)]:
                changed=copy.deepcopy(OPTIONS)
                (changed if location=='root' else changed[location])[key]=value
                invalid.append(('fit_model',dict(data_path='binary.rds',options=changed,seed=1702)))
            for tool in ('diagnostics','summarise_fit'):
                for rows in (True,1.0,'1'):
                    invalid.append((tool,dict(fit_id=fit['fit_id'],max_rows=rows)))
            invalid.append(('validate_data',dict(data_path='binary.rds',occCovariates='X_psi.EnvCov.1',collCovariates=[])))
            invalid.append(('validate_data',dict(data_path='binary.rds',occCovariates=[1],collCovariates=[])))
            for tool,args in invalid:
                result=await client.call_tool(tool,args,raise_on_error=False)
                assert result.is_error, f'Coercible payload unexpectedly accepted: {tool} {args}'
                text=' '.join(x.text for x in result.content if hasattr(x,'text')).lower()
                assert 'validation' in text, f'Error came from tool body instead of strict transport validation: {text}'
                assert published(folder)==before
            for tool in ('diagnostics','summarise_fit'):
                result=(await client.call_tool(tool,dict(fit_id=fit['fit_id'],max_rows=1))).data
                assert len(result['rows'])==1
                assert any('unsuitable for scientific interpretation' in q for q in result['fit_qualifications'])
    asyncio.run(scenario())


def test_module_server_strict_validation_precedes_file_access(monkeypatch):
    # Both isolated module use and top-level mounted stdio use must enforce the same contract.
    sys.path.insert(0,str(ROOT/'src'))
    from tools.quickstart import quickstart_mcp
    monkeypatch.setenv('P2A_PROJECT_ROOT',str(ROOT))
    async def scenario():
        async with Client(quickstart_mcp) as client:
            result=await client.call_tool('fit_model',dict(data_path='missing-unreadable.rds',options=OPTIONS,seed=True),raise_on_error=False)
            text=' '.join(x.text for x in result.content if hasattr(x,'text')).lower()
            assert result.is_error and 'validation' in text and 'int' in text
    asyncio.run(scenario())
