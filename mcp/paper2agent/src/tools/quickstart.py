"""Four native R tools derived from the pinned occJSDM quickstart."""
import os
import shutil
import sys
from importlib.metadata import version
from pathlib import Path
from typing import Annotated
import uuid
from fastmcp import FastMCP
from pydantic import BaseModel, ConfigDict, Field
from .runtime import (artifact, confined_path, digest, fresh_output,
                      resolve_fit, run_r, write_json_atomic)

SOURCE_REVISION = 'b7b7001e56cea8e3931b09c917ea0e29e6cef3c6'
REFERENCE = f'https://github.com/AlexDiana/occJSDM/blob/{SOURCE_REVISION}/vignettes/occJSDM.Rmd'
quickstart_mcp = FastMCP(name='quickstart', strict_input_validation=True)


class StrictOptions(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)


class ListParams(StrictOptions):
    n_factors: Annotated[int, Field(ge=0, description='Number of native latent factors')]


class MCMCParams(StrictOptions):
    nchain: Annotated[int, Field(ge=2, description='Number of chains, at least two for diagnostics')]
    nburn: Annotated[int, Field(ge=0, description='Discarded burn-in iterations per chain')]
    niter: Annotated[int, Field(ge=1, description='Retained draws per chain; native very-short-fit errors are preserved')]
    nthin: Annotated[int, Field(ge=1, description='Iterations between retained draws')]


class FitOptions(StrictOptions):
    occCovariates: Annotated[list[str], Field(description='Selected occupancy columns in info, or an explicit empty list')]
    collCovariates: Annotated[list[str], Field(description='Selected collection columns in info, or an explicit empty list')]
    threshold: Annotated[int, Field(ge=1, description='Integer read threshold for replicated observations')]
    listParams: Annotated[ListParams, Field(description='Native latent-factor settings')]
    MCMCparams: Annotated[MCMCParams, Field(description='Explicit native chain settings')]


def response(message: str, native: dict, out: Path) -> dict:
    payload = {**native, 'schema_version':1, 'message':message, 'reference':REFERENCE,
               'artifacts':[artifact('log', out/'stdout.log','Captured native R console output'),
                            artifact('log', out/'stderr.log','Captured native R startup diagnostics')]}
    return payload


def covariates(occ: list[str], coll: list[str]) -> None:
    for values in (occ,coll):
        if any(not x.strip() for x in values) or len(values) != len(set(values)):
            raise ValueError('Covariate names must be nonempty and unique within each selection')


@quickstart_mcp.tool()
async def validate_data(
    data_path: Annotated[str, 'Project-local RDS containing native list(info, OTU, traits); no automatic repairs'],
    occCovariates: Annotated[list[str], 'Occupancy covariate column names in info, or an explicit empty list'],
    collCovariates: Annotated[list[str], 'Collection covariate column names in info, or an explicit empty list'],
) -> dict:
    """Check native input structure and eligibility for a non-spatial occJSDM fit.
    Input is a project-local RDS and covariate names; output is an eligibility report without observation values.
    """
    covariates(occCovariates,collCovariates)
    path = confined_path(data_path)
    out = fresh_output('validate_data')
    native = await run_r('validate_data',dict(data_path=str(path),occCovariates=occCovariates,collCovariates=collCovariates),out)
    return response('Data are eligible for the supported non-spatial workflow.' if native['eligible'] else 'Resolve the reported issues before fitting.',native,out)


@quickstart_mcp.tool()
async def fit_model(
    data_path: Annotated[str, 'Project-local RDS containing native occJSDM-ready data'],
    options: Annotated[FitOptions, 'Explicit covariates, threshold, latent factors and MCMC settings; unknown and spatial options rejected'],
    seed: Annotated[int, Field(ge=0, le=2147483647, strict=True, description='Seed passed to native R set.seed before fitting')],
) -> dict:
    """Fit a supported non-spatial occJSDM model with the native R sampler and default priors.
    Input is native RDS data, explicit options and seed; output is a completed fit identifier and hashed local artifacts.
    """
    covariates(options.occCovariates, options.collCovariates)
    path = confined_path(data_path)
    out = fresh_output('fit_model')
    # The hash belongs to the exact immutable bytes passed to R, even if the original changes during the fit.
    snapshot = out/'input.rds'
    shutil.copyfile(path,snapshot)
    input_hash = digest(snapshot)
    requested = options.model_dump()
    native = await run_r('fit_model',dict(data_path=str(snapshot),options=requested,seed=seed),out)
    pending = out/'fit.rds.pending'
    if not pending.is_file() or not native.get('readback_valid'):
        raise RuntimeError('Native fit did not pass RDS read-back; no fit identifier was published')
    fit_path = out/'fit.rds'
    os.replace(pending,fit_path)
    fit_id = 'fit_'+uuid.uuid4().hex
    effective = {**requested,'listParams':{'n_factors':native['effective_n_factors']},'spatCovariates':[],'summarisedLatentPresences':True,'listPriors':{}}
    manifest = dict(schema_version=1,status='complete',fit_id=fit_id,fit_path=str(fit_path),fit_sha256=digest(fit_path),
                    input_path=str(path),input_sha256=input_hash,seed=seed,model=native['model'],
                    effective_settings=effective,requested_settings=requested,source_revision=SOURCE_REVISION,runtime=native['runtime'],
                    python_version=sys.version,fastmcp_version=version('fastmcp'),pydantic_version=version('pydantic'),thread_settings=native['runtime']['thread_env'],warnings=native.get('warnings',[]),qualifications=native.get('qualifications',[]))
    manifest_path = out/'manifest.json'
    write_json_atomic(manifest_path,manifest)
    registry = confined_path('artifacts/registry',exists=False)
    registry.mkdir(parents=True,exist_ok=True)
    write_json_atomic(registry/(fit_id+'.json'),dict(manifest_path=str(manifest_path),manifest_sha256=digest(manifest_path)))
    result=response('Native non-spatial fit saved; inspect convergence before interpretation.',native,out)
    result.update(fit_id=fit_id,settings=effective,input_sha256=input_hash)
    result['artifacts'].extend([artifact('fit',fit_path,'Native occJSDM fit, reload with readRDS'),artifact('manifest',manifest_path,'Fit provenance and effective settings')])
    return result


@quickstart_mcp.tool()
async def diagnostics(
    fit_id: Annotated[str, 'Server-owned identifier returned by fit_model'],
    max_rows: Annotated[int, Field(ge=1,le=20,strict=True,description='Maximum diagnostic preview rows, capped at 20')] = 20,
) -> dict:
    """Inspect native convergence diagnostics for selected occupancy and detection coefficients.
    Input is a completed fit identifier; output is bounded screening results and the full native diagnostic CSV.
    """
    path, manifest = resolve_fit(fit_id)
    out = fresh_output('diagnostics')
    native = await run_r('diagnostics',dict(fit_path=str(path),max_rows=max_rows),out)
    result = response('Screening covers selected coefficient blocks; unavailable diagnostics are inconclusive.',native,out)
    result.update(fit_id=fit_id,model=manifest['model'],settings=manifest['effective_settings'],fit_qualifications=manifest.get('qualifications',[]))
    result['artifacts'].append(artifact('diagnostics',out/'diagnostics.csv','Complete native convergence diagnostic table'))
    return result


@quickstart_mcp.tool()
async def summarise_fit(
    fit_id: Annotated[str, 'Server-owned identifier returned by fit_model'],
    max_rows: Annotated[int, Field(ge=1,le=20,strict=True,description='Maximum baseline occupancy preview rows, capped at 20')] = 20,
) -> dict:
    """Summarise native baseline occupancy draws with convergence qualifications.
    Input is a completed fit identifier; output is bounded per-species baseline summaries and full local CSV tables.
    """
    path, manifest = resolve_fit(fit_id)
    out = fresh_output('summarise_fit')
    native = await run_r('summarise_fit',dict(fit_path=str(path),max_rows=max_rows),out)
    result = response('Baseline occupancy summarizes intercept-based draws; inspect diagnostic qualifications.',native,out)
    result.update(fit_id=fit_id,model=manifest['model'],settings=manifest['effective_settings'],fit_qualifications=manifest.get('qualifications',[]))
    result['artifacts'].extend([artifact('summary',out/'summary.csv','Full per-species baseline occupancy summary'),artifact('diagnostics',out/'diagnostics.csv','Full native diagnostic qualifications')])
    return result
