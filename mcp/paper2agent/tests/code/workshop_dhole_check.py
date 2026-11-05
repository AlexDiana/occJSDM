"""Run Session 2's dhole MCP steps through the real tools and time them.

Not collected by pytest (no test_ prefix). Usage, from mcp/paper2agent:
  ./occJSDM-env/bin/python tests/code/workshop_dhole_check.py PATH/TO/dhole_input.rds OUTPUT.json
The input is copied into tmp/ because the server only reads files inside its project root.
"""
import asyncio
import json
import os
from pathlib import Path
import shutil
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
os.environ.setdefault('P2A_PROJECT_ROOT', str(ROOT))
os.environ.setdefault('P2A_R_PROJECT', str(ROOT / 'r-runtime'))
os.environ.setdefault('P2A_RSCRIPT', '/usr/local/bin/Rscript')
from fastmcp import Client
from tools.quickstart import quickstart_mcp

# The covariates and schedule Session 2 asks students to give the assistant.
OCC = ['annual_precipitation', 'mean_temperature', 'Human_modification_1km', 'Percent_of_forest_cover_1km',
       'Patch_density_1km', 'Edge_density_1km', 'streamOrigin', 'Regional_Division']
COLL = ['Filter_vol', 'PrevDayRain']
MCMC = dict(nchain=2, nburn=5000, niter=5000, nthin=1)


async def main(source: Path, report: Path) -> None:
    dest = ROOT / 'tmp' / 'workshop' / 'dhole_input.rds'
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, dest)
    out = {'input': str(source), 'settings': MCMC, 'occCovariates': OCC, 'collCovariates': COLL}
    async with Client(quickstart_mcp) as client:
        async def call(tool, args):
            start = time.monotonic()
            result = (await client.call_tool(tool, args)).data
            return result, round(time.monotonic() - start, 1)
        valid, out['validate_seconds'] = await call('validate_data', dict(data_path=str(dest), occCovariates=OCC, collCovariates=COLL))
        out['validate'] = {k: valid.get(k) for k in ('eligible', 'model', 'issues', 'notes', 'dimensions')}
        fit, out['fit_seconds'] = await call('fit_model', dict(data_path=str(dest), seed=2026, options=dict(
            occCovariates=OCC, collCovariates=COLL, threshold=1, listParams={'n_factors': 0}, MCMCparams=MCMC)))
        out['fit'] = {k: fit.get(k) for k in ('fit_id', 'model', 'readback_valid', 'kept_draws_per_chain', 'warning_count', 'warnings')}
        diag, out['diagnostics_seconds'] = await call('diagnostics', dict(fit_id=fit['fit_id'], max_rows=20))
        out['screening'] = diag['screening']
        summ, out['summary_seconds'] = await call('summarise_fit', dict(fit_id=fit['fit_id'], max_rows=5))
        out['summary_rows'] = summ['rows']
        out['fit_path'] = next(a['path'] for a in fit['artifacts'] if a['kind'] == 'fit')
    report.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    asyncio.run(main(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()))
