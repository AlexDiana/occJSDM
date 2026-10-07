"""Execute actual acceptance cases through installed Client(module) before handoff."""
import asyncio
import json
import os
from pathlib import Path
import sys
import time
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,'/Users/douglasyu/.codex/skills/paper2agent/paper2mcp/scripts')
from verify_mcp_server import check_subset,check_artifacts,validate_cases
from fastmcp import Client
from tools.quickstart import quickstart_mcp
os.environ.update(P2A_PROJECT_ROOT=str(ROOT),P2A_RSCRIPT='/usr/local/bin/Rscript',P2A_R_PROJECT=str(ROOT/'r-runtime'))
cases=json.loads((ROOT/'reports/mcp-acceptance-quickstart.json').read_text())
validate_cases(cases,['validate_data','fit_model','diagnostics','summarise_fit'],True)
report={'mode':'installed FastMCP Client(module)','cutoff':time.time(),'cases':[]}
async def main():
    seen=set()
    async with Client(quickstart_mcp) as client:
        for c in cases:
            row={'name':c['name'],'success':False};report['cases'].append(row)
            try:
                r=await client.call_tool(c['tool'],c['arguments'],raise_on_error=False)
                if 'error_contains' in c:
                    msg='\n'.join(getattr(x,'text','') for x in r.content)
                    assert r.is_error and c['error_contains'] in msg,msg
                else:
                    assert not r.is_error,r
                    check_subset(r.data,c['expected_subset'])
                    paths=check_artifacts(r.data,c,seen);seen.update(paths)
                    row['result']=r.data
                row['success']=True
            except Exception as e:row['error']=str(e)
asyncio.run(main())
report['success']=all(r['success'] for r in report['cases'])
(ROOT/'reports/verification/quickstart/acceptance_module.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps({'success':report['success'],'cases':len(cases),'cutoff':report['cutoff']}))
raise SystemExit(0 if report['success'] else 1)
