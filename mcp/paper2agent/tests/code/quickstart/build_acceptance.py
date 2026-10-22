"""Build stable acceptance expectations from fresh direct-R oracle tables."""
import csv
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'reports/verification/quickstart'
specs=json.loads((OUT/'independent_inputs.json').read_text())
actual=json.loads((OUT/'mcp_exercised_results.json').read_text())
constant=json.loads((OUT/'constant_diagnostics_result.json').read_text())
required={'validate_data':['data_path','occCovariates','collCovariates'],'fit_model':['data_path','options','seed'],'diagnostics':['fit_id'],'summarise_fit':['fit_id']}
cases=[]
def positive(name,tool,args,expected,count):
    cases.append(dict(name=name,tool=tool,arguments=args,required_inputs=required[tool],expected_subset=expected,min_artifacts=count,artifact_root=str(ROOT/'artifacts/calls'),unique_artifacts=True))
def negative(name,tool,args,error):
    cases.append(dict(name=name,tool=tool,arguments=args,required_inputs=required[tool],error_contains=error))
def screening(path):
    rows=list(csv.DictReader(Path(path).open()))
    finite=lambda v:v not in ('NA','NaN','Inf','-Inf','')
    rh=[finite(r['rhat']) and float(r['rhat'])>1.01 for r in rows]
    es=[finite(r['ess']) and float(r['ess'])<400 for r in rows]
    return dict(total_rows=len(rows),rhat_flagged=sum(rh),ess_flagged=sum(es),flagged_rows=sum(a or b for a,b in zip(rh,es)),
                unavailable_rhat=sum(not finite(r['rhat']) for r in rows),unavailable_ess=sum(not finite(r['ess']) for r in rows),
                unavailable_rows=sum(not finite(r['rhat']) or not finite(r['ess']) for r in rows),thresholds={'rhat':1.01,'ess':400})
c=specs['binary_single'];o=c['options']
for suffix in ('first','repeat'):
    positive('validate-binary-single-'+suffix,'validate_data',dict(data_path=c['data_path'],occCovariates=o['occCovariates'],collCovariates=o['collCovariates']),
             dict(eligible=True,model='binary',dimensions=[40,1],traits_present=False,issues=[],design={'species_count':1,'n_sites':40,'missing_observations':0}),2)
for name in ('binary_single','occupancy_changed','two_stage_changed','binary_single'):
    c=specs[name]
    positive('fit-'+name+('-repeat' if name=='binary_single' and any(x['name']=='fit-binary_single' for x in cases) else ''),'fit_model',
             dict(data_path=c['data_path'],options=c['options'],seed=c['seed']),
             dict(model=c['model'],readback_valid=True,effective_n_factors=c['options']['listParams']['n_factors'],kept_draws_per_chain=20,
                  total_iterations_per_chain=20+20*c['options']['MCMCparams']['nthin'],dimensions=c['dimensions'],runtime={'source_revision':'b7b7001e56cea8e3931b09c917ea0e29e6cef3c6'}),4)
fid=actual['binary_single']['fit_model']['fit_id']
for suffix in ('first','repeat'):
    positive('diagnostics-single-'+suffix,'diagnostics',dict(fit_id=fid,max_rows=2),
             dict(model='binary',total_rows=2,screening=screening(OUT/'binary_single_oracle_diagnostics.csv')),3)
    positive('summary-single-'+suffix,'summarise_fit',dict(fit_id=fid,max_rows=20),
             dict(model='binary',total_rows=1,rows=[{'species':'verifier_species_1'}],screening=screening(OUT/'binary_single_oracle_diagnostics.csv')),4)
positive('diagnostics-constant-chain','diagnostics',dict(fit_id=constant['fit_id'],max_rows=20),
         dict(screening=screening(ROOT/'notebooks/quickstart/reference/constant_species_diagnostics.csv')),3)
negative('path-escape','validate_data',dict(data_path='/etc/passwd',occCovariates=[],collCovariates=[]),'Path escapes')
negative('unknown-fit','summarise_fit',dict(fit_id='../fit.rds'),'Unknown fit_id')
negative('diagnostics-max21','diagnostics',dict(fit_id=fid,max_rows=21),'validation')
negative('summary-max21','summarise_fit',dict(fit_id=fid,max_rows=21),'validation')
c=specs['binary_single'];o={**c['options'],'spatCovariates':['x','y']}
negative('unsupported-spatial','fit_model',dict(data_path=c['data_path'],options=o,seed=c['seed']),'validation')
c=specs['two_stage_single']
negative('native-single-two-stage-dimension-error','fit_model',dict(data_path=c['data_path'],options=c['options'],seed=c['seed']),c['error'])
c=specs['occupancy_changed']
negative('composed-row-identity-failure','fit_model',dict(data_path=str(OUT/'row_mismatch.rds'),options=c['options'],seed=c['seed']),'row identifiers')
positive('traits-alias-report','validate_data',dict(data_path=str(OUT/'traits_alias.rds'),occCovariates=['habitat_changed'],collCovariates=[]),dict(eligible=False,model=None),2)
(ROOT/'reports/mcp-acceptance-quickstart.json').write_text(json.dumps(cases,indent=2)+'\n')
mapping={c['name']:{'input_path':c['arguments']['data_path'],'oracle_fit':str(OUT/(c['name'].removeprefix('fit-').removesuffix('-repeat')+'_oracle_fit.rds'))} for c in cases if c['tool']=='fit_model' and 'error_contains' not in c}
(OUT/'acceptance_oracle_mapping.json').write_text(json.dumps({'fits':mapping,'known_fit_ids':{n:r['fit_model']['fit_id'] for n,r in actual.items() if 'fit_model' in r},'constant_fit_id':constant['fit_id'],'static_dependency':'Downstream cases explicitly consume independently tested registry IDs in this same project root; they do not depend on pytest or runtime case ordering.','scientific_comparator':'tests/code/quickstart/audit_acceptance_artifacts.R PROJECT_ROOT CUTOFF_UNIX_SECONDS R_PROJECT'},indent=2)+'\n')
print(json.dumps({'cases':len(cases),'positive':sum('error_contains' not in c for c in cases)}))
