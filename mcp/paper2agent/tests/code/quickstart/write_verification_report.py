"""Seal the fresh verifier report only after all independent checks pass."""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'reports/verification/quickstart'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
production=read(OUT/'producer_hashes.json')
assert all(sha(ROOT/p)==h for p,h in production.items())
source=read(ROOT/'notebooks/quickstart/source_hashes_before.json')
assert all(sha(ROOT/'repo/occJSDM'/p)==h for p,h in source.items())
immutable=read(OUT/'immutable_before.json')
assert all(sha(p)==h for p,h in immutable.items())
assert read(OUT/'numerical_audit.json')['all_passed']
assert read(OUT/'acceptance_module.json')['success']
assert read(OUT/'acceptance_artifact_audit_guarded.json')['success']
assert '37 passed' in (OUT/'pytest-attempt2.log').read_text()
commands=[
 ('Rscript --vanilla tests/code/quickstart/independent_prepare.R PROJECT_ROOT',0,'prepare.log'),
 ('Rscript --vanilla tests/code/quickstart/independent_prepare.R PROJECT_ROOT',0,'prepare-attempt2.log'),
 ('occJSDM-env/bin/python -m pytest -q tests/code/quickstart',1,'pytest-attempt1.log'),
 ('Independent numerical_audit fixture invoked against retained MCP exercise outputs',0,'audit-run-attempt1.log'),
 ('occJSDM-env/bin/python -m pytest -q tests/code/quickstart',0,'pytest-attempt2.log'),
 ('Independent numerical_audit fixture invoked after adding empty nonspatial Xs comparison',0,'audit-final-command.log'),
 ('occJSDM-env/bin/python tests/code/quickstart/exercise_acceptance_module.py',0,'acceptance-module.log'),
 ('Rscript --vanilla tests/code/quickstart/audit_acceptance_artifacts.R PROJECT_ROOT CUTOFF R_PROJECT OUTPUT_JSON',0,'acceptance-artifact-audit.log'),
 ('occJSDM-env/bin/python tests/code/quickstart/audit_acceptance_artifacts.py PROJECT_ROOT PRIOR_DIRS_JSON R_PROJECT OUTPUT_JSON',0,'acceptance-artifact-final-command.log'),
 ('Same comparator with exact options/input/publication/producer hashes and isolated-source guards',0,'acceptance-artifact-guarded-command.log')]
report={
 'schema_version':1,'module':'quickstart','status':'verified','implementation_run_id':'extract-quickstart-attempt-1',
 'verified_at':datetime.now(timezone.utc).isoformat(),'project_root':str(ROOT),'source_revision':'b7b7001e56cea8e3931b09c917ea0e29e6cef3c6',
 'verifier_identity':'Fresh independent /root/mcp_verifier, distinct from extraction implementer',
 'exposed_tools':['validate_data','fit_model','diagnostics','summarise_fit'],'exclusions':[],
 'tested_files':production,'tested_production_hashes':production,
 'source_call_review':{
  'validate_data':'Fixed R dispatcher calls native inferDataModel, process_covariates and create_covariates_matrix; cheap structural/alignment/invariance guards prevent native errors or misleading first-row preparation. No sampler or scientific formula is copied.',
  'fit_model':'Native runOccJSDM called once per invocation, set.seed applied first, no spatial covariates, summarisedLatentPresences TRUE and unchanged default priors. Input snapshot, atomic publication, RDS read-back and hashed registry are transport adaptations.',
  'diagnostics':'Native returnConvergenceDiagnostics table is written in full; selected-block screen only counts Rhat>1.01, ESS<400 and unavailable entries. JSON NA becomes null; no convergence statistic is reimplemented.',
  'summarise_fit':'Native returnOccupancyRates returns pooled draws by species; approved R mean and quantile output adaptation is validated on all draws, full summaries and labels. Estimand is inverse-logit species intercept, not a site prediction or observed detection count.',
  'runtime.py':'Explicit isolated activation under --vanilla, pinned source guard, fixed argv/data JSON, project-confined resolved paths, captured diagnostics, timeout and cancellation reaping, no untracked user-library fallback.',
  '__init__.py':'Documentation-only package marker inspected; imported during every production call.',
  'installed_package':'All six relevant installed function bodies and formals exactly match the pinned R source. Existing recorded source installation, origin and compile evidence were reviewed; every pinned source file remains unchanged.'},
 'scientific_comparisons':read(OUT/'numerical_audit.json'),
 'changed_input_oracles':{
  'binary_single':'40 observations, one renamed species, no traits, changed numeric covariate, seed2401, zero factors, threshold1,20burn/20retained/2chains/thin2. Full native fit and 40x1 baseline draw matrix match direct R.',
  'occupancy_changed':'Three renamed species, numeric and categorical occupancy covariates with explicit factor levels, changed collection covariate, seed2402, one factor, threshold2 and thin2; full native fit/classes/labels/diagnostics/baseline match direct R.',
  'two_stage_changed':'Three renamed species, changed numeric/categorical occupancy and collection covariates, seed2404, one factor, threshold3 and thin3; full native fit/classes/labels/diagnostics/baseline match direct R.',
  'two_stage_single':'Direct native one-species two_stage/zero-factor call fails with incorrect number of dimensions. MCP preserves that error and publishes no fit. The failure is retained as a documented native boundary, not silently skipped.'},
 'independent_tests':{'test_files':[str(p.relative_to(ROOT)) for p in sorted((ROOT/'tests/code/quickstart').glob('test_*.py'))],
  'count':37,'passed':37,'elapsed_seconds':77.31,'framework':'project Python3.14.6/FastMCP4.0.3 Client(module)',
  'checks':['full native arrays and all six top-level fit components','all18savednativecalls independently rerun in one activation','immutable source and reference hashes','source functions/classes/dimensions/design/species labels','strict schemas and native error transport','closed unsupported options and no spatial/custom priors','changed settings and seeds','single-species dimensional simplification','all diagnostic and baseline CSV values','constant chains Rhat null/ESS0/unavailable counts','preview1/2/20 and errors0/21/noninteger','exact row identity and trait-alias semantics','composed preflight failure/no publication','unchanged original inputs','repeated fit/diagnostic/summary artifact separation','publication and integrity hashes','unavailable configured Rscript','persistent20/20 scientific-unsuitability qualifications']},
 'reused_inspected_evidence':{'report':'reports/implementation-quickstart.json','frozen_suite':'reports/implementation-green-suite.json','count':66,
  'qualification':'These are inspected earlier broad tests, not additional fresh-verifier tests. Matching frozen runtime.py/R hashes and only a Path import repair permit reuse for unchanged logic.',
  'coverage':['17 native validation/error or misleading-success guards','native two_stage NA support','traits reordering/categorical encoding','local/global/character and unbalanced alignment','timeout/cancellation process-group reaping','symlink confinement','effective factor cap','two fresh changed-settings native comparisons']},
 'repairs':[{'file':'src/tools/quickstart.py','attempt':1,'change':'Import pathlib.Path for the private response annotation, which Python3.14 lazy annotations masked during calls. get_type_hints and FastMCP schemas now pass. No scientific behavior changed.'}],
 'attempts':{'pytest_runs':2,'production_repairs_per_tool':{'validate_data':1,'fit_model':1,'diagnostics':1,'summarise_fit':1},'maximum_per_tool':6},
 'retained_failures':[{'log':'reports/verification/quickstart/pytest-attempt1.log','outcome':'30passed,7setup errors: verifier selected historical manifests lacking requested_settings. Corrected to require complete metadata and exact input/options/seed/hash match for every case. No assertion relaxed or case excluded.'},
  {'log':'reports/verification/quickstart/prepare.log','outcome':'Fresh direct R identified native single-species two_stage dimension error; preserved in independent_inputs.json and actual MCP acceptance.'}],
 'runtime_acceptance':{'cases':'reports/mcp-acceptance-quickstart.json','count':19,'positive':12,'negative':7,
  'module_execution':'reports/verification/quickstart/acceptance_module.json','module_success':True,
  'scientific_artifacts':'reports/verification/quickstart/acceptance_artifact_audit_guarded.json','full_native_fits_compared':4,'full_native_tables_compared':5,'publication_hash_checks':4,
  'dependency':'Downstream cases intentionally use explicit independently tested server-owned registry IDs in the same project root. They do not rely on pytest or case ordering; this is scoped runtime acceptance, not a portable fixture import API.',
  'oracle_mapping':'reports/verification/quickstart/acceptance_oracle_mapping.json',
  'comparator_command':'PROJECT_PYTHON tests/code/quickstart/audit_acceptance_artifacts.py PROJECT_ROOT COORDINATOR_PRIOR_DIR_NAMES_JSON SELECTED_R_PROJECT COORDINATOR_OUTPUT_JSON',
  'comparator_output_note':'Coordinator should save snapshots and outputs under reports/coordinator-stdio-*.json, outside the immutable verifier namespace.',
  'remaining_coordinator_work':['root stdio integration acceptance','fresh Python plus separately restored R runtime acceptance','scientific artifact comparisons after each acceptance run']},
 'commands':[{'command':c,'exit_code':e,'log':'reports/verification/quickstart/'+log,'observed_log_mtime_utc':datetime.fromtimestamp((OUT/log).stat().st_mtime,timezone.utc).isoformat()} for c,e,log in commands],
 'read_only_shared_dependencies':{str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'r-runtime/activate.R',ROOT/'r-runtime/renv/activate.R',ROOT/'r-runtime/renv.lock',ROOT/'reports/tool-bindings.json',ROOT/'reports/executed_notebook_quickstart.json',ROOT/'reports/implementation-quickstart.json']},
 'immutable_evidence':{'pinned_source_files_checked':len(source),'reference_artifacts_checked':len(immutable),'reference_hashes':'reports/verification/quickstart/immutable_before.json','all_unchanged':True},
 'limitations':['Nonspatial binary/occupancy/two_stage only; continuous/no-stage counts/spatial/prediction/conversion out of scope.','Native single-species two_stage fails for tested settings and is not suitable for a positive classroom example; binary single-species is verified.','Faithful beta execution does not prove unbiased inference, convergence or calibrated intervals.20/20 fits are explicitly unsuitable for scientific interpretation.','Python must be>=3.11 because hashlib.file_digest is used; only3.14.6 is tested here. FastMCP metadata alone does not establish Python3.10 wrapper support.','Mac arm64 R4.5.0/project renv only; Linux, Posit Cloud, paid calls and classroom trial untested.','Module execution is not stdio or clean-runtime validation; those remain coordinator responsibilities.']}
report['verification_evidence_hashes']={str(p.relative_to(ROOT)):sha(p) for p in OUT.glob('*') if p.is_file()}
report['test_code_hashes']={str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'tests/code/quickstart').glob('*') if p.is_file()}
with (ROOT/'reports/verification-quickstart.json').open('x') as stream:stream.write(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps({'status':'verified','report':str(ROOT/'reports/verification-quickstart.json'),'sha256':sha(ROOT/'reports/verification-quickstart.json'),'tested_files':production,'verified_at':report['verified_at']}))
