#!/usr/bin/env python3
"""Curate an explicitly identified Cloud project copy; never alter the source project.

The default action stages and reports. --apply-copy 13076906 publishes the staged
base in that copy and archives its prior contents in a private /tmp directory.
The numeric confirmation is operator confirmation, not automatic Cloud identity
discovery. Run only in the new project at the URL supplied by Doug.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
import uuid

PINS = {
    'source_sha256':'eb9f9004aec5c2b24e404aa7afd795792d9344083a1a441d1262b2ada57bed05',
    'source_revision':'b7b7001e56cea8e3931b09c917ea0e29e6cef3c6',
    'fit_id':'fit_fa9d078129f9450cabc029f0af695958',
    'fit_sha256':'446556b0ab908541034af3fb5369864e8164531086e25d3256d19f0078fbd5f8',
    'input_sha256':'86a183af12fa6c8e7f2d87abd01decc30f6280bb8d4002c8ccc9dee4fef7b584',
    'R':'R version 4.6.1 (2026-06-24)',
    'settings':dict(occCovariates=['X_psi.EnvCov.1','X_psi.EnvCov.2'],
                    collCovariates=['X_theta'], threshold=1, listParams={'n_factors':2},
                    MCMCparams={'nchain':2,'nburn':20,'niter':20,'nthin':1},
                    spatCovariates=[], summarisedLatentPresences=True, listPriors={}),
}


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def inside(root, value):
    path = Path(value).resolve(strict=True)
    if not path.is_relative_to(root.resolve()):
        raise ValueError('Reference artifact escapes the copied project')
    return path


def reference(root, pins):
    index_path = root / 'artifacts/registry' / (pins['fit_id'] + '.json')
    index = json.loads(index_path.read_text())
    manifest_path = inside(root, index['manifest_path'])
    if sha(manifest_path) != index['manifest_sha256']:
        raise ValueError('Reference manifest integrity check failed')
    m = json.loads(manifest_path.read_text())
    if (m['status'] != 'complete' or m['fit_id'] != pins['fit_id']
            or m['source_revision'] != pins['source_revision']):
        raise ValueError('Unexpected reference identity or source revision')
    fit = inside(root, m['fit_path'])
    data = inside(root, fit.parent / 'input.rds')
    if sha(fit) != pins['fit_sha256'] or m['fit_sha256'] != pins['fit_sha256']:
        raise ValueError('Reference fit integrity check failed')
    if sha(data) != pins['input_sha256'] or m['input_sha256'] != pins['input_sha256']:
        raise ValueError('Reference input integrity check failed')
    if 'settings' in pins and m['effective_settings'] != pins['settings']:
        raise ValueError('Reference settings differ from the frozen trial')
    if 'R' in pins and m['runtime']['R'] != pins['R']:
        raise ValueError('Reference R version differs from the frozen trial')
    return m, (index_path, manifest_path, fit, data)


def source_allowed(rel):
    p = Path(rel)
    if rel in {'DESCRIPTION','NAMESPACE','LICENSE','README.md','CITATION.cff',
               'data/sampledata.rda','data/sampleresults.rda',
               'vignettes/occJSDM.Rmd','vignettes/lesson-links.R','vignettes/teaching.css'}:
        return True
    if len(p.parts) == 2:
        return ((p.parts[0] == 'R' and p.suffix == '.R')
                or (p.parts[0] == 'man' and p.suffix == '.Rd')
                or (p.parts[0] == 'src' and (p.suffix in {'.cpp','.h'}
                                            or p.name in {'Makevars','Makevars.win'})))
    return rel.startswith('vignettes/occJSDM_files/figure-gfm/') and p.suffix == '.png'


def copy_file(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def stage_project(root, stage, pins, route):
    if route not in {'base','direct-R'}:
        raise ValueError('Unknown comparison route')
    if stage.exists():
        raise FileExistsError('Staging directory must be new')
    root = root.resolve()
    stage = stage.parent.resolve() / stage.name
    m, artifacts = reference(root, pins)
    runtime = root / 'mcp/paper2agent/cloud-r-runtime'
    archives = runtime / 'renv/cellar/occJSDM'
    archive = next((p for p in sorted(archives.glob('*.tar.gz'))
                    if sha(p) == pins['source_sha256']), None)
    if archive is None:
        raise ValueError('Pinned source archive integrity check failed')
    # Complete the checks before creating any staging output.
    for rel in ('activate.R','renv/activate.R','renv.lock',
                'bootstrap-library','library','sandbox','state'):
        if not (runtime / rel).exists():
            raise FileNotFoundError('Missing isolated runtime component: ' + rel)
    settings_path = root / '.posit/assistant/settings.json'
    config = json.loads(settings_path.read_text())
    if config.get('mcpServers',{}).get('occJSDM',{}).get('environment',{}).get('P2A_RSCRIPT'):
        raise ValueError('Remove the temporary/custom R launcher before comparison setup')
    stage.mkdir()
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            parts = Path(member.name).parts
            if not parts or parts[0] != 'occJSDM' or '..' in parts:
                continue
            rel = '/'.join(parts[1:])
            if not source_allowed(rel):
                continue
            if not member.isfile():
                raise ValueError('Curated source member is not a regular file: ' + rel)
            destination = stage / rel
            destination.parent.mkdir(parents=True, exist_ok=True)
            with tar.extractfile(member) as stream:
                destination.write_bytes(stream.read())
    target_runtime = stage / 'mcp/paper2agent/cloud-r-runtime'
    for rel in ('activate.R','renv/activate.R','renv/settings.json','renv.lock','source-provenance.json',
                'bootstrap-library','library','sandbox','state'):
        source = runtime / rel
        destination = target_runtime / rel
        if source.is_dir():
            shutil.copytree(source, destination, symlinks=True)
        elif source.is_file():
            copy_file(source, destination)
    for name in ('assistant-instructions.md','glossary.md'):
        copy_file(root / 'mcp/paper2agent/class' / name, stage / 'teaching' / name)
    if route == 'base':
        for name in ('src','occJSDM-env'):
            shutil.copytree(root / 'mcp/paper2agent' / name,
                            stage / 'mcp/paper2agent' / name, symlinks=True,
                            ignore=shutil.ignore_patterns('__pycache__','.pytest_cache'))
    else:
        config.get('mcpServers',{}).pop('occJSDM', None)
    settings_target = stage / '.posit/assistant/settings.json'
    settings_target.parent.mkdir(parents=True, mode=0o700)
    settings_target.write_text(json.dumps(config,indent=2)+'\n')
    settings_target.chmod(0o600)
    for path in artifacts:
        copy_file(path, stage / path.relative_to(root))
    copy_file(artifacts[2], stage / 'trial/fit-short.rds')
    copy_file(artifacts[3], stage / 'trial/data.rds')
    # Only the synthetic teaching inputs get agent-readable permissions.
    for path in (stage / 'trial').iterdir():
        path.chmod(0o644)
    (stage / 'occJSDM.Rproj').write_text('Version: 1.0\n\nRestoreWorkspace: No\n'
                                      'SaveWorkspace: No\nAlwaysSaveHistory: No\n')
    (stage / '.Rprofile').write_text('Sys.setenv(P2A_R_PROJECT="/cloud/project/mcp/paper2agent/cloud-r-runtime")\n'
        'source(file.path(Sys.getenv("P2A_R_PROJECT"), "activate.R"))\n'
        'stopifnot(R.version.string == "R version 4.6.1 (2026-06-24)",\n'
        '  startsWith(find.package("occJSDM"), Sys.getenv("P2A_R_PROJECT")))\n')
    return m


VARIANTS_R = r'''
args <- commandArgs(trailingOnly=TRUE)
Sys.setenv(P2A_R_PROJECT=args[[1]])
source(file.path(args[[1]], "activate.R"))
stopifnot(R.version.string == args[[4]])
library(occJSDM)
revision <- as.character(read.dcf(file.path(find.package("occJSDM"),"DESCRIPTION"),fields="RemoteSha")[,"RemoteSha"])
stopifnot(identical(unique(revision[!is.na(revision)]), "b7b7001e56cea8e3931b09c917ea0e29e6cef3c6"))
trial <- file.path(args[[2]], "trial")
grader <- args[[3]]
d <- readRDS(file.path(trial, "data.rds"))
stopifnot(nrow(d$info)==160L, nrow(d$OTU)==160L, ncol(d$OTU)==4L)
malformed <- d
malformed$OTU <- malformed$OTU[-nrow(malformed$OTU), , drop=FALSE]
saveRDS(malformed, file.path(trial, "malformed.rds"))
no_traits <- d
no_traits$traits <- NULL
saveRDS(no_traits, file.path(trial, "no-traits.rds"))
counts <- d
take <- !duplicated(d$info$Site)
counts$info <- d$info[take, , drop=FALSE]
counts$OTU <- d$OTU[take, , drop=FALSE] + 5L
stopifnot(!anyDuplicated(counts$info$Site), !anyDuplicated(counts$info$Sample), max(counts$OTU)>1)
msg <- tryCatch(occJSDM:::inferDataModel(counts),error=function(e)conditionMessage(e))
stopifnot(grepl("Counts model not supported",msg,fixed=TRUE))
saveRDS(counts, file.path(trial, "unreplicated-counts.rds"))
set.seed(1704)
one <- runOccJSDM(d, listParams=list(n_factors=2L), threshold=1L,
  occCovariates=c("X_psi.EnvCov.1","X_psi.EnvCov.2"), collCovariates="X_theta",
  spatCovariates=NULL, MCMCparams=list(nchain=1L,nburn=20L,niter=20L,nthin=1L),
  summarisedLatentPresences=TRUE)
unavailable <- returnConvergenceDiagnostics(one)
stopifnot(all(is.na(unavailable$rhat)))
write.csv(unavailable,file.path(trial,"diagnostics-unavailable.csv"),row.names=FALSE,na="")
fit <- readRDS(file.path(trial,"fit-short.rds"))
diagnostics <- returnConvergenceDiagnostics(fit)
write.csv(diagnostics,file.path(grader,"diagnostics.csv"),row.names=FALSE,na="")
draws <- returnOccupancyRates(fit)
summary <- data.frame(species=colnames(draws),mean=colMeans(draws),
  q2.5=apply(draws,2,quantile,0.025),q97.5=apply(draws,2,quantile,0.975),row.names=NULL)
write.csv(summary,file.path(grader,"baseline-occupancy.csv"),row.names=FALSE)
jsonlite::write_json(list(R=R.version.string,package_path=find.package("occJSDM"),
  diagnostics=list(rows=nrow(diagnostics),rhat_above_1_01=sum(diagnostics$rhat>1.01,na.rm=TRUE),
    ess_below_400=sum(diagnostics$ess<400,na.rm=TRUE),rhat_unavailable=sum(is.na(diagnostics$rhat)),
    ess_unavailable=sum(is.na(diagnostics$ess)))),file.path(grader,"native-reference.json"),
  pretty=TRUE,auto_unbox=TRUE)
'''


def generate_variants(root, stage, backup, manifest):
    script = backup / 'generate-variants.R'
    script.write_text(VARIANTS_R)
    grader = backup / 'grader'
    grader.mkdir()
    with (backup / 'native-setup.log').open('w') as log:
        result = subprocess.run([manifest['runtime']['Rscript'],'--vanilla',str(script),
            str(root/'mcp/paper2agent/cloud-r-runtime'),str(stage),str(grader),
            manifest['runtime']['R']], stdout=log,stderr=subprocess.STDOUT,timeout=600,cwd=root)
    if result.returncode:
        raise RuntimeError('Native fixture preparation failed; original project unchanged. Log: '
                           + str(backup/'native-setup.log'))
    for path in (stage / 'trial').iterdir():
        path.chmod(0o644)


def check_retained_links(root, stage):
    """Reject links into excluded project content before publication."""
    root, stage = root.resolve(), stage.resolve()
    count = 0
    for path in stage.rglob('*'):
        if not path.is_symlink():
            continue
        count += 1
        target = Path(os.readlink(path))
        if not target.is_absolute():
            target = path.parent / target
        target = Path(os.path.abspath(target))
        # /tmp and /var aliases need canonical parents without dereferencing the
        # final link into content that is about to be excluded.
        if target.is_relative_to(stage):
            retained = target
        elif target.is_relative_to(root):
            retained = stage / target.relative_to(root)
        else:
            retained = target
        if not retained.exists():
            raise ValueError('Runtime link is broken or points to excluded content: '
                             + str(path.relative_to(stage)))
    return count


def tree_record(path):
    """Private backup verification; never print credential-bearing file hashes."""
    if path.is_symlink():
        return ('symlink', os.readlink(path))
    if path.is_file():
        return ('file', sha(path), path.stat().st_mode & 0o777)
    if path.is_dir():
        return ('directory', {p.name:tree_record(p) for p in sorted(path.iterdir())})
    raise ValueError('Unsupported filesystem entry: ' + str(path))


def archive_copy(source, destination):
    if source.is_symlink():
        destination.symlink_to(os.readlink(source))
    elif source.is_dir():
        shutil.copytree(source, destination, symlinks=True)
    else:
        shutil.copy2(source, destination)
    if tree_record(source) != tree_record(destination):
        raise RuntimeError('Backup verification failed; copied project unchanged')


def publish_stage(root, stage, backup):
    root = root.resolve()
    stage = stage.resolve()
    if backup.resolve().is_relative_to(root):
        raise ValueError('The backup must be outside the model project')
    if stage.parent != root or stage.stat().st_dev != root.stat().st_dev:
        raise ValueError('Staging must be directly inside the copied project')
    original = backup / 'original'
    original.mkdir(mode=0o700)
    old_paths = [p for p in sorted(root.iterdir()) if p != stage]
    # Cross-filesystem archival is copy-only. Verify every complete copy before
    # changing any source entry; never use copy-and-delete shutil.move here.
    for path in old_paths:
        archive_copy(path, original/path.name)
    holding = root / ('.comparison-original-' + uuid.uuid4().hex)
    holding.mkdir(mode=0o700)
    moved_old, moved_new = [], []
    try:
        for path in old_paths:
            os.rename(path, holding/path.name)
            moved_old.append(path.name)
        for path in sorted(stage.iterdir()):
            os.rename(path, root/path.name)
            moved_new.append(path.name)
    except BaseException:
        # Same-filesystem renames are atomic. The complete private archive also
        # remains available independently of rollback.
        for name in reversed(moved_new):
            os.rename(root/name, stage/name)
        for name in reversed(moved_old):
            if os.path.lexists(root/name):
                raise RuntimeError('Rollback conflict; original retained at ' + str(holding/name))
            os.rename(holding/name, root/name)
        holding.rmdir()
        raise
    stage.rmdir()
    try:
        # renv sandbox directories can be read-only. Change only the discarded
        # duplicate after its complete archive has been verified; never follow
        # links or change permissions in the retained runtime or private backup.
        for directory, _, _ in os.walk(holding, followlinks=False):
            path = Path(directory)
            if not path.is_symlink():
                path.chmod(path.stat().st_mode | 0o700)
        shutil.rmtree(holding)
    except OSError as error:
        raise RuntimeError('Publication finished but excluded-material cleanup failed. '
            'Do not start trials. Complete private backup: ' + str(original)
            + '; remaining excluded material: ' + str(holding)) from error


def verify_deployment(root, backup, manifest):
    """Read-only native activation and MCP metadata check; no model/provider call."""
    script = backup / 'verify-published.R'
    script.write_text('''args <- commandArgs(TRUE)
Sys.setenv(P2A_R_PROJECT=args[1]); source(file.path(args[1],"activate.R"));library(occJSDM)
revision <- as.character(read.dcf(file.path(find.package("occJSDM"),"DESCRIPTION"),fields="RemoteSha")[,"RemoteSha"])
stopifnot(R.version.string==args[2], startsWith(find.package("occJSDM"),args[1]),
 identical(unique(revision[!is.na(revision)]),args[3]))
fit <- readRDS(args[4]); stopifnot(is.list(fit),"results_output" %in% names(fit))
cat("Published native runtime and saved fit loaded successfully\\n")
''')
    with (backup/'published-runtime.log').open('w') as log:
        result = subprocess.run([manifest['runtime']['Rscript'],'--vanilla',str(script),
            str(root/'mcp/paper2agent/cloud-r-runtime'),manifest['runtime']['R'],
            PINS['source_revision'],str(root/'trial/fit-short.rds')],
            stdout=log,stderr=subprocess.STDOUT,cwd=root,timeout=60)
    if result.returncode:
        raise RuntimeError('Published runtime check failed; do not start trials. Private log: '
                           + str(backup/'published-runtime.log'))
    script = backup/'verify-tools.py'
    script.write_text('''import asyncio,json,sys
from pathlib import Path
from fastmcp import Client
from fastmcp.client.transports import StdioTransport
root=Path(sys.argv[1])
server=json.loads((root/".posit/assistant/settings.json").read_text())["mcpServers"]["occJSDM"]
command=server["command"]
async def check():
    transport=StdioTransport(command=command[0],args=command[1:],
        env=server.get("environment",{}),cwd=str(root),keep_alive=False,
        log_file=Path(sys.argv[2]))
    async with Client(transport) as client:
        names=sorted(t.name for t in await client.list_tools())
        assert names==["diagnostics","fit_model","summarise_fit","validate_data"],names
        print(json.dumps(names))
asyncio.run(check())
''')
    with (backup/'published-tool-inventory.json').open('w') as log:
        result = subprocess.run([str(root/'mcp/paper2agent/occJSDM-env/bin/python'),
            str(script),str(root),str(backup/'published-mcp.log')],stdout=log,
            stderr=subprocess.STDOUT,cwd=root,timeout=60)
    if result.returncode:
        raise RuntimeError('Published tool inventory failed; do not start trials. Private log: '
                           + str(backup/'published-tool-inventory.json'))
    return json.loads((backup/'published-tool-inventory.json').read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply-copy', choices=['13076906'])
    args = parser.parse_args()
    root = Path('/cloud/project')
    stage = root / ('.comparison-stage-' + uuid.uuid4().hex)
    backup = Path(tempfile.mkdtemp(prefix='occjsdm-comparison-',dir='/tmp'))
    manifest = stage_project(root,stage,PINS,'base')
    symlinks = check_retained_links(root,stage)
    generate_variants(root,stage,backup,manifest)
    if sha(stage/'trial/fit-short.rds') != PINS['fit_sha256']:
        raise ValueError('Saved fit changed during preparation')
    # The grader inventory excludes credential-bearing settings and dependency caches.
    inventory = {str(p.relative_to(stage)):sha(p) for p in sorted(stage.rglob('*'))
        if p.is_file() and p.relative_to(stage).parts[0] in
        {'R','src','man','vignettes','data','teaching','trial','artifacts'}}
    report = dict(status='staged_pending_publication',project_url=
        'https://posit.cloud/spaces/828411/content/13076906',source_revision=PINS['source_revision'],
        fit_sha256=PINS['fit_sha256'],input_sha256=PINS['input_sha256'],
        student_files_sha256=inventory,retained_symlinks_checked=symlinks,teaching_review='draft_pending_review',
        grader_location=str(backup/'grader'),backup_location=str(backup),
        original_source_project_untouched=True,condition='base_not_a_scored_trial')
    (backup/'grader/student-project-manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    if args.apply_copy:
        publish_stage(root,stage,backup)
        report['tool_inventory']=verify_deployment(root,backup,manifest)
        if sha(root/'trial/fit-short.rds') != PINS['fit_sha256']:
            raise ValueError('Published saved fit integrity failed; do not start trials')
        report['status']='curated_base_published_restart_R_and_Assistant_required'
    (backup/'grader/student-project-manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:value for key,value in report.items() if key!='student_files_sha256'},indent=2))
    print('Student files:',len(inventory))
    print('Native reference:',(backup/'grader/native-reference.json').read_text())
    print('Private backup is temporary; the original Cloud project remains the durable backup.')


if __name__ == '__main__':
    main()
