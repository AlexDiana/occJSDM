#!/usr/bin/env python3
"""Resumable local study. Workers and outputs remain inside the chosen run root."""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'md5').hexdigest()


def complete(job, root):
    output = root / 'jobs' / job['job']
    if not (output / 'status.json').exists():
        return False
    status = json.loads((output / 'status.json').read_text())
    if status['input_md5'] != job['input_md5']:
        raise RuntimeError(f"Input identity changed for {job['job']}")
    if digest(output / 'result.rds') != status['result_md5']:
        raise RuntimeError(f"Incomplete or modified result for {job['job']}")
    return True  # Failed fits count as attempted, not as permission to silently redraw/retry.


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('root', type=Path)
    parser.add_argument('code', type=Path)
    parser.add_argument('--launcher', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--score', action='store_true')
    parser.add_argument('--jobs', help='Comma-separated manifest row numbers for implementation checks')
    args = parser.parse_args()
    root, code = args.root.resolve(), args.code.resolve()
    all_jobs = list(csv.DictReader((root / 'manifest.csv').open()))
    selected = set(args.jobs.split(',')) if args.jobs else None
    jobs = [j for j in all_jobs if selected is None or j['job_number'] in selected]
    def needs_work(job):
        if not complete(job, root): return True
        output = root / 'jobs' / job['job']
        status = json.loads((output / 'status.json').read_text())
        return args.score and status['ok'] and not (output / 'score.rds').exists()
    pending = [j for j in jobs if needs_work(j)]
    # An atomic lock prevents concurrent launchers from fitting the same rows.
    lock = root / 'runner.lock'
    fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    os.write(fd, str(os.getpid()).encode()); os.close(fd)
    print(f"{len(jobs) - len(pending)} complete; {len(pending)} remaining", flush=True)

    def run(job):
        if shutil.disk_usage(root).free < 8 * 1024**3:
            raise RuntimeError('Less than 8 GiB free; stop scheduling and retain existing results.')
        if digest(root / 'inputs' / (job['input'] + '.rds')) != job['input_md5']:
            raise RuntimeError('Changed training input: ' + job['job'])
        output = root / 'jobs' / job['job']
        log = root / 'logs' / (job['job'] + '.log')
        return_code = 0
        if not complete(job, root):
            if (output / 'result.rds').exists():
                raise RuntimeError('Result exists without verified status; inspect before resuming: ' + job['job'])
            print('FIT ' + job['job'], flush=True)
            with log.open('a') as stream:
                result = subprocess.run([str(args.launcher), str(code / 'extension' / 'fit-job.R'),
                    str(root), str(code), job['job_number']], stdout=stream, stderr=subprocess.STDOUT)
                return_code = result.returncode
        if args.score and complete(job, root):
            status = json.loads((output / 'status.json').read_text())
            if status['ok'] and not (output / 'score.rds').exists():
                print('SCORE ' + job['job'], flush=True)
                with log.open('a') as stream:
                    scored = subprocess.run([str(args.launcher), str(code / 'extension' / 'score-job.R'),
                        str(root), str(code), job['job_number']], stdout=stream, stderr=subprocess.STDOUT)
                if scored.returncode:
                    raise RuntimeError('Scoring failed; inspect ' + str(log))
        valid = complete(job, root)
        print(f"DONE {job['job']} exit={return_code} recorded={valid}", flush=True)
        return dict(job=job['job'], exit=return_code, recorded=valid)

    try:
        # Keep only one task per worker submitted, so a failure does not launch the full queue.
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            active = {}
            iterator = iter(pending)
            for _ in range(args.workers):
                job = next(iterator, None)
                if job is not None: active[pool.submit(run, job)] = job
            while active:
                done, _ = concurrent.futures.wait(active, return_when=concurrent.futures.FIRST_COMPLETED)
                for future in done:
                    future.result()
                    del active[future]
                    job = next(iterator, None)
                    if job is not None: active[pool.submit(run, job)] = job
        print('All requested jobs have finished. Fit failures and provisional diagnostics remain in their records.', flush=True)
    finally:
        lock.unlink(missing_ok=True)


if __name__ == '__main__':
    main()
