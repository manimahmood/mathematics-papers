"""Replay the three existing conditioning verifiers in a fresh external copy."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    here = Path(__file__).resolve().parent
    out = args.out.resolve()
    if out.exists() or out.is_relative_to(here.parent):
        raise ValueError('NEW_EXTERNAL_OUTPUT_REQUIRED')
    shutil.copytree(here/'conditioning', out)
    tasks = [
        ('p6_t1', 'P6_T1_COND0_VERIFY_REPORT_20260813.json'),
        ('p6_t1c', 'P6_T1C_DESIGN0_VERIFY_REPORT_20260813.json'),
        ('p6_t1b', 'P6_T1B_INCUMBENT0_VERIFY_REPORT_20260814.json'),
    ]

    def run(task):
        folder, report = task
        expected = (here/'conditioning'/folder/'results'/report).read_bytes()
        start = time.monotonic()
        p = subprocess.run([sys.executable, '-B', str(out/folder/(folder+'_verify.py'))],
                           cwd=out, capture_output=True, timeout=240)
        (out/(folder+'.stdout.txt')).write_bytes(p.stdout)
        (out/(folder+'.stderr.txt')).write_bytes(p.stderr)
        if p.returncode or p.stderr:
            raise ValueError('REPLAY_FAILED:'+folder)
        actual = (out/folder/'results'/report).read_bytes()
        if actual != expected:
            raise ValueError('REPORT_BYTES_DIFFER:'+folder)
        return dict(name=folder, seconds=time.monotonic()-start,
                    report_sha256=hashlib.sha256(actual).hexdigest(),
                    stdout=p.stdout.decode('ascii'))

    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(run, tasks))
    receipt = dict(result='P6_PORTABLE_CONDITIONING_REPLAY_PASS', rows=results,
                   algorithm_changes=False,
                   verifier_edits='public vocabulary and hash pins updated; arithmetic algorithms unchanged; public reports regenerated')
    (out/'REPLAY.json').write_text(json.dumps(receipt, indent=1, sort_keys=True)+'\n',
                                  encoding='ascii', newline='\n')
    print(json.dumps(receipt, sort_keys=True))


if __name__ == '__main__':
    main()
