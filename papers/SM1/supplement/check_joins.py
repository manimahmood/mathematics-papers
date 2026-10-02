"""Recheck retained sign summaries and fresh-input binding, without integrating."""
import argparse
import hashlib
import json
from pathlib import Path
from join_checks import sign_checks, dependency_checks

HERE = Path(__file__).resolve().parent

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    out = a.out.resolve()
    if out == HERE or HERE in out.parents:
        raise SystemExit('OUTPUT_MUST_BE_EXTERNAL')
    r = HERE/'retained_replay'
    raw96 = (r/'PF96B_REPLAY.json').read_bytes()
    d96 = json.loads(raw96)
    d97 = json.loads((r/'PF97B_REPLAY.json').read_bytes())
    d98 = json.loads((r/'PF98B_REPLAY.json').read_bytes())
    checks = sign_checks(d96, 'ii') | sign_checks(d98, 'i') | dependency_checks(d97, raw96)
    receipt = json.loads((r/'REPLAY.json').read_bytes())
    checks['full_receipt_21_checks'] = len(receipt['checks']) == 21 and all(receipt['checks'].values())
    checks['full_receipt_pass'] = receipt['result'] == 'SM1_FRESH_SOLVER_REPLAY_PASS'
    sources = {f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (HERE/'src').glob('*.py')}
    checks['source_hashes'] = receipt['sources'] == sources
    for key, name in [('driver_sha256','replay.py'), ('join_checker_sha256','join_checks.py')]:
        checks[key] = receipt[key] == hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    passed = all(checks.values())
    report = {'result':'SM1_JOIN_RECHECK_PASS' if passed else 'SM1_JOIN_RECHECK_FAIL',
              'checks':checks, 'scope':'Retained joins and source binding, not fresh integration.'}
    with out.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(report, f, sort_keys=True, indent=1); f.write('\n')
    print(json.dumps(report, indent=1))
    return 0 if passed else 1

if __name__ == '__main__':
    raise SystemExit(main())
