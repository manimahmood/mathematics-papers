"""Fresh replay of the four known interval calculations, never a new blind test.

Usage: python -B replay.py --out NEW_EXTERNAL_DIRECTORY
Requires mpmath 1.3.0. Only the new output directory is written.
Arithmetic is the retained implementation; PF99 additionally stores successes.
Near-zero continuum charge and weighted-tail proof are in PROOF.md.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from join_checks import sign_checks, dependency_checks

HERE = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    out = a.out.resolve()
    if out == HERE or HERE in out.parents:
        raise SystemExit("OUTPUT_MUST_BE_EXTERNAL")
    out.mkdir(parents=True, exist_ok=False)
    (out / "evidence").mkdir()
    shutil.copytree(HERE / "src", out / "src")
    for f in (HERE / "evidence").iterdir():
        if f.is_file():
            shutil.copy2(f, out / "evidence" / f.name)
    started = time.monotonic()

    def run(name):
        with (out / (name + ".stdout.txt")).open("xb") as so, (out / (name + ".stderr.txt")).open("xb") as se:
            t = time.monotonic()
            pr = subprocess.run([sys.executable, "-B", str(out / "src" / name), "REPLAY"],
                                cwd=out, stdout=so, stderr=se, timeout=1200)
        return {"source": name, "exit": pr.returncode, "seconds": time.monotonic() - t}

    def sequence(names):
        result = []
        for name in names:
            row = run(name)
            result.append(row)
            if row["exit"] != 0: raise RuntimeError("SOURCE_PROCESS_REFUSED:"+name)
            if name.startswith("pf96"):
                fresh = json.loads((out/"evidence/PF96B_REPLAY.json").read_bytes())
                if fresh["verdict"]!="ALL_PASS" or not all(sign_checks(fresh,"ii").values()):
                    raise RuntimeError("FRESH_PF96_JOIN_REFUSED")
        return result

    # Two workers; the dependency is now PF96 -> PF97, never retained -> PF97.
    with ThreadPoolExecutor(max_workers=2) as pool:
        groups = list(pool.map(sequence, (("pf96b_certified_late_sign.py","pf97b_certified_unique_zero.py"),
                                         ("pf98b_strong_geometry_certified.py","pf99b_initial_stretch_certified.py"))))
    runs = [row for group in groups for row in group]
    data = {n: json.loads((out / "evidence" / ("PF" + n + "B_REPLAY.json")).read_text())
            for n in ("96", "97", "98", "99")}
    checks = {"processes_exit_0": all(x["exit"] == 0 for x in runs),
              "four_retained_gate_predicates": all(x["verdict"] == "ALL_PASS" for x in data.values()),
              "slow_zero_bracket": data["97"]["zero_bracket"] == ["2199/128", "10997/640"],
              "strong_zero_bracket": data["98"]["zero_bracket_sub"] == ["1307/2400", "109/200"],
              "slow_slabs": data["96"]["steps"] == 1600,
              "strong_slabs": data["98"]["steps"] == 1200}
    checks.update(sign_checks(data["96"], "ii"))
    checks.update(sign_checks(data["98"], "i"))
    checks.update(dependency_checks(data["97"], (out/"evidence/PF96B_REPLAY.json").read_bytes()))
    for key, start, end, sign, count, c, safe in (
            ("ii", F(3,40), F(19,40), 1, 2752, F(1,4), F(3,100)),
            ("i", F(1,300), F(13,600), -1, 1472, F(8), F(299,10))):
        part = data["99"]["parts"][key]
        rows = part["all_subcells"]
        checks[key + "_complete_grid"] = (len(rows) == count and
            F(rows[0]["t0"]) == start and F(rows[-1]["t1"]) == end and
            all(F(x["t0"]) < F(x["t1"]) and x["sgnC"] == sign for x in rows) and
            all(F(x["t1"]) == F(y["t0"]) for x,y in zip(rows, rows[1:])))
        # Reporting check against a deliberately widened binary64 endpoint.
        # The independent exact coefficient proof is supplied in PROOF.md.
        from math import factorial, nextafter
        early = next(x for x in part["early_rows"] if F(x["tau1"]) == start)
        low = F.from_float(nextafter(early["bracket_lo"], float("-inf")))
        charge = 4 * (2*c)**24 * start**22 / (factorial(24)*(1-2*c*start))
        checks[key + "_continuum_charge_keeps_safe_coefficient"] = low - charge > safe
    report = {"result": "SM1_FRESH_SOLVER_REPLAY_PASS" if all(checks.values()) else "SM1_FRESH_SOLVER_REPLAY_FAIL",
              "checks": checks, "runs": runs, "seconds": time.monotonic() - started,
              "scope": "Same-code known-result replay, not independent arithmetic or a new prediction.",
              "sources": {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in (HERE / "src").glob("*.py")},
              "driver_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "join_checker_sha256": hashlib.sha256((HERE/"join_checks.py").read_bytes()).hexdigest()}
    with (out / "REPLAY.json").open("x", encoding="utf-8", newline="\n") as f:
        json.dump(report, f, indent=1, sort_keys=True)
        f.write("\n")
    print(json.dumps(report, indent=1))
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
