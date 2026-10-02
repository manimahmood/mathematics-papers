"""Exact grid/summary joins on replayed enclosures; not a new interval solver."""
from fractions import Fraction as F
import hashlib
import math

SPECS = {
    "ii": (1600, F(1, 40), [(0,"0","19/40"),(1,"19/40","137/8"),
                              (0,"137/8","689/40"),(-1,"689/40","40")]),
    "i": (1200, F(1, 600), [(0,"0","13/600"),(-1,"13/600","163/300"),
                              (0,"163/300","41/75"),(1,"41/75","2")]),
}

def sign_checks(data, label):
    count, h, spec = SPECS[label]
    expected = [{"sign":s,"from":a,"to":b} for s,a,b in spec]
    checks = {label+"_step_grid":False, label+"_sign_runs":False,
              label+"_straddle_rows":False, label+"_initial_join":False}
    if label == "i": checks["i_zero_join"] = False
    try:
        rows = data["steps_lohi"]
        grid = len(rows)==count and data["steps"]==count
        runs, straddles = [], []
        for k, row in enumerate(rows):
            t, lo, hi = row
            grid = grid and F(t)==k*h and math.isfinite(lo) and math.isfinite(hi) and lo<=hi
            s = 1 if lo>0 else (-1 if hi<0 else 0)
            start, end = str(k*h), str((k+1)*h)
            if runs and runs[-1]["sign"]==s: runs[-1]["to"] = end
            else: runs.append({"sign":s,"from":start,"to":end})
            if s==0: straddles.append({"t0":start,"t1":end,"lo":lo,"hi":hi})
        checks[label+"_step_grid"] = bool(grid)
        checks[label+"_sign_runs"] = runs==expected and data["sign_runs"]==expected
        checks[label+"_straddle_rows"] = data["straddle_steps"]==straddles
        early = F(spec[0][2])
        late0, late1 = F(spec[2][1]), F(spec[2][2])
        wanted = [(str(k*h), str((k+1)*h)) for k in range(count)
                  if k*h<early or late0<=k*h<late1]
        checks[label+"_initial_join"] = (data["early_straddle_end"]==str(early)
            and [(x["t0"],x["t1"]) for x in straddles]==wanted)
        if label=="i": checks["i_zero_join"] = data["zero_run"]==[str(late0),str(late1)]
    except (KeyError, TypeError, ValueError, OverflowError):
        pass
    return checks

def dependency_checks(pf97, fresh_pf96_bytes):
    binding = pf97.get("pf96_input", {})
    return {"pf97_fresh_pf96_filename":binding.get("file")=="PF96B_REPLAY.json",
            "pf97_fresh_pf96_sha256":binding.get("sha256")==hashlib.sha256(fresh_pf96_bytes).hexdigest()}
