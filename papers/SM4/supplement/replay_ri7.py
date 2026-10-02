"""Recheck the reduced stored-array study, not a replay of raw time captures."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    out=a.out.resolve()
    if out==HERE or HERE in out.parents:raise SystemExit("OUTPUT_MUST_BE_EXTERNAL")
    out.mkdir(parents=True,exist_ok=False)
    data=HERE/"data"
    def run(calibration=data/"CALIBRATION.json",prediction=data/"PREDICTION.json",
            targets=data/"TARGETS.json",owner=data/"OWNER.json"):
        cmd=[sys.executable,"-B",str(HERE/"ri7/verify.py")]
        for k,v in (("calibration",calibration),("prediction",prediction),("targets",targets),("owner",owner)):
            cmd += ["--"+k,str(v)]
        p=subprocess.run(cmd,capture_output=True,timeout=90,cwd=out)
        return p,json.loads(p.stdout)
    p,result=run()
    (out/"baseline.stdout.json").write_bytes(p.stdout)
    (out/"baseline.stderr.txt").write_bytes(p.stderr)
    checks={"baseline_pass":p.returncode==0 and result["result"]=="RI7_PREDICTIONS_PASS",
            "forecasts_16":result["counts"]["forecasts"]==16,
            "refusals_10":result["counts"]["refusals"]==10,
            "maps_108":result["counts"]["observed_role_maps"]==108}
    mutants={}
    for name,label in (("missing_calibration","OBSERVATION_COUNT"),
                       ("zero_sigma","SIGMA_INVALID"),("tampered_forecast","PREDICTION_MISMATCH")):
        if name=="tampered_forecast":
            value=json.loads((data/"PREDICTION.json").read_bytes())
            value["rows"][0]["targets"][0]["Rmeasured"]="1"
            kw="prediction"
        else:
            value=json.loads((data/"CALIBRATION.json").read_bytes())
            if name=="missing_calibration":value["rows"].pop()
            else:value["sigma"]="0"
            kw="calibration"
        f=out/(name+".json")
        with f.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(value,stream,sort_keys=True,indent=1);stream.write("\n")
        pr,rr=run(**{kw:f})
        (out/(name+".stdout.json")).write_bytes(pr.stdout)
        (out/(name+".stderr.txt")).write_bytes(pr.stderr)
        ok=pr.returncode==1 and rr.get("owner")==label
        checks[name+"_killed"]=ok
        mutants[name]={"exit":pr.returncode,"owner":rr.get("owner"),"expected_owner":label}
    metrics=result["rows"]
    vals={"forecast_error_ppm":[min(r["measured_relative_error"] for r in metrics)*1e6,
                                max(r["measured_relative_error"] for r in metrics)*1e6],
          "radius_percent":[min(r["relative_radius"] for r in metrics)*100,
                             max(r["relative_radius"] for r in metrics)*100],
          "same_data_comparator_max":max(r["comparator_relative_error"] for r in metrics)}
    report={"result":"SM4_REDUCED_ARRAY_REPLAY_PASS" if all(checks.values()) else "SM4_REDUCED_ARRAY_REPLAY_FAIL",
            "checks":checks,"values":vals,"mutations":mutants,
            "scope":"Reduced response arrays and bound reconstruction, not raw-capture execution, blindness re-establishment, or physical data."}
    with (out/"REPLAY.json").open("x",encoding="utf-8",newline="\n") as f:
        json.dump(report,f,sort_keys=True,indent=1);f.write("\n")
    print(json.dumps(report,indent=1))
    return 0 if all(checks.values()) else 1
if __name__=="__main__":
    raise SystemExit(main())
