"""Read-only finite-distribution regression of the two chain-rule losses.

Copyright (c) 2026 Luqman Mahmood. Software rights reserved.
See ../RIGHTS.md. This checks bookkeeping, not a CTMC or a physical experiment.
"""

from decimal import Decimal as D, localcontext
import json
import sys


def probability(rows):
    total = sum(map(sum, rows))
    return [[D(v)/D(total) for v in row] for row in rows]


def kl(p, q):
    return sum((a*(a/b).ln() for a, b in zip(p, q)), D(0))


def check(p, q):
    # theta has two values, X four, and Y=f(X) two: f=(0,0,1,1).
    assert all(sum(row) == D("0.5") for row in p+q)
    px = [p[0][j]+p[1][j] for j in range(4)]
    qx = [q[0][j]+q[1][j] for j in range(4)]
    py = [sum(px[2*k:2*k+2]) for k in range(2)]
    qy = [sum(qx[2*k:2*k+2]) for k in range(2)]
    phase = sum((px[j]*kl([p[i][j]/px[j] for i in range(2)],
                         [q[i][j]/qx[j] for i in range(2)])
                 for j in range(4)), D(0))
    obs = sum((py[k]*kl([px[j]/py[k] for j in range(2*k, 2*k+2)],
                       [qx[j]/qy[k] for j in range(2*k, 2*k+2)])
               for k in range(2)), D(0))
    joint = kl(sum(p, []), sum(q, []))
    observed = kl(py, qy)
    tol = D("1e-60")
    assert abs(joint-observed-phase-obs) < tol
    assert phase >= -tol and obs >= -tol
    return joint, observed, phase, obs


def main():
    if len(sys.argv) != 1:
        raise SystemExit("usage: python -B supplement/check_information_loss.py")
    cases = {
        "both_losses": ([[8, 2, 4, 6], [2, 6, 8, 4]],
                        [[3, 7, 6, 4], [7, 3, 2, 8]]),
        "phase_only": ([[3, 1, 3, 1], [1, 3, 1, 3]], [[2]*4, [2]*4]),
        "observation_only": ([[3, 1, 3, 1]]*2, [[1, 3, 1, 3]]*2),
        "observable_only": ([[3, 3, 1, 1]]*2, [[1, 1, 3, 3]]*2),
        "identical": ([[3, 1, 2, 2]]*2, [[3, 1, 2, 2]]*2),
    }
    with localcontext() as context:
        context.prec = 75
        rows = {name: check(probability(p), probability(q))
                for name, (p, q) in cases.items()}
        tol = D("1e-60")
        assert rows["both_losses"][2] > tol and rows["both_losses"][3] > tol
        assert rows["phase_only"][2] > tol and abs(rows["phase_only"][3]) < tol
        assert rows["observation_only"][3] > tol and abs(rows["observation_only"][2]) < tol
        assert rows["observable_only"][1] > tol
        assert abs(rows["observable_only"][2])+abs(rows["observable_only"][3]) < tol
        assert all(abs(v) < tol for v in rows["identical"])
        # Omitting either strictly positive loss must not reconstruct the joint KL.
        joint, observed, phase, obs = rows["both_losses"]
        assert abs(joint-observed-phase) > tol and abs(joint-observed-obs) > tol
        print(json.dumps({
            "result": "EP_INFORMATION_LOSS_REGRESSION_PASS",
            "cases": len(rows), "decimal_precision": context.prec,
            "omitted_loss_controls_rejected": 2,
            "rows": {name: dict(zip(("joint", "observed", "phase_loss", "path_loss"),
                                     map(str, values))) for name, values in rows.items()},
            "scope": "finite joint distributions; the manuscript proves the general identity",
        }, sort_keys=True))


if __name__ == "__main__":
    main()
