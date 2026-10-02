"""Read-only rational regression of the conditional adaptive-pivot bound.

Copyright (c) 2026 Luqman Mahmood. Software rights reserved.
See the package LICENSE.md. The manuscript contains the general proof.
"""

from fractions import Fraction as F
from itertools import product
import json
import sys


def main():
    if len(sys.argv) != 1:
        raise SystemExit("usage: python -B supplements/check_adaptive_pivot.py")
    directions = [tuple(map(F, w)) for w in product(range(-2, 3), repeat=3)
                  if any(w)]
    errors = list(product(map(F, (-1, 0, 1)), repeat=2))
    checks = 0
    sharp = False
    for w in directions:
        # Normalizing w cancels in the reconstruction map; no square root needed.
        p = max(range(3), key=lambda i: abs(w[i]))
        assert p == min(i for i in range(3) if abs(w[i]) == max(map(abs, w)))
        factor_squared = sum(x*x for x in w) / (w[p]*w[p])
        assert factor_squared <= 3
        q = [i for i in range(3) if i != p]
        for e in errors:
            n2 = sum(x*x for x in e)
            pivot_error = -sum(w[i]*e[j] for j, i in enumerate(q))/w[p]
            recovered_n2 = n2 + pivot_error*pivot_error
            assert recovered_n2 <= factor_squared*n2 <= 3*n2
            sharp |= n2 > 0 and recovered_n2 == 3*n2
            checks += 1
    # The earlier first-nonzero rule reconstructs exactly but has no such bound.
    w = (F(1, 10**6), F(1), F(0))
    old_p, new_p = 0, max(range(3), key=lambda i: abs(w[i]))
    old_factor = sum(x*x for x in w)/(w[old_p]*w[old_p])
    new_factor = sum(x*x for x in w)/(w[new_p]*w[new_p])
    assert old_factor > 3 and new_factor <= 3 and sharp
    assert checks == 1116
    print(json.dumps({
        "result": "P6_ADAPTIVE_PIVOT_REGRESSION_PASS",
        "directions": len(directions), "error_rows": checks,
        "sharp_bound_attained": sharp,
        "first_nonzero_counterexample_factor_squared": str(old_factor),
        "maximum_pivot_factor_squared": str(new_factor),
        "scope": "exact first column; second-column recovery only",
    }, sort_keys=True))


if __name__ == "__main__":
    main()
