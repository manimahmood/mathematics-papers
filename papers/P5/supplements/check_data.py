"""Read-only checks of the explicitly stated finite-data relationships."""
from collections import Counter
from fractions import Fraction as Q
import json
from math import isqrt
from pathlib import Path


def need(ok, name):
    if not ok:
        raise ValueError(name)


def sqrt_bounds(q):
    need(q >= 0, 'NEGATIVE_RADICAND')
    scale = 10**80
    k = isqrt(q.numerator*scale*scale//q.denominator)
    lo, hi = Q(k, scale), Q(k+1, scale)
    need(lo*lo <= q < hi*hi, 'SQRT_ENCLOSURE')
    return lo, hi


def main():
    d = json.loads(Path(__file__).with_name('triad_data.json').read_bytes())
    need(d['schema'] == 'coexact-triads.public-data.v1', 'SCHEMA')
    c = d['tower_census']
    need(len(c['census']) == c['type_triples_classified'] == 816, 'CENSUS_COUNT')
    need(len({row[0] for row in c['census']}) == 816, 'CENSUS_UNIQUENESS')
    counts = Counter(band+':'+label for _, band, label in c['census'])
    need(dict(counts) == c['classification_counts'], 'CLASSIFICATION_ACCOUNTING')
    need(counts['CROSS:NONZERO'] == 171 and counts['H2:NONZERO'] == 58, 'NONZERO_TYPES')
    need(len(c['entries']) == 658, 'REDUCED_ENTRY_COUNT')
    nonzero = Counter()
    for e in c['entries']:
        sq = Q(e['reduced_sq'])
        need(sq == Q(e['lambda2'])*Q(e['threejL_sq']), 'REDUCED_PRODUCT')
        if sq:
            nonzero[e['band']] += 1
    need(nonzero['CROSS'] == c['cross_nonzero_reduced_entries'] == 385, 'CROSS_COUNT')
    need(nonzero['H2'] == c['h2_nonzero_reduced_entries'] == 160, 'HIGH_COUNT')
    need(Q(180,12467)**2 == Q(28503900,237409081)*Q(556,320229), 'ANCHOR_NORMALIZATION')
    for name, r in d['band_ratio_framing']['ratios'].items():
        lo, hi = map(Q, r['amplitude_ratio_enclosure'])
        need(0 <= lo < hi and lo*lo < Q(r['reduced_sq_ratio_exact']) < hi*hi, 'RATIO:'+name)
    windows = {
        'SUM0': (28, '-0.215211944909057098265281950100', '-0.215211944909057098265281950099',
                 '0.000000172881536303333133674600', '0.000000172881536303333133674601'),
        'SUMX': (149, '40.310908737094023421275640522452', '40.310908737094023421275640522453',
                 '0.444138023269106169431820390834', '0.444138023269106169431820390835'),
        'SUMH2': (73, '15.446453832334867491005362912980', '15.446453832334867491005362912981',
                 '-0.074211304316484249049721645542', '-0.074211304316484249049721645541'),
    }
    for name, (n, *bounds) in windows.items():
        terms = d['complete_coherent_sum_term_lists'][name]
        need(len(terms) == n, 'SUM_TERMS:'+name)
        real = [Q(0), Q(0)]
        imag = [Q(0), Q(0)]
        for rad, a, b in terms:
            lo, hi = sqrt_bounds(Q(rad))
            for coeff, total in ((a, real), (b, imag)):
                need(type(coeff) is int, 'INTEGER_COEFFICIENT')
                x, y = sorted((coeff*lo, coeff*hi))
                total[0] += x
                total[1] += y
        rlo, rhi, ilo, ihi = map(Q, bounds)
        need(rlo <= real[0] <= real[1] <= rhi and ilo <= imag[0] <= imag[1] <= ihi, 'SUM_WINDOW:'+name)
    print(json.dumps(dict(result='FINITE_TRIAD_DATA_RELATIONSHIPS_PASS',
        census_rows=816, product_rows=658, coherent_terms=250,
        independent_mode_integral_derivation=False), sort_keys=True))


if __name__ == '__main__':
    main()
