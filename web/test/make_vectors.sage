# Makes web/data/knots.json (KnotInfo/LinkInfo names -> braid words) and web/test/vectors.json (Sage's Seifert
# matrix, Alexander polynomial, signature and component count for a list of braids, to test seifert.js against).
#   sage web/test/make_vectors.sage
import json, random, os
from sage.knots.link import Link
from sage.knots.knotinfo import KnotInfo
here = os.path.dirname(os.path.abspath(__file__))
R = LaurentPolynomialRing(ZZ, 't'); t = R.gen()

def alex_coeffs(p):                       # Laurent polynomial -> [min degree, [coefficients]]
    p = R(p)
    if p == 0: return [int(0), [int(0)]]
    d = p.valuation(); q = p.shift(-d).polynomial_construction()[0] if hasattr(p, 'polynomial_construction') else None
    coeffs = [int(p[i]) for i in range(d, p.degree() + 1)]
    return [int(d), coeffs]

def sage_data(word):
    # Sage's Link freely reduces the braid word and drops the strands no crossing touches, so the word it actually
    # uses (L.braid().Tietze()) is recorded as 'braid', the given one as 'input'
    n = max(abs(x) for x in word) + 1
    L = Link(BraidGroup(n)(list(word)))
    V = L.seifert_matrix()
    return { 'input': [int(x) for x in word], 'braid': [int(x) for x in L.braid().Tietze()], 'components': int(L.number_of_components()),
             'seifert': [[int(v) for v in row] for row in V.rows()], 'alexander': alex_coeffs(L.alexander_polynomial()),
             'signature': int(L.signature()), 'sage_genus': int(L.genus()) }

# ---- the table: every KnotInfo knot through 12 crossings and every LinkInfo link through 11 crossings
table = []
for K in KnotInfo:
    try:
        cr = int(K.crossing_number())
        if K.is_knot() and cr > 12: continue
        if not K.is_knot() and cr > 11: continue
        b = [int(x) for x in K.braid_notation()]
    except Exception as e:
        continue
    rec = { 'name': K.name, 'braid': b, 'crossings': cr, 'knot': bool(K.is_knot()) }
    if K.is_knot():
        try: rec['genus'] = int(K.three_genus())
        except Exception: pass
    table.append(rec)
json.dump(table, open(os.path.join(here, '..', 'data', 'knots.json'), 'w'), separators=(',', ':'))
print('table entries', len(table))

# ---- the vectors: named examples, random braids, and the knots through 9 crossings
examples = { 'trefoil': [1,1,1], 'mirror trefoil': [-1,-1,-1], 'figure-eight': [1,-2,1,-2], 'hopf': [1,1], 'hopf mirror': [-1,-1],
             '5_1': [1,1,1,1,1], '6_1': [1,1,2,-1,-3,2,-3], '8_19 = T(3,4)': [1,2]*4, '10_124 = T(3,5)': [1,2]*5,
             '12n_242 = P(-2,3,7)': [1,2,2,1,1,2,2,2,2,2,2,2], 'whitehead': [1,1,2,-1,-1,2], 'borromean': [1,-2,1,-2,1,-2],
             'T(2,5)': [1]*5, 'T(4,5)': [1,2,3]*5, 'unknot 2 strands': [1], 'split': [1,3,1,3], '4 strands': [1,2,3,-1,2,-3,1,2] }
vectors = { 'examples': [] }
for name, w in examples.items():
    try:
        d = sage_data(w); d['name'] = name; vectors['examples'].append(d)
    except Exception as e:
        print('skip', name, e)
random.seed(int(7))
vectors['random'] = []
while len(vectors['random']) < 160:
    n = int(random.choice([2, 3, 3, 4, 4, 5])); c = int(random.randint(2, 14))
    w = [int(random.choice([1, -1])) * int(random.randint(1, n - 1)) for _ in range(c)]
    if set(abs(x) for x in w) != set(range(1, n)): continue        # connected Bennequin surface
    try: vectors['random'].append(sage_data(w))
    except Exception as e: print('skip', w, e)
vectors['knots'] = []
for rec in table:
    if rec['knot'] and rec['crossings'] <= 9:
        d = sage_data(rec['braid']); d['name'] = rec['name']; d['knotinfo_genus'] = rec.get('genus'); vectors['knots'].append(d)
json.dump(vectors, open(os.path.join(here, 'vectors.json'), 'w'), separators=(',', ':'))
print('vectors', len(vectors['examples']), len(vectors['random']), len(vectors['knots']))
