# The braid word of a knot or link, for the app's input line, from SnapPy/spherogram (Vogel's algorithm) or Sage:
#   sage web/data/braid.sage K12n242            a KnotInfo/LinkInfo name or a Rolfsen name (3_1, L6a4)
#   sage web/data/braid.sage "[[1,5,2,4],[3,1,4,6],[5,3,6,2]]"    a PD code
#   sage web/data/braid.sage "DT:[4,6,2]"       a DT code
# SnapPy's braid words use the same convention as Sage and the app (σ_i is the positive crossing; the closure of
# σ_1^3 has signature -2), but the chirality of a *named* knot depends on the table's diagram: SnapPy's 3_1 is the
# mirror of KnotInfo's (signature +2), its K12n242 is not.  Negate the word for the mirror image.
import sys, ast
from sage.knots.link import Link
arg = sys.argv[1] if len(sys.argv) > 1 else '3_1'
if arg.startswith('DT:'):
    import snappy
    L = snappy.Link('DT:' + arg[3:])
elif arg.startswith('['):
    import snappy
    L = snappy.Link(ast.literal_eval(arg))
else:
    import snappy
    L = snappy.Link(arg)
word = [int(x) for x in L.braid_word()]
n = max(abs(x) for x in word) + 1
S = Link(BraidGroup(n)(word))
print(' '.join(str(x) for x in word))
print('# %d strands, %d crossings, %d component(s), signature %d, Alexander %s' % (n, len(word), S.number_of_components(), S.signature(), S.alexander_polynomial()))
