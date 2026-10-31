"""Deterministic exact rational regression checks for committed NPMS v1.2 certificates.
Requires SymPy; uses committed rational coefficients (no LP or candidate search).
NOTE: not a general arbitrary-topology spline proof.
"""
from pathlib import Path
import pickle,json
import sympy as S
from authentic_ps12 import rows
P=Path(__file__).parent
load=lambda fn:json.loads((P/fn).read_text())
e=load('exact_edge_certificate.json');t=load('exact_refinement_certificate.json');Q=S.Matrix([[S.Rational(x) for x in row] for row in e['coefficient_matrix_8x7']]);sc=S.diag(*[S.Rational(x) for x in e['basis_Bernstein_scalings']]);B=pickle.loads((P/'positive_experiment.pkl').read_bytes())['BB']
assert Q.shape==(8,7) and Q.rank()==7
F=B*Q*sc
assert F.shape==(504,7) and min(F)>=0 and all(max(F[:,i])==1 for i in range(7))
original_rows,_=rows();L=S.Matrix([[S.Rational(x) for x in r.split()] for r in (P/'full_local_42.txt').read_text().splitlines()]).T
assert L.shape==(252,42)
for row in original_rows:
 for j in range(42):
  assert sum(S.Rational(v.numerator,v.denominator)*L[k,j] for k,v in row.items())==0
obj=pickle.loads((P/'edge_exact.pkl').read_bytes());BL=obj['BL'];SS=obj['SS'];Z=obj['Z'];assert BL.shape==(84,8) and SS*Z==S.zeros(SS.rows,8)
q1=pickle.loads((P/'positive_experiment.pkl').read_bytes())['q1'];q2=pickle.loads((P/'positive_experiment.pkl').read_bytes())['q2'];assert q1.rank()==1 and q1*Q==S.zeros(1,7)
assert min(B*Q)>=0
M7=S.Matrix([[S.Rational(z) for z in r] for r in e['exact_sample_matrix']]);assert M7.rank()==7
bound7=1/max(sum(abs(M7.inv()[i,j]) for j in range(7)) for i in range(7))
assert str(bound7)==e['certified_coefficient_stability_lower_bound']
X=S.Matrix([[S.Rational(x) for x in r] for r in t['coarse_to_fine_transfer_10x7']]);assert X.shape==(10,7) and min(X)>=0 and X.rank()==7
assert t['coarse_basis_size']==7 and t['fine_basis_size']==10
assert S.Rational(t['corrected_affine_min_bernstein'])>0
assert S.Rational(t['fixed_patch_fine_coefficient_stability_bound'])>0
print('PASS exact 42D local continuity, 8D edge-star, 7 positive independent modes, 3,528 nonnegative Bernstein tests, seam-jet obstruction, and certified collocation bound')
print('PASS exact 10x7 nonnegative transfer and strictly positive affine-remainder certificate')
print('Bounds:',float(bound7),t['bound_float'])
