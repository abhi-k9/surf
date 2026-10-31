"""Exact 8D = 7 positive modes + 1 necessarily signed first-normal edge mode.
The signed complement is an analytic jet mode, NOT a positive geometry control.
"""
from pathlib import Path
import json,pickle,sympy as s
D=Path(__file__).resolve().parent/'v12_reference'
cert=json.loads((D/'exact_edge_certificate.json').read_text())
a=pickle.loads((D/'positive_experiment.pkl').read_bytes())
BB=a['BB'];L=a['q1'];Q=s.Matrix([[s.Rational(x) for x in row] for row in cert['coefficient_matrix_8x7']])
assert Q.shape==(8,7) and Q.rank()==7
assert L.rank()==1 and L*Q==s.zeros(1,7)
j=next(j for j in range(8) if L[0,j]!=0)
h=s.zeros(8,1);h[j]=1/L[0,j]
C=s.Matrix.hstack(Q,h)
assert C.det()!=0 and (L*h)[0]==1
signed_BB=BB*h
assert min(signed_BB)<0 and max(signed_BB)>0
result={'compact_two_macro_dim':8,'positive_kernel_dim':7,'first_normal_jet_index':j,'unit_first_normal_signed_mode_coefficient':str(h[j]),'exact_direct_sum':True,'signed_mode_min_bernstein':str(min(signed_BB)),'signed_mode_max_bernstein':str(max(signed_BB)),'all_positive_modes_first_normal_zero':True}
(Path(__file__).resolve().parent/'signed_edge_results.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS direct sum rank 8 and unit first-normal signed mode, with both Bernstein signs')
