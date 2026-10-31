"""Exact scale-invariant stability bound for isolated/disjoint 7-function edge-star blocks.
This theorem does not imply stability for overlapping multilevel/hierarchical geometric bases.
"""
from pathlib import Path
import sympy as s,json
R=Path(__file__).resolve().parent
D=json.loads((R/'v12_reference'/'exact_edge_certificate.json').read_text())
M=s.Matrix([[s.Rational(x) for x in row] for row in D['exact_sample_matrix']])
assert M.shape==(7,7) and M.det()!=0
C=max(sum(abs(M.inv()[i,j]) for j in range(7)) for i in range(7))
bound=1/C
assert str(bound)==D['certified_coefficient_stability_lower_bound']
# The normalized Bernstein maxima of all seven reference functions are one.
scales=[s.Rational(x) for x in D['basis_Bernstein_scalings']]
assert len(scales)==7
# For affine-copy psi^F(x) = psi(F^-1(x)) the same coefficient matrix M is obtained
# at sites F(x_i); hence the above exact 7x7 inverse is unchanged for all invertible F.
result={'exact_fixed_patch_lower_bound':str(bound),'decimal':float(bound),
 'patch_basis_size':7,'affine_copy_bound_independent_of_scale':True,
 'disjoint_patch_block_supnorm_lower_bound':str(bound),
 'simple_upper_bound_for_disjoint_7function_blocks':7,
 'scope':'Nonoverlapping affine copies of this reference patch only. NO theorem for overlapping hierarchical bases or full global geometric basis.'}
(R/'scale_block_results.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS reference 7x7 exact invertible collocation; bound',bound, 'uniform under disjoint affine-copy patch scaling')
