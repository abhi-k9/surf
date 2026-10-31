"""Exact proof that the authentic 42D PS12 reference element reproduces all P5.
Uses local polynomial jet functionals and independent Bernstein conversion of 21 monomials.
"""
from pathlib import Path
import sys,math,json
import sympy as s
root=Path(__file__).resolve().parent
sys.path.insert(0,str(root/'base_v17'))
import physical_nodal as p
from authentic_ps12 import cells,V,monos
Q=s.Rational
assert p.N.shape==(42,42) and p.N*p.Ninv==s.eye(42)

# Direct polynomial trace/jet construction from each reference monomial.
def nodals_for_monomial(index):
 out=[]
 for pt in p.verts:
  for a,b in p.jetkeys:
   out.append(p.mon_jet(pt,[(1,0)]*a+[(0,1)]*b)[index])
 for i,j in ((0,1),(1,2),(2,0)):
  A=p.verts[i];B=p.verts[j];mid=tuple((A[k]+B[k])/2 for k in (0,1));t=(B[0]-A[0],B[1]-A[1]);n=(t[1],-t[0])
  for ds in ([n],[n,n],[t,n,n],[n,n,n]):out.append(p.mon_jet(mid,ds)[index])
 return s.Matrix(out)

checks=[]
for j,(a,b) in enumerate(monos):
 # Reconstruct from real D-Y nodal dual and compare with Bernstein coefficients calculated independently.
 q=p.BREF*nodals_for_monomial(j)
 expected=[]
 for cell in cells:
  tri=[V[z] for z in cell]
  values=[]
  for l in p.pts:
   x=sum(l[k]*tri[k][0] for k in range(3));y=sum(l[k]*tri[k][1] for k in range(3))
   values.append(x**a*y**b)
  expected.extend(p.BMi*s.Matrix(values))
 assert q==s.Matrix(expected),(j,a,b)
 checks.append((a,b))
result={'space':'Authentic 42D local Davydov–Yeo quintic PS12 element', 'tested_monomials':[[a,b] for a,b in checks],'monomials_count':len(checks),'exact_P5_reproduction':True,'max_degree':5,'global_scope':'Local polynomial reproduction; global reproduction follows when physical Hermite data are assembled consistently. This is not reproduction by restricted positive geometric controls alone.'}
(root/'p5_analytic_results.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS exact 21 monomials 0<=a+b<=5 in full 42D PS12 reference element')
