"""Independent physical-coordinate seam C0/C1/C2 pointwise polynomial certification.
No matching of local normal coordinates without the inverse physical Jacobian.
"""
import json
from pathlib import Path
import sympy as s
import physical_nodal as X
from authentic_ps12 import deriv_row
from fan_assembly import points5,points9
Q=s.Rational;P=Path(__file__).resolve().parent
C=json.loads((P/'positive_physical_coefficients.json').read_text())

def checked_fan(n,pts):
 # 42-DOF physical element polynomial reference; matrices (252x42) per sector
 raw=[]
 for i in range(n):raw.append(X.L*X.Ninv*X.nodal_map(pts[i],pts[(i+1)%n]).inv())
 center=s.Matrix([[Q(z) for z in row] for row in C[f'fan{n}']['center_coefficients']]);count=0;cases_count=0
 for i in range(n):
  prev=(i-1)%n;post=i;edge=s.Matrix([[Q(z) for z in row] for row in C[f'fan{n}']['edge_coefficients_by_radial_index'][i]])
  cases=[('center-3',raw[prev][:,:10]*center,raw[post][:,:10]*center),
         ('central-value',raw[prev][:,0],raw[post][:,0]),
         ('outer-value',raw[prev][:,20],raw[post][:,10]),
         ('edge-3',raw[prev][:,38:42]*s.diag(-1,1,-1,-1)*edge,raw[post][:,30:34]*edge)]
  pi=s.Matrix(list(map(Q,pts[i])));normal=s.Matrix([pi[1],-pi[0]])
  def invj(p0,p1):
   p0=s.Matrix(list(map(Q,p0)));p1=s.Matrix(list(map(Q,p1)))
   return s.Matrix.hstack((p1-p0)/3,(-p0)/3).inv()*normal
  d1=invj(pts[prev],pts[i]);d2=invj(pts[i],pts[(i+1)%n]);cases_count+=len(cases)
  for _name,P0,P1 in cases:
   for a,b in ((Q(0),Q(1,2)),(Q(1,2),Q(1))):
    for k in range(3):
     for j in range(6-k):
      t=a+(b-a)*Q(j+1,7-k)
      uv0=(3*t,3*(1-t));uv1=(Q(0),3*(1-t));c0=X.cell_for(uv0);c1=X.cell_for(uv1)
      r0=s.Matrix([deriv_row(uv0[0],uv0[1],d1[0],d1[1],k)])*P0[21*c0:21*(c0+1),:]
      r1=s.Matrix([deriv_row(uv1[0],uv1[1],d2[0],d2[1],k)])*P1[21*c1:21*(c1+1),:]
      assert r0==r1,(n,i,_name,k,t,r0-r1)
      count+=r0.cols
 return {'radial_edges':n,'families_per_edge':4,'scalar_exact_value_normal_1_2_tests':count,'verified_C2':True}

if __name__=='__main__':
 out={}
 for n,pts in ((5,points5),(9,points9)):
  v=checked_fan(n,pts);out[f'fan{n}']=v
  print('fan',n,'exact physical seam scalar equalities',v['scalar_exact_value_normal_1_2_tests'],flush=True)
 (P/'physical_C2_seam_results.json').write_text(json.dumps(out,indent=2))
