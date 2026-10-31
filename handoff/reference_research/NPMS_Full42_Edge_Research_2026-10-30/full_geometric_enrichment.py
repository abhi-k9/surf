"""Experimental stronger positive vertex-Hessian + radial-edge split of v1.7 carrier-free basis.
Native D-Y quintic, exact rational nonnegative Bernstein certificates.
No full P2 reproduction or spatial subdivision is claimed without independent tests.
"""
from pathlib import Path
import sys, math, json, itertools
import numpy as np
import sympy as s
root=Path(__file__).resolve().parent
base=root/'base_v17';sys.path[:0]=[str(base),str(base/'verification')]
from verify_v17 import run,points5,points9
from authentic_ps12 import V,cells
import physical_nodal as X
Q=s.Rational
cert=json.loads((root/'positive_physical_coefficients.json').read_text())

def ratio_abs(p,h):
 if any(p[i]==0 and h[i]!=0 for i in range(p.rows)):return Q(0)
 return min((p[i]/abs(h[i]) for i in range(p.rows) if h[i]!=0),default=Q(0))
def ratio_pos(p,e):
 if any(p[i]==0 and e[i]>0 for i in range(p.rows)):return Q(0)
 return min((p[i]/e[i] for i in range(p.rows) if e[i]>0),default=Q(0))
def best_perm(rs):
 feasible=[p for p in itertools.permutations(range(3)) if all(rs[p[k],k]>0 for k in range(3))]
 assert feasible
 return max(feasible,key=lambda p:min(rs[p[k],k] for k in range(3)))

def target_coefficients(points,deg):
 n=len(points);all_b=[]
 for i in range(n):
  a,b=[s.Matrix([Q(z) for z in z0]) for z0 in (points[i],points[(i+1)%n])]
  for cell in cells:
   tri=[s.Matrix(V[z]) for z in cell]
   vals=[]
   for lam in X.pts:
    z=sum((lam[j]*tri[j] for j in range(3)),s.zeros(2,1))
    pt=a+(b-a)*z[0]/3-a*z[1]/3
    vals.append(pt[0]**deg[0]*pt[1]**deg[1])
   all_b.extend(X.BMi*s.Matrix(vals))
 return s.Matrix(all_b)

def build(n,pts):
 print('BASE',n,flush=True)
 base_result,B=run(pts,'fan'+str(n),True)
 M=[X.physical_bernstein(pts[i],pts[(i+1)%n]) for i in range(n)]
 n0=len(B);m=252*n
 # 3 second-derivative nodal modes per vertex, exact assembled physical-coordinate jets
 H=[]
 for j in (3,4,5):H.append(s.Matrix.vstack(*[A[:,j] for A in M]))
 for i in range(n):
  prev=(i-1)%n
  for j in (3,4,5):
   h=s.zeros(m,1)
   h[252*prev:252*(prev+1),:]=M[prev][:,20+j]
   h[252*i:252*(i+1),:]=M[i][:,10+j]
   H.append(h)
 # center parent assigned one second-jet each
 h_assign={};h_scales={};edge_assign={};edge_scales={}
 rH=s.Matrix([[ratio_abs(B[a],H[k]) for k in range(3)] for a in range(3)])
 ph=best_perm(rH)
 for k in range(3):
  parent=ph[k];h_assign[parent]=k;h_scales[parent]=min(Q(1,2),rH[parent,k]/2)
 print('CENTER Hessian scales',[str(h_scales[k]) for k in range(3)],flush=True)
 E=[]
 for i,rows in enumerate(cert['fan'+str(n)]['edge_coefficients_by_radial_index']):
  coefficient=s.Matrix([[Q(z) for z in row] for row in rows]);prev=(i-1)%n
  col=s.zeros(m,3)
  col[252*prev:252*(prev+1),:]=M[prev][:,38:42]*s.diag(-1,1,-1,-1)*coefficient
  col[252*i:252*(i+1),:]=M[i][:,30:34]*coefficient
  assert min(col)>=0 and coefficient.rank()==3
  E.extend([col[:,k] for k in range(3)])
 for i in range(n):
  parents=[3*(i+1)+k for k in range(3)]
  rsE=s.Matrix([[ratio_pos(B[parents[a]],E[3*i+k]) for k in range(3)] for a in range(3)])
  rsH=s.Matrix([[ratio_abs(B[parents[a]],H[3+3*i+k]) for k in range(3)] for a in range(3)])
  pe=best_perm(rsE);ph=best_perm(rsH)
  for k in range(3):
   ia=parents[pe[k]];ja=parents[ph[k]]
   edge_assign[ia]=3*i+k;edge_scales[ia]=min(Q(1,2),rsE[pe[k],k]/2)
   h_assign[ja]=3+3*i+k;h_scales[ja]=min(Q(1,2),rsH[ph[k],k]/2)
  print('OUTER',i,'E',[str(edge_scales[j]) for j in parents], 'H',[str(h_scales[j]) for j in parents],flush=True)
 # 2 child Hessian polarities for each center parent; 3 children per outer parent incl edge mode
 F=[];parent_for_fine=[]
 for j,Pp in enumerate(B):
  h=H[h_assign[j]];delta=h_scales[j]
  if j<3:
   kids=[(Pp+delta*h)/2,(Pp-delta*h)/2]
  else:
   e=E[edge_assign[j]];eps=edge_scales[j]
   kids=[(Pp-eps*e+delta*h)/2,(Pp-eps*e-delta*h)/2,eps*e]
  assert sum(kids,s.zeros(m,1))==Pp
  for kid in kids:
   assert min(kid)>=0,('positivity',n,j,min(kid))
   F.append(kid);parent_for_fine.append(j)
 n1=len(F);assert n1==9*n+6
 mat=s.Matrix.hstack(*F)
 T=s.zeros(n1,n0)
 for j,pa in enumerate(parent_for_fine):T[j,pa]=1
 assert mat*T==s.Matrix.hstack(*B)
 assert mat*s.ones(n1,1)==s.ones(m,1)
 assert min(mat)>=0
 ar=np.asarray(mat,float);sv=np.linalg.svd(ar/np.linalg.norm(ar,axis=0),compute_uv=False)
 assert np.linalg.matrix_rank(ar)>=n1
 # Same coarse parameter sites copied to each fine child: affine precision inherited exactly.
 # P2 test uses full exact degree-5 Bernstein coefficient vectors as targets.
 quad={}
 for p in ((2,0),(1,1),(0,2)):
  target=target_coefficients(pts,p)
  ar_t=np.asarray(target,float).ravel()
  coeff,_,_,_=np.linalg.lstsq(ar,ar_t,rcond=1e-11)
  residual=float(np.linalg.norm(ar@coeff-ar_t)/max(np.linalg.norm(ar_t),1e-10))
  quad['x^'+str(p[0])+'y^'+str(p[1])]=residual
 res={'valence':n,'coarse_controls':n0,'enriched_controls':n1,'positive':True,'exact_fixed_mesh_transfer':True,'exact_P1_from_parent_copy':True,'support_max_outer_sectors':4,'center_support_sectors':n,'edge_support_sectors':2,'center_hessian_modes':3,'outer_hessian_modes':3*n,'radial_edge_curvature_modes':3*n, 'minimum_outer_edge_scale':str(min(edge_scales.values())), 'minimum_hessian_scale':str(min(h_scales.values())),'column_scaled_bernstein_condition':float(sv[0]/sv[-1]),'quadratic_relative_lstsq_residuals':quad,'quadratic_reproduction_not_certified':True,'geometric_completeness':'Not a full P2/P5 positive basis; this test only adds Hessian degrees of freedom without the full independent edge/vertex jet set.'}
 print('FINISH',n,'controls',n1,'epsEmin',res['minimum_outer_edge_scale'],'epsHmin',res['minimum_hessian_scale'],'cond',res['column_scaled_bernstein_condition'],'P2',quad,flush=True)
 return res
if __name__=='__main__':
 results={}
 for n,pts in ((5,points5),(9,points9)):
  results['fan'+str(n)]=build(n,pts)
 (root/'full_geometric_results.json').write_text(json.dumps(results,indent=2)+'\n')
 print('PASS full positive vertex Hessian and edge-curvature enrichment on both fans',flush=True)
