"""NPMS v1.7 exact local carrier-free affine quintic fan experiment.
Verified on 2 fixed fan geometries only. Does not prove uniform multilevel stability.
"""
import sys,json, math
from pathlib import Path
from functools import reduce
sys.path.insert(0,str(Path(__file__).resolve().parent))
sys.path.insert(0,str(__import__('pathlib').Path(__file__).resolve().parent/'verification'))
import sympy as s
import numpy as np
from moment_jet_probe import build,points5,points9
from authentic_ps12 import V,cells
import physical_nodal as P
Q=s.Rational
ALPHA=Q(1,32);THETA=Q(1,5);L=Q(8)

def run(points,name,return_basis=False):
 n=len(points)
 print('BUILD exact D-Y jets',name,flush=True)
 VV,DX,DY,EX,N=build(points)
 # W_outer_i = (1-alpha) [ (1-2theta)V_i + theta V_i-1 + theta V_i+1 ]
 A=s.zeros(n,n)
 for i in range(n):
  A[i,i]=1-2*THETA
  A[i,(i-1)%n]=THETA
  A[i,(i+1)%n]=THETA
 assert A.det()!=0
 # W0 = V0 + alpha sum outer V; W_outer = (1-alpha) A V_outer
 W=[VV[0]+ALPHA*sum(VV[1:],s.zeros(252*n,1))]
 W.extend((1-ALPHA)*sum((A[i,j]*VV[j+1] for j in range(n)),s.zeros(252*n,1)) for i in range(n))
 assert sum(W,s.zeros(252*n,1))==s.ones(252*n,1)
 # site base for center is center; solve A^T q adjusted by 1-alpha
 p=s.Matrix([[Q(c) for c in t] for t in points]); base=s.zeros(n+1,2)
 q=A.T.inv()*p/(1-ALPHA)
 for i in range(n):base[i+1,0]=q[i,0];base[i+1,1]=q[i,1]
 # moments J_{v,dx}, J_{v,dy}, and first normal edge jets distribute to outer groups
 Hx=[DX[0]]; Hy=[DY[0]]
 for i in range(n):
  x=DX[i+1]+N[i][0]*EX[i]+N[n+i][0]*EX[n+i]
  y=DY[i+1]+N[i][1]*EX[i]+N[n+i][1]*EX[n+i]
  Hx.append(x);Hy.append(y)
 # verify target affine coefficients by exact interpolation on physical cells (avoid redundant interpolation)
 pt=[]
 for k in range(n):
  a,b=map(lambda x:s.Matrix([Q(xx) for xx in x]),(points[k],points[(k+1)%n]))
  for cell in cells:
   verts=[s.Matrix(V[z]) for z in cell]
   for lam in P.pts:
    z=sum((lam[j]*verts[j] for j in range(3)),s.zeros(2,1))
    x=a+(b-a)*z[0]/3-a*z[1]/3
    pt.append((x[0],x[1]))
 tmx=s.Matrix([r[0] for r in pt]);tmy=s.Matrix([r[1] for r in pt]);
 for axis,target,hs in ((0,tmx,Hx),(1,tmy,Hy)):
  val=sum((base[i,axis]*W[i]+hs[i] for i in range(n+1)),s.zeros(252*n,1))
  assert val==target
 print('global moments exact',flush=True)
 # B_i0,1,2 independent, with 3D barycentric locally about q_i; 3 sites q_i+(L,0),(0,L),(-L,-L)
 basis=[];sites=[];min_per_group=[];nzero=0
 for i in range(n+1):
  w,hx,hy=W[i],Hx[i],Hy[i]
  numerators=[L*w+2*hx-hy,L*w+2*hy-hx,L*w-hx-hy]
  mn=min(min(q) for q in numerators)
  assert mn>=0,('positivity failure',name,i,mn)
  Z=[z/(3*L) for z in numerators]
  assert sum(Z,s.zeros(252*n,1))==w
  sx,sy=base[i,0],base[i,1]
  ps=[(sx+L,sy),(sx,sy+L),(sx-L,sy-L)]
  assert sum((ps[j][0]*Z[j] for j in range(3)),s.zeros(252*n,1))==sx*w+hx
  assert sum((ps[j][1]*Z[j] for j in range(3)),s.zeros(252*n,1))==sy*w+hy
  basis.extend(Z);sites.extend(ps);min_per_group.append(str(mn/(3*L)))
  support=[t for t in range(n) if any(z[t*252:(t+1)*252,:]!=s.zeros(252,1) for z in Z)]
  expect=n if i==0 else (4 if n>4 else n)
  assert len(support)<=expect,(name,i,support)
  print(' cert',i,'min',str(mn/(3*L)),'macro sectors',support,flush=True)
 assert sum(basis,s.zeros(252*n,1))==s.ones(252*n,1)
 assert sum((sites[i][0]*basis[i] for i in range(len(basis))),s.zeros(252*n,1))==tmx
 assert sum((sites[i][1]*basis[i] for i in range(len(basis))),s.zeros(252*n,1))==tmy
 # Analytical independence: center jets remove first block; outer derivatives remove 2 DOF per group;
 # remaining group sum funcs are W_i and circulant A invertible at outer vertices.
 assert (A.det()!=0 and L!=0 and ALPHA!=1)
 # Bernstein coefficient sampling (4 points per PS12 microcell, etc)
 B=np.empty((252*n,3*(n+1)))
 for i,fn in enumerate(basis):B[:,i]=np.asarray(fn,float).reshape(-1)
 coeff=B.reshape(n*12,21,-1)
 from math import factorial
 sample_lam=[(1/3,1/3,1/3),(.5,.25,.25),(.25,.5,.25),(.25,.25,.5)]
 weights=np.array([[float(factorial(5)/math.prod(factorial(t) for t in ijk))*np.prod([p[k]**ijk[k] for k in range(3)]) for ijk in P.ijk] for p in sample_lam])
 sampled=np.einsum('ab,kbc->kac',weights,coeff).reshape(-1,3*(n+1))
 norms=np.linalg.norm(sampled,axis=0)
 sv=np.linalg.svd(sampled/norms,compute_uv=False)
 rank=int(np.linalg.matrix_rank(sampled/norms))
 assert rank==3*(n+1)
 result={'fan':name,'valence':n,'basis_functions':len(basis),'outer_support_max_sectors':4,'center_support_sectors':n,'external_affine_carrier_count':0,'exact_positive':True,'exact_partition_unity':True,'exact_affine_precision':True,'exact_physical_C2_inherited_from_DY':True,'coefficients_checked':int(B.size),'min_coefficient_global':str(min(min(f) for f in basis)),'min_per_group':min_per_group,'column_scaled_sampled_condition_number':float(sv[0]/sv[-1]),'sample_rank':rank,'alpha':str(ALPHA),'theta':str(THETA),'moment_site_offset':str(L),'outer_anchor_sites':[[str(a) for a in base[i+1,:]] for i in range(n)],'refinement_uniform_stability':'not established','box_spline_transition':'not established'}
 print('PASS',name,'basis',len(basis),'cond',result['column_scaled_sampled_condition_number'],flush=True)
 return (result,basis) if return_basis else result
if __name__=='__main__':
 results={}
 for name,pts in [('fan5',points5),('fan9',points9)]:results[name]=run(pts,name)
 out=Path(__file__).parent/'v17_exact_results.json'
 out.write_text(json.dumps(results,indent=2)+'\n')
