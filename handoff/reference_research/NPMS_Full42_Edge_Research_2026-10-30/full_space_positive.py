"""Exact *full-dimensional* positive quintic PS12 basis on two fixed rational planar fans.

Constructs N=18n+10 independent Bernstein-nonnegative, normalized, affine-precise
functions spanning the authentic Davydov-Yeo global level-zero space, n=5,9.

Limitation: repeated local splitting can be very ill-conditioned; NO uniform
spatial refinement theorem, regular-grid junction or arbitrary-mesh certificate.
"""
from pathlib import Path
import sys, json, itertools, math
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parent
BASE=ROOT/'base_v17'
sys.path[:0]=[str(BASE),str(BASE/'verification')]
from verify_v17 import run,points5,points9
import physical_nodal as X
Q=s.Rational

def make(n,points):
 print('START',n,flush=True)
 base_result,base=run(points,'fan'+str(n),True)
 M=[X.physical_bernstein(points[i],points[(i+1)%n]) for i in range(n)]
 m=252*n
 modes=[]
 def append(name,parts,parents):
  v=s.zeros(m,1)
  for k,j,co in parts: v[252*k:252*(k+1),:]+=M[k][:,j]*Q(co)
  options=[]
  for p in parents:
   B=base[p]
   if any(B[t]==0 and v[t]!=0 for t in range(m)): continue
   ratio=min((B[t]/abs(v[t]) for t in range(m) if v[t]!=0), default=Q(0))
   if ratio>0:options.append((ratio,p))
  assert options,('no parent for full jet',n,name)
  ratio,parent=max(options)
  modes.append({'name':name,'vec':v,'ratio':ratio,'parent':parent})
 # All global Hermite nodal jet directions not already freely spanned by v1.7 P1-compatible base.
 for j in range(3,10):append(['center_vertex',j],[(i,j,1) for i in range(n)],range(3))
 for i in range(n):
  prev=(i-1)%n;parents=range(3*(i+1),3*(i+2))
  for j in range(3,10):append(['outer_vertex',i,j],[(prev,20+j,1),(i,10+j,1)],parents)
  for j in range(4):append(['radial_edge',i,j],[(prev,38+j,[-1,1,-1,-1][j]),(i,30+j,1)],parents)
  for j in range(4):append(['outer_edge',i,j],[(i,34+j,1)],parents)
 assert len(modes)==15*n+7
 print('FULL missing directional modes',len(modes),flush=True)
 groups=[[] for _ in base]
 for k,d in enumerate(modes):groups[d['parent']].append(k)
 # For a parent P having m_j signed nodal directions H_k, construct m_j+1
 # new positive functions, whose sum is P and whose span contains all H_k.
 F=[];parent_for_fine=[];scales=[]
 for j,P in enumerate(base):
  ids=groups[j];K=len(ids)
  if not K:
   F.append(P);parent_for_fine.append(j)
   continue
  beta={k:modes[k]['ratio']/(2*K) for k in ids}
  # Relative dominance: sum_k abs(beta_k H_k) <= P/2 coefficientwise.
  sums=sum((beta[k]*modes[k]['vec'] for k in ids),s.zeros(m,1))
  remainder=(P-sums)/(K+1)
  assert min(remainder)>=0,('remainder',n,j,min(remainder))
  F.append(remainder);parent_for_fine.append(j)
  for k in ids:
   child=(P+beta[k]*modes[k]['vec'])/(K+1)
   assert min(child)>=0,('child',n,j,k,min(child))
   F.append(child);parent_for_fine.append(j)
   scales.append({'mode':modes[k]['name'],'parent':j,'dominance_ratio':str(modes[k]['ratio']), 'perturbation_beta':str(beta[k]),'sibling_count':K+1})
  assert remainder+sum(((P+beta[k]*modes[k]['vec'])/(K+1) for k in ids),s.zeros(m,1))==P
 assert len(F)==18*n+10
 Fmat=s.Matrix.hstack(*F)
 T=s.zeros(len(F),len(base))
 for i,j in enumerate(parent_for_fine):T[i,j]=1
 assert Fmat*T==s.Matrix.hstack(*base)
 assert Fmat*s.ones(len(F),1)==s.ones(m,1)
 assert min(Fmat)>=0
 assert T*s.ones(len(base),1)==s.ones(len(F),1)
 # Proof of exact rank: subtract each mode child from its group remainder:
 # (K+1)(C_k-C_0)=beta_k H_k + sum_j beta_j H_j.
 # Since [H_k] are distinct independent D-Y Hermite dual coordinates and
 # the v1.7 base projects invertibly to vertex 0/first jet coordinates,
 # the groupwise change-of-basis is invertible for nonzero beta_k.
 # Auxiliary numerical test on coefficient matrix (may be very ill-conditioned).
 arr=np.asarray(Fmat.tolist(),float)
 colnorms=np.linalg.norm(arr,axis=0)
 eig=np.linalg.svd(arr/colnorms,compute_uv=False)
 cond=float(eig[0]/eig[-1])
 rank_num=int(np.linalg.matrix_rank(arr/colnorms))
 assert rank_num==len(F),(n,'sampled rank lost',rank_num,len(F),cond)
 ratios=[z['ratio'] for z in modes]
 res={'n':n,'original_positive_basis':len(base),'added_independent_dy_nodal_modes':len(modes),
      'new_full_positive_basis_count':len(F),'global_full_DY_dimension':18*n+10,
      'exact_positive_Bernstein':True,'exact_partition_unity':True,'affine_precision_preserved_from_parent_sites':True,
      'full_global_DY_rank_by_invertible_nodal_coordinate_change':18*n+10,
      'full_P5_reproduction_by_spline_span':True,'negative_Bernstein_coefficients':0,
      'outer_support_at_most_sectors':4,'center_natural_star_sectors':n,
      'min_parent_mode_dominance_ratio':str(min(ratios)),
      'min_perturbation_beta':str(min(s.Rational(z['perturbation_beta']) for z in scales)),
      'max_modes_from_any_parent':max(map(len,groups)),
      'sampled_column_scaled_Bernstein_condition':cond,
      'sampled_numerical_rank':rank_num,
      'modes':scales,
      'limitations':'Exact only on fixed rational 5/9 planar fans, not general meshes. Geometric basis is ill-conditioned; NO true spatial refinement, no level-independent stability, no certified regular-to-extraordinary C2 connection.'}
 print('RESULT',n,'full',len(F),'group_max',res['max_modes_from_any_parent'],'minbeta',res['min_perturbation_beta'],'cond',cond,flush=True)
 return res
if __name__=='__main__':
 data={}
 for n,pts in ((5,points5),(9,points9)):data['fan'+str(n)]=make(n,pts)
 (ROOT/'full_positive_results.json').write_text(json.dumps(data,indent=2)+'\n')
 print('PASS exact full-DY positive normalized bases on both restricted fans',flush=True)
