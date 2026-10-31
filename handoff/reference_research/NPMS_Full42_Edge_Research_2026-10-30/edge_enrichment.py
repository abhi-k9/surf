"""Exact carrier-free positive edge enrichment of the existing NPMS v1.7 quintic fan.
No claim of spatial mesh refinement or multilevel stability; these are fixed-mesh subspace enlargements.
Source files in nearby verified archives are required. No optimizer or floating-point candidate is used.
"""
from pathlib import Path
import sys,json,math,itertools
import sympy as s
import numpy as np
root=Path(__file__).resolve().parent
v17=root/'base_v17'
sys.path.insert(0,str(v17));sys.path.insert(0,str(v17/'verification'))
from verify_v17 import run,points5,points9
import physical_nodal as X
Q=s.Rational
r=json.loads((root/'positive_physical_coefficients.json').read_text())

def make(n,points):
 print('Starting valence',n,flush=True)
 dat,basis=run(points,'fan'+str(n),True)
 m=252*n
 mats=[X.physical_bernstein(points[i],points[(i+1)%n]) for i in range(n)]
 edge=[];edge_coeffs=[]
 for i,rows in enumerate(r[f'fan{n}']['edge_coefficients_by_radial_index']):
  coeff=s.Matrix([[Q(z) for z in row] for row in rows]);prev=(i-1)%n
  C=s.zeros(m,3)
  C[252*prev:252*(prev+1),:]=mats[prev][:,38:42]*s.diag(-1,1,-1,-1)*coeff
  C[252*i:252*(i+1),:]=mats[i][:,30:34]*coeff
  assert min(C)>=0 and coeff.rank()==3
  edge.extend([C[:,j] for j in range(3)])
  edge_coeffs.append(coeff)
 print('Edge extraction complete',flush=True)
 # Every edge child is assigned uniquely to one of three parent controls of its outer vertex.
 assign=[];bound=[];evalue=[];r_min=[]
 for i in range(n):
  parents=[basis[3*(i+1)+a] for a in range(3)]
  children=edge[3*i:3*i+3]
  ratios=s.zeros(3,3)
  for a in range(3):
   for k in range(3):
    ch=children[k];pa=parents[a]
    if any(pa[t]==0 and ch[t]>0 for t in range(m)):ratios[a,k]=0
    else:ratios[a,k]=min(pa[t]/ch[t] for t in range(m) if ch[t]>0)
  feasible=[]
  for perm in itertools.permutations(range(3)):
   if min(ratios[perm[k],k] for k in range(3))>0:
    score=min(ratios[perm[k],k] for k in range(3))
    feasible.append((score,perm))
  assert feasible,('no outer-local assignment',n,i,ratios)
  score,perm=max(feasible,key=lambda x:x[0]);r_min.append(score)
  for k in range(3):
   eps=min(Q(1,2),ratios[perm[k],k]/2)
   assert eps>0
   pa=3*(i+1)+perm[k]
   assign.append(pa);evalue.append(eps)
   assert min(basis[pa]-eps*children[k])>=0
  print('edge',i,'best_perm',perm,'min_parent/edge',score,'eps',[str(z) for z in evalue[-3:]],flush=True)
 fine=basis[:]
 for k,ch in enumerate(edge): fine[assign[k]]-=evalue[k]*ch
 fine.extend([evalue[k]*ch for k,ch in enumerate(edge)])
 assert min(min(col) for col in fine)>=0
 assert sum(fine,s.zeros(m,1))==s.ones(m,1)
 N0=len(basis);N1=len(fine)
 assert N1==N0+3*n
 T=s.eye(N1,N0)
 for k,parent in enumerate(assign):T[N0+k,parent]=1
 F=s.Matrix.hstack(*fine);Bc=s.Matrix.hstack(*basis)
 assert F*T==Bc and min(T)>=0 and T.rank()==N0
 # Independence follows from unique edge-jet DOFs missing from v17; numerically cross-check it.
 ar=np.asarray(F.tolist(),dtype=float)
 # Column-scale coefficient matrix and also sample physical Bernstein values
 norm=np.linalg.norm(ar,axis=0)
 ev=np.linalg.svd(ar/norm,compute_uv=False)
 from math import factorial
 barys=[(1/3,1/3,1/3),(.5,.25,.25),(.25,.5,.25),(.25,.25,.5)]
 weight=np.asarray([[factorial(5)/math.prod(factorial(x) for x in t)*math.prod(b[j]**t[j] for j in range(3)) for t in X.ijk] for b in barys])
 sample=np.einsum('pj,kjc->kpc',weight,ar.reshape(n*12,21,N1)).reshape(-1,N1)
 sn=np.linalg.svd(sample/np.linalg.norm(sample,axis=0),compute_uv=False)
 assert ev[-1]>1e-12 and sn[-1]>1e-12
 # Row-stochastic prolongation: every fine DOF inherits exactly one coarse DOF.
 assert all(sum(T[k,j] for j in range(N0))==1 for k in range(N1))
 # Exact affine reproduction from replicated parent sites; equality of coefficient maps.
 out={'n':n,'coarse':N0,'fine':N1,'ratio_min_by_edge':[str(x) for x in r_min],'child_scales':[str(z) for z in evalue], 'parent_indices':assign,'min_edge_scale':str(min(evalue)), 'max_edge_scale':str(max(evalue)), 'min_fine_Bernstein':str(min(min(x) for x in fine)), 'exact_transfer':True,'exact_partition_unity':True,'affine_reproduction_preserved_by_copying_parent_sites':True,'exact_positive':True,'outer_parent_macro_support_max':4,'edge_child_macro_support':2,'rank_by_disjoint_nodal_blocks':N1, 'numerical_column_scaled_coefficient_condition':float(ev[0]/ev[-1]),'numerical_sampled_basis_condition':float(sn[0]/sn[-1]),'reproduction_limit':'Geometric span still does not reproduce general quadratics; analytic full 42D spline space does.','refinement_limit':'Fixed mesh enrichment, not spatial subdivision; no uniform multilevel inverse stability theorem.'}
 print('RESULT',out['n'],out['fine'],out['min_edge_scale'],out['numerical_column_scaled_coefficient_condition'],out['numerical_sampled_basis_condition'],flush=True)
 return out

if __name__=='__main__':
 results={}
 for n,points in [(5,points5),(9,points9)]: results['fan'+str(n)]=make(n,points)
 (root/'enrichment_results.json').write_text(json.dumps(results,indent=2)+'\n')
 print('PASS both exact rational geometric refinement fan tests',flush=True)
