"""Exact finite-field dimensions of the Davydov-Yeo S0 quintic PS12 space
on arbitrary-valence planar fans (5, 9). Published space, NOT Speleers positive basis.
No sampled/surface-plot claim is used as proof of dimension.
"""
import numpy as np
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
from authentic_ps12 import rows, rank_mod, V, cells, monos

@lru_cache(None)
def nullspace_mod(p):
 R,_=rows();n=252
 M=np.zeros((len(R),n),dtype=np.int64)
 for k,row in enumerate(R):
  for j,q in row.items():M[k,j]=int(q.numerator%p)*pow(int(q.denominator%p),-1,p)%p
 rr=0;pivot=[]
 for col in range(n):
  ids=np.where(M[rr:,col]!=0)[0]
  if not len(ids):continue
  k=rr+int(ids[0]);M[[rr,k]]=M[[k,rr]]
  M[rr,col:]=M[rr,col:]*pow(int(M[rr,col]),-1,p)%p
  below=np.flatnonzero(M[rr+1:,col])+rr+1
  if len(below):M[below,col:]=(M[below,col:]-M[below,col,None]*M[rr,None,col:])%p
  pivot.append(col);rr+=1
  if rr==min(M.shape):break
 # above pivot backwards for RREF
 for ix in range(rr-1,-1,-1):
  col=pivot[ix];ind=np.flatnonzero(M[:ix,col])
  if len(ind):M[ind,col:]=(M[ind,col:]-M[ind,col,None]*M[ix,None,col:])%p
 free=[j for j in range(n) if j not in set(pivot)]
 N=np.zeros((n,len(free)),dtype=np.int64)
 for k,f in enumerate(free):
  N[f,k]=1
  N[pivot,k]=-M[np.arange(rr),f]%p
 assert N.shape==(252,42) and np.count_nonzero(M@N%p)==0
 return N

points5=[(2,0),(1,2),(-1,2),(-2,-1),(0,-2)]
points9=[(3,0),(3,2),(1,3),(-1,3),(-3,2),(-3,0),(-2,-3),(0,-3),(2,-2)]

def poly_jet_at(x,y,dirs,p):
 # polynomial p in reference (x,y), directional derivative sequence 'dirs' in reference
 k=len(dirs);dcoefs=[1]
 for vx,vy in dirs:
  new=[0]*(len(dcoefs)+1)
  for a,c in enumerate(dcoefs):
   new[a]=(new[a]+c*vy)%p
   new[a+1]=(new[a+1]+c*vx)%p
  dcoefs=new
 out=np.zeros(21,dtype=np.int64)
 for m,(i,j) in enumerate(monos):
  val=0
  for a,c in enumerate(dcoefs):
   b=k-a
   if i>=a and j>=b:
    v= factorial(i)//factorial(i-a)*factorial(j)//factorial(j-b)
    val+=(c*v*pow(x,i-a,p)*pow(y,j-b,p))
  out[m]=val%p
 return out

def affine_inverse(A,B,C,p):
 # ref vertices: (0,3)-> A, (0,0)-> B, (3,0)-> C
 ax=(C[0]-B[0]);ay=(C[1]-B[1]);bx=(A[0]-B[0]);by=(A[1]-B[1])
 d=ax*by-ay*bx
 assert d%p
 inv=pow(d%p,-1,p)
 # coordinates ref = 3/d * [by -bx; -ay ax] (world-B)
 return np.array([[3*by*inv,-3*bx*inv],[-3*ay*inv,3*ax*inv]],dtype=np.int64)%p

def to_ref(xy,B,J,p):
 v=np.array([int(F(k-b).numerator%p)*pow(int(F(k-b).denominator%p),-1,p)%p for k,b in zip(xy,B)],dtype=np.int64)
 return (J@v)%p

def basis_functional(N,A,B,C,xy,global_dirs,p):
 J=affine_inverse(A,B,C,p);u,v=map(int,to_ref(xy,B,J,p))
 dd=[tuple(map(int,(J@np.array(d,dtype=np.int64))%p)) for d in global_dirs]
 # select local piece whose closure contains the ref point via reference barycentric triangle test
 # choose cells containing feature: because evaluated point is one of main vertex/midpoint, contained in cells.
 chosen=None
 for k,tri in enumerate(cells):
  xyz=[V[z] for z in tri]
  # check barycentric in rationals from u,v cast to p; locate via global coords using exact preimage in ref
  # here test integer feature x,y matching names in reference modulo prime
  if any((int(V[z][0].numerator*pow(V[z][0].denominator,-1,p))%p==u and int(V[z][1].numerator*pow(V[z][1].denominator,-1,p))%p==v) for z in tri):
   chosen=k;break
 if chosen is None:
  raise ValueError(('feature not PS node',u,v))
 mon=poly_jet_at(u,v,dd,p)
 return (mon @ N[chosen*21:(chosen+1)*21,:])%p

def modrank(A,p):
 A=A.copy()%p;r=0
 for c in range(A.shape[1]):
  ix=np.where(A[r:,c]!=0)[0]
  if not len(ix):continue
  j=r+int(ix[0]);A[[r,j]]=A[[j,r]]
  A[r,c:]=A[r,c:]*pow(int(A[r,c]),-1,p)%p
  idx=np.where(A[r+1:,c]!=0)[0]+r+1
  if len(idx):A[idx,c:]=(A[idx,c:]-A[idx,c,None]*A[r,None,c:])%p
  r+=1
  if r>=A.shape[0]:break
 return r

def fan_rows(points,p):
 n=len(points);N=nullspace_mod(p);O=(0,0)
 # each macro T_i=(center,P_i,P_i+1)
 macro=[(O,points[i],points[(i+1)%n]) for i in range(n)]
 globalrows=[]
 def append_functional(a,b,xy,dirs):
  ra=basis_functional(N,*macro[a],xy,dirs,p)
  rb=basis_functional(N,*macro[b],xy,dirs,p)
  r=np.zeros(42*n,dtype=np.int64)
  r[a*42:(a+1)*42]=ra
  r[b*42:(b+1)*42]=-rb%p
  globalrows.append(r)
 # each radial edge has two adjacent macrotriangles
 for i in range(n):
  a=(i-1)%n;b=i;P=points[i];mid=(F(P[0],2),F(P[1],2))
  # C3 Taylor jet continuity at both endpoints (10 each)
  for xy in (O,P):
   for k in range(4):
    for dx in range(k+1):
     dy=k-dx
     append_functional(a,b,xy,[(1,0)]*dx+[(0,1)]*dy)
  # 4 C2-edge nodal MDS degrees of freedom at the midpoint
  tangent=P;normal=(P[1],-P[0]);
  for dirs in ([normal],[normal,normal],[tangent,normal,normal],[normal,normal,normal]):
   append_functional(a,b,mid,dirs)
 return np.vstack(globalrows)

if __name__=='__main__':
 for n,pts in [(5,points5),(9,points9)]:
  for p in (1009,2017):
   A=fan_rows(pts,p)
   r=modrank(A,p); expected=42*n-(10*(n+1)+4*2*n)
   print('fan',n,'prime',p,'matrix',A.shape,'rank',r,'dimension',A.shape[1]-r,'published dimension',10*(n+1)+8*n,'expected row rank',expected,flush=True)


def local_nodal_matrix(p):
 N=nullspace_mod(p);A=(0,3);B=(0,0);C=(3,0);M=[]
 for name in ('v1','v2','v3'):
  xy=V[name]
  for k in range(4):
   for a in range(k+1):
    M.append(basis_functional(N,A,B,C,xy,[(1,0)]*a+[(0,1)]*(k-a),p))
 for va,vb in [('v1','v2'),('v2','v3'),('v3','v1')]:
  pa,pb=V[va],V[vb];mid=tuple((a+b)/2 for a,b in zip(pa,pb));tan=tuple(b-a for a,b in zip(pa,pb));nor=(tan[1],-tan[0])
  for ds in [[nor],[nor,nor],[tan,nor,nor],[nor,nor,nor]]:
   M.append(basis_functional(N,A,B,C,mid,ds,p))
 return np.vstack(M)
