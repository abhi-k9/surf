"""Independent exact-rational assembly of the published Davydov-Yeo PS12 C2/C3 local space.
It reproduces the degree-five local *space*, NOT the Speleers normalized positive basis.
"""
from fractions import Fraction as F
from collections import defaultdict
from math import comb, factorial
import numpy as np

V={'v1':(F(0),F(3)), 'v2':(F(0),F(0)), 'v3':(F(3),F(0))}
mid=lambda a,b:tuple((u+v)/2 for u,v in zip(a,b))
V['w1']=mid(V['v2'],V['v3']);V['w2']=mid(V['v1'],V['v3']);V['w3']=mid(V['v1'],V['v2'])
V['u1']=mid(V['w3'],V['w2']);V['u2']=mid(V['w1'],V['w3']);V['u3']=mid(V['w1'],V['w2']);V['z']=(F(1),F(1))
# 12 PS12 sub-triangles, avoiding any overlay approximations
cells=[('v1','w3','u1'),('v1','u1','w2'),('v2','w1','u2'),('v2','u2','w3'),('v3','w2','u3'),('v3','u3','w1'),
       ('z','w3','u1'),('z','u1','w2'),('z','w1','u2'),('z','u2','w3'),('z','w2','u3'),('z','u3','w1')]
edge_cells=defaultdict(list)
for k,cell in enumerate(cells):
 for a,b in [(cell[0],cell[1]),(cell[1],cell[2]),(cell[2],cell[0])]:edge_cells[tuple(sorted((a,b)))].append(k)
type3={tuple(sorted(x)) for x in [('v1','u1'),('v2','u2'),('v3','u3'),('u1','z'),('u2','z'),('u3','z'),('w1','z'),('w2','z'),('w3','z')]}
monos=[(i,j) for i in range(6) for j in range(6-i)] # 21
assert len(monos)==21 and len(cells)==12 and len(V)==10

def deriv_row(x,y,nx,ny,k):
 ans=[]
 for i,j in monos:
  v=F(0)
  for a in range(k+1):
   b=k-a
   if i>=a and j>=b:
    v+=F(comb(k,a)*factorial(i)//factorial(i-a)*factorial(j)//factorial(j-b))*nx**a*ny**b*x**(i-a)*y**(j-b)
  ans.append(v)
 return ans

def rows():
 out=[]; descriptions=[]
 for (a,b),adj in edge_cells.items():
  if len(adj)!=2:continue
  ca,cb=adj;xa,ya=V[a];xb,yb=V[b]
  dx=xb-xa;dy=yb-ya;nx=dy;ny=-dx
  for order in range(4 if (a,b) in type3 else 3):
   # degree of restricted k-th derivative is at most 5-k
   for kk in range(6-order):
    t=F(kk+1,7-order);x=xa+t*dx;y=ya+t*dy
    r=deriv_row(x,y,nx,ny,order)
    row={ca*21+i:q for i,q in enumerate(r) if q}
    row.update({cb*21+i:-q for i,q in enumerate(r) if q})
    out.append(row);descriptions.append(('edge',a,b,order,str(t)))
 # Additional C3 vertex jets at original triangle vertices and edge midpoints.
 # C3 across designated type-1/2/3 edges alone leaves a 6-DOF excess.
 for name in ['v1','v2','v3','w1','w2','w3']:
  touching=[ci for ci,c in enumerate(cells) if name in c]
  ref=touching[0];x,y=V[name]
  for other in touching[1:]:
   for aa in range(4):
    bb=3-aa
    jr=[]
    for i,j in monos:
     jr.append(F(factorial(i)//factorial(i-aa)*factorial(j)//factorial(j-bb))*x**(i-aa)*y**(j-bb) if i>=aa and j>=bb else F(0))
    row={ref*21+i:q for i,q in enumerate(jr) if q}
    row.update({other*21+i:-q for i,q in enumerate(jr) if q})
    out.append(row);descriptions.append(('vertex_C3',name,other,aa,bb))
 return out,descriptions

def rank_mod(rows,p):
 n=12*21
 M=np.zeros((len(rows),n),dtype=np.int64)
 for k,r in enumerate(rows):
  for j,q in r.items(): M[k,j]=int(q.numerator%p)*pow(int(q.denominator%p),-1,p)%p
 r=0;piv=[]
 for c in range(n):
  ids=np.where(M[r:,c]!=0)[0]
  if len(ids)==0:continue
  i=r+int(ids[0]);M[[r,i]]=M[[i,r]]
  M[r,c:]=(M[r,c:]*pow(int(M[r,c]),-1,p))%p
  below=np.flatnonzero(M[r+1:,c])+r+1
  if len(below):M[below,c:]=(M[below,c:]-M[below,c,None]*M[r,None,c:])%p
  piv.append(c);r+=1
  if r==min(M.shape):break
 return r, piv
if __name__=='__main__':
 R,desc=rows()
 print('vertices',len(V),'cells',len(cells),'internal edges',sum(len(q)==2 for q in edge_cells.values()),'constraints',len(R),'DOF',len(cells)*21)
 print('types',sum(len(q)==2 and k in type3 for k,q in edge_cells.items()))
 for p in (1009,2017):
  rank,piv=rank_mod(R,p)
  print('rank_mod',p,rank,'dimension upper bound',252-rank)
