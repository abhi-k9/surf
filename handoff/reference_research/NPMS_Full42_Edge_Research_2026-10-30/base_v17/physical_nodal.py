"""Reference PS12 local nodal and physical-coordinate Bernstein extraction.
Independent exact rational construction; not a general Speleers positive basis.
"""
from pathlib import Path
import sympy as s
from math import factorial
from authentic_ps12 import V, cells, monos
Q=s.Rational
ROOT=Path(__file__).resolve().parent
L=s.Matrix([[Q(q) for q in line.split()] for line in (ROOT/'full_local_42.txt').read_text().splitlines()]).T
assert L.shape==(252,42)
jetkeys=[(a,k-a) for k in range(4) for a in range(k+1)]

def bary(pt, tri):
 (x,y)=pt; (ax,ay),(bx,by),(cx,cy)=tri
 den=(bx-ax)*(cy-ay)-(by-ay)*(cx-ax)
 q=((x-ax)*(cy-ay)-(y-ay)*(cx-ax))/den
 r=((bx-ax)*(y-ay)-(by-ay)*(x-ax))/den
 return 1-q-r,q,r

def cell_for(pt):
 return next(i for i,c in enumerate(cells) if min(bary(pt,[V[k] for k in c]))>=0)

def mon_jet(pt,dirs):
 x,y=pt;cs=[Q(1)]
 for dx,dy in dirs:
  cc=[Q(0)]*(len(cs)+1)
  for a,c in enumerate(cs):cc[a]+=c*dy;cc[a+1]+=c*dx
  cs=cc
 k=len(dirs);out=[]
 for i,j in monos:
  v=Q(0)
  for a,c in enumerate(cs):
   b=k-a
   if i>=a and j>=b:
    v+=c*(factorial(i)//factorial(i-a))*(factorial(j)//factorial(j-b))*x**(i-a)*y**(j-b)
  out.append(v)
 return s.Matrix([out])

def jet(pt,dirs):
 ci=cell_for(pt)
 return mon_jet(pt,dirs)*L[ci*21:(ci+1)*21,:]

# Global standard reference nodal degrees: 3x10 vertex jets and 3x4 edge jets
verts=[V[k] for k in ('v1','v2','v3')]
rows=[]
for pt in verts:
 for a,b in jetkeys:rows.append(jet(pt,[(1,0)]*a+[(0,1)]*b))
for i,j in ((0,1),(1,2),(2,0)):
 A=verts[i];B=verts[j];mid=tuple((A[k]+B[k])/2 for k in (0,1));t=(B[0]-A[0],B[1]-A[1]);n=(t[1],-t[0])
 for dirs in ([n],[n,n],[t,n,n],[n,n,n]):rows.append(jet(mid,dirs))
N=s.Matrix.vstack(*rows)
assert N.shape==(42,42) and N.rank()==42
Ninv=N.inv()
ijk=[(i,j,5-i-j) for i in range(6) for j in range(6-i)]
pts=[[Q(k,5) for k in triple] for triple in ijk]
BM=s.Matrix([[Q(factorial(5), factorial(r[0])*factorial(r[1])*factorial(r[2]))*s.prod(z**v for z,v in zip(p,r)) for r in ijk] for p in pts]);BMi=BM.inv()
BREF=[]
for ci,cell in enumerate(cells):
 tri=[V[z] for z in cell]
 r=[]
 for p in pts:
  x=sum(p[k]*tri[k][0] for k in range(3));y=sum(p[k]*tri[k][1] for k in range(3))
  r.append([x**i*y**j for i,j in monos])
 BREF.append(BMi*s.Matrix(r)*L[21*ci:21*(ci+1),:]*Ninv)
BREF=s.Matrix.vstack(*BREF)

def nodal_map(pi,pj):
 pi=s.Matrix(list(map(Q,pi)));pj=s.Matrix(list(map(Q,pj)));O=s.zeros(2,1)
 world=[O,pi,pj];local=[V['v1'],V['v2'],V['v3']]
 J=s.Matrix.hstack((pj-pi)/3,(-pi)/3);inv=J.inv();ex=s.Matrix([1,0]);ey=s.Matrix([0,1]);out=[]
 for p in local:
  for a,b in jetkeys:
   dirs=[inv*ex]*a+[inv*ey]*b
   out.append(jet(tuple(p),[tuple(v) for v in dirs])*Ninv)
 for i,j in ((0,1),(1,2),(2,0)):
  t=world[j]-world[i];n=s.Matrix([t[1],-t[0]]);tt=inv*t;nn=inv*n
  mid=tuple((local[i][k]+local[j][k])/2 for k in (0,1))
  for ds in ([nn],[nn,nn],[tt,nn,nn],[nn,nn,nn]):out.append(jet(mid,[tuple(v) for v in ds])*Ninv)
 Z=s.Matrix.vstack(*out)
 assert Z.shape==(42,42) and Z.rank()==42
 return Z

def physical_bernstein(pi,pj):
 return BREF*nodal_map(pi,pj).inv()

if __name__=='__main__':
 print('42D physical nodal reference rank:',N.rank(),'Bernstein extraction:',BREF.shape)
