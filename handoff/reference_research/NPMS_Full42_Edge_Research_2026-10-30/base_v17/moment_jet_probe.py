import sys,time
sys.path.insert(0,str(__import__('pathlib').Path(__file__).resolve().parent/'verification'))
import numpy as np
import sympy as s
from fan_assembly import points5,points9
from authentic_ps12 import V,cells
import physical_nodal as P
Q=s.Rational

def build(points):
 n=len(points)
 mats=[]
 for i in range(n):
  print('matrix',n,i,flush=True)
  mats.append(P.physical_bernstein(points[i],points[(i+1)%n]))
 def global_col(parts):
  out=[s.zeros(252,1) for _ in range(n)]
  for i,col,c in parts:out[i]+=mats[i][:,col]*c
  return s.Matrix.vstack(*out)
 VV=[]; DX=[];DY=[]
 VV.append(global_col([(i,0,Q(1)) for i in range(n)]))
 DX.append(global_col([(i,2,Q(1)) for i in range(n)]));DY.append(global_col([(i,1,Q(1)) for i in range(n)]))
 for i in range(n):
  p=(i-1)%n
  VV.append(global_col([(p,20,Q(1)),(i,10,Q(1))]))
  DX.append(global_col([(p,22,Q(1)),(i,12,Q(1))]))
  DY.append(global_col([(p,21,Q(1)),(i,11,Q(1))]))
 EX=[]; normals=[]
 for i in range(n):
  p=(i-1)%n
  # global normal on center to outer P_i
  normal=(Q(points[i][1]),-Q(points[i][0]));normals.append(normal)
  EX.append(global_col([(p,38,Q(-1)),(i,30,Q(1))]))
 for i in range(n):
  # global edge from P_i to P_i+1
  a,b=points[i],points[(i+1)%n]
  normal=(Q(b[1]-a[1]),Q(a[0]-b[0]));normals.append(normal)
  EX.append(global_col([(i,34,Q(1))]))
 # check affine decomposition
 unit=sum(VV,s.zeros(252*n,1));assert unit==s.ones(252*n,1)
 A=[]
 for i in range(n):
  a,b=points[i],points[(i+1)%n]
  for cell in cells:
   verts=[V[z] for z in cell]
   for lam in P.pts:
    z=sum((lam[j]*s.Matrix(verts[j]) for j in range(3)),s.zeros(2,1))
    x=s.Matrix(a)+(s.Matrix(b)-s.Matrix(a))*z[0]/3-s.Matrix(a)*z[1]/3
    A.append((x[0],x[1]))
 for axis in [0,1]:
  approx=sum((Q(points[i][axis])*VV[i+1] for i in range(n)),s.zeros(252*n,1))
  approx+=sum((DX if axis==0 else DY),s.zeros(252*n,1))
  approx+=sum((normals[e][axis]*EX[e] for e in range(2*n)),s.zeros(252*n,1))
  want=s.Matrix([v[axis] for v in A]);assert approx==want,('affine failed',axis)
 print('verified exact complete affine Hermite decomposition',n,flush=True)
 return VV,DX,DY,EX,normals

if __name__=='__main__':
 for name,p in [('five',points5),('nine',points9)]:
  VV,DX,DY,EX,normals=build(p)
  n=len(p)
  for i in range(n+1):
   # moment group i=0 center; else outer vertex i-1 with its own radial edge and outer edge
   if i==0:
    Hx=DX[0];Hy=DY[0];parent=VV[0]
   else:
    k=i-1
    Hx=DX[i]+normals[k][0]*EX[k]+normals[n+k][0]*EX[n+k]
    Hy=DY[i]+normals[k][1]*EX[k]+normals[n+k][1]*EX[n+k]
    parent=VV[i]
   arr=np.array(parent,float).ravel()
   hx=np.array(Hx,float).ravel();hy=np.array(Hy,float).ravel()
   mask=abs(arr)<1e-10
   print(name,i,'pure parent zeros',int(mask.sum()),'max H on pure zero',round(max(np.maximum(abs(hx[mask]),abs(hy[mask]))),7),flush=True)
   if i>0:
    k=i-1
    neighbors=[VV[1+(k-1)%n],VV[1+(k+1)%n]]
    W=parent*Q(7,8)+sum((z*Q(1,16) for z in neighbors),s.zeros(252*n,1))
    w=np.array(W,float).ravel()
    m=abs(w)<1e-10
    print(' mixed zeros',int(m.sum()),'maxHmixedzero',round(max(np.maximum(abs(hx[m]),abs(hy[m]))),7),'maxratio',max(np.maximum(abs(hx[w>1e-10]),abs(hy[w>1e-10]))/w[w>1e-10]),flush=True)
