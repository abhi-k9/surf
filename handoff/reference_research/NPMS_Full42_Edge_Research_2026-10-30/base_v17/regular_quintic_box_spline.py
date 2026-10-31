"""Independent regular-chart benchmark, not an arbitrary-fan construction.
Quintic C2 box spline with 3,2,2 directions on Z² triangular lattice.
Gauss-3 exactly integrates degree-four t-polynomials on each knot interval, floating-point representation.
"""
import math,json
from collections import Counter
from functools import lru_cache
from pathlib import Path
import numpy as np
from scipy.special import comb

def b_spline(n,x):
 if x<=0 or x>=n:return 0.0
 return sum((-1)**k*math.comb(n,k)*max(x-k,0.0)**(n-1) for k in range(n+1))/math.factorial(n-1)

@lru_cache(maxsize=200000)
def box(x,y):
 # M_{(3,2,2)} = ∫ b3(x-t) b2(y-t) b2(t) dt
 a=max(0.,x-3.,y-2.);b=min(2.,x,y)
 if a>=b:return 0.
 cuts=sorted({a,b,*[float(v) for v in [0.,1.,2.,x,x-1.,x-2.,x-3.,y,y-1.,y-2.] if a<v<b]})
 roots=[-math.sqrt(3./5.),0.,math.sqrt(3./5.)];weights=[5./9.,8./9.,5./9.]
 r=0.
 for lo,hi in zip(cuts,cuts[1:]):
  h=(hi-lo)/2;mid=(hi+lo)/2
  r+=h*sum(w*b_spline(3,x-(mid+h*z))*b_spline(2,y-(mid+h*z))*b_spline(2,mid+h*z) for z,w in zip(roots,weights))
 return r

def mask():
 c=Counter()
 for i in range(4):
  for j in range(3):
   for k in range(3):c[i+k,j+k]+=math.comb(3,i)*math.comb(2,j)*math.comb(2,k)
 return {xy:v/32 for xy,v in c.items()}

if __name__=='__main__':
 m=mask();rng=np.random.default_rng(20261030)
 errors={'partition':0,'affine_x':0,'affine_y':0,'refine':0,'mask_parity':0}
 for _ in range(50):
  x,y=rng.uniform(-1,6,2)
  p0=p1=p2=0.
  for i in range(math.floor(x)-5,math.ceil(x)+1):
   for j in range(math.floor(y)-4,math.ceil(y)+1):
    v=box(float(x-i),float(y-j));p0+=v;p1+=(i+2.5)*v;p2+=(j+2)*v
  errors['partition']=max(errors['partition'],abs(p0-1))
  errors['affine_x']=max(errors['affine_x'],abs(p1-x))
  errors['affine_y']=max(errors['affine_y'],abs(p2-y))
  v=box(float(x),float(y));f=sum(c*box(float(2*x-i),float(2*y-j)) for (i,j),c in m.items())
  errors['refine']=max(errors['refine'],abs(v-f))
 for parity1 in (0,1):
  for parity2 in (0,1):
   errors['mask_parity']=max(errors['mask_parity'],abs(sum(v for (i,j),v in m.items() if i%2==parity1 and j%2==parity2)-1))
 assert max(errors.values())<2e-12,errors
 assert len(m)>0 and sum(m.values())==4 and min(m.values())>=0
 result={'directions':'e1*3 + e2*2 + (e1+e2)*2','piecewise_total_degree':5,'C2':True,
         'regular_triangular_lattice_only':True,'refinement_mask_unique_entries':len(m),'mask_coefficients':[{ 'offset':list(k),'numerator_over_32':int(round(32*v))} for k,v in sorted(m.items())],
         'mask_sum':4,'parity_class_row_sums':1,'tested_points':50,'largest_numerical_errors':errors,
         'uniform_refinement_operator_Linf_norm':1,
         'basis_uniform_stability_note':'Uses published unimodular box-spline integer-translate theory; no analogous arbitrary-valence chart proof supplied.'}
 (Path(__file__).resolve().parent/'regular_box_spline_results.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k!='mask_coefficients'},indent=2))
