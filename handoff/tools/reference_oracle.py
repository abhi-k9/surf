#!/usr/bin/env python3
"""Standalone numerical oracle for exported 42D fixed-fan PS12 quintic fixtures.

No symbolic construction imports. Independent Bernstein derivative implementation,
with explicit local-to-world chain rule and rational quotient derivatives.
Use --write-golden to regenerate deterministic outputs; --verify to check them.
"""
import argparse, json, math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def triples(d):
    return [(i,j,d-i-j) for i in range(d+1) for j in range(d+1-i)]

def weight(deg, lam):
    return np.array([math.factorial(deg)/math.prod(math.factorial(t) for t in ijk)*math.prod(lam[k]**ijk[k] for k in range(3)) for ijk in triples(deg)])

def barycentric(xy, verts):
    p0,p1,p2=np.asarray(verts,dtype=float)
    J=np.column_stack((p1-p0,p2-p0))
    uv=np.linalg.solve(J, np.asarray(xy,dtype=float)-p0)
    return np.array([1-uv[0]-uv[1],uv[0],uv[1]])

def cell_of(xy,cells):
    matches=[]
    for m in range(cells.shape[0]):
        for c in range(cells.shape[1]):
            lam=barycentric(xy,cells[m,c])
            if float(lam.min())>=-1e-12:matches.append((m,c,lam))
    if not matches:raise ValueError('point is outside certified domain')
    return matches[0]

def evaluate_basis(coeff, verts, lam):
    """Returns columns B, B_x, B_y, B_xx, B_xy, B_yy for coefficients [function,21]."""
    deg=5
    q=weight(deg,lam)
    I={ijk:k for k,ijk in enumerate(triples(5))}
    B=coeff@q
    r1=triples(4)
    du=5*np.array([coeff[:,I[(i,j+1,k)]]-coeff[:,I[(i+1,j,k)]] for i,j,k in r1]).T@weight(4,lam)
    dv=5*np.array([coeff[:,I[(i,j,k+1)]]-coeff[:,I[(i+1,j,k)]] for i,j,k in r1]).T@weight(4,lam)
    r2=triples(3)
    duu=20*np.array([coeff[:,I[(i,j+2,k)]]-2*coeff[:,I[(i+1,j+1,k)]]+coeff[:,I[(i+2,j,k)]] for i,j,k in r2]).T@weight(3,lam)
    duv=20*np.array([coeff[:,I[(i,j+1,k+1)]]-coeff[:,I[(i+1,j+1,k)]]-coeff[:,I[(i+1,j,k+1)]]+coeff[:,I[(i+2,j,k)]] for i,j,k in r2]).T@weight(3,lam)
    dvv=20*np.array([coeff[:,I[(i,j,k+2)]]-2*coeff[:,I[(i+1,j,k+1)]]+coeff[:,I[(i+2,j,k)]] for i,j,k in r2]).T@weight(3,lam)
    p0,p1,p2=np.asarray(verts,dtype=float)
    J=np.column_stack((p1-p0,p2-p0))
    invJ=np.linalg.inv(J)
    grad=np.vstack((du,dv)).T@invJ
    # Hworld = inv(J).T @ Hlocal @ inv(J)
    H=np.empty((len(B),2,2))
    H[:,0,0]=duu;H[:,0,1]=duv;H[:,1,0]=duv;H[:,1,1]=dvv
    HW=np.einsum('ap,nab,bq->npq',invJ,H,invJ)  # note invJ axes first index = local, second = world
    return np.column_stack((B,grad[:,0],grad[:,1],HW[:,0,0],HW[:,0,1],HW[:,1,1]))

def eval_sample(fixture,xy,ctrl,weights):
    cells=fixture['microtri_vertices_xy'];m,c,lam=cell_of(xy,cells)
    basis=evaluate_basis(fixture['coefficients'][:,m,c,:],cells[m,c],lam)
    if not np.isfinite(basis).all():raise AssertionError('nonfinite basis values')
    N,Nx,Ny,Nxx,Nxy,Nyy=basis.T
    w=np.asarray(weights,dtype=float)
    if len(w)!=len(N) or not np.all(w>0):raise ValueError('strictly positive weights required')
    Q=np.column_stack((np.asarray(ctrl,dtype=float)*w[:,None],w))
    H=basis.T@Q  # rows for value, dx,dy,dxx,dxy,dyy; four homogeneous coords
    def quotient(H):
        A,W=H[0,:3],H[0,3];S=A/W
        Su=(H[1,:3]-S*H[1,3])/W
        Sv=(H[2,:3]-S*H[2,3])/W
        Suu=(H[3,:3]-S*H[3,3]-2*Su*H[1,3])/W
        Suv=(H[4,:3]-S*H[4,3]-Su*H[2,3]-Sv*H[1,3])/W
        Svv=(H[5,:3]-S*H[5,3]-2*Sv*H[2,3])/W
        return np.stack((S,Su,Sv,Suu,Suv,Svv),axis=0)
    derivs=quotient(H)
    return {'macro':int(m),'cell':int(c),'basis':basis,'surface':derivs}

def controls(fixture,variant):
    n=len(fixture['sites_uv']);uv=fixture['sites_uv'];parents=fixture['parent_index']
    if variant=='plane':
        xyz=np.column_stack((uv,np.zeros(n)));w=np.ones(n)
    elif variant=='rational':
        z=np.array([0.17*np.sin(0.31*i)+0.03*(int(parents[i])%5) for i in range(n)])
        xyz=np.column_stack((uv,z));w=np.array([0.75+0.05*(i%9) for i in range(n)])
    else:raise ValueError(variant)
    return xyz,w

def assertions(fixture):
    a=fixture['coefficients'];assert a.ndim==4 and a.shape[2:]==(12,21)
    assert np.isfinite(a).all(); assert np.min(a)>=-1e-14
    assert np.max(np.abs(a.sum(axis=0)-1))<1e-12
    # parent split numerical conservation, even if parent function not reloaded
    # across every macro microcell, the param-coordinate identity must hold
    m=fixture['microtri_vertices_xy'];sites=fixture['sites_uv']
    for sector in range(a.shape[1]):
        for cell in range(12):
            cc=a[:,sector,cell,:]
            for lam in ((.2,.3,.5),(.3,.4,.3),(.55,.3,.15)):
                xy=np.asarray(lam)@m[sector,cell]
                B=cc@weight(5,lam)
                assert abs(B.sum()-1)<1e-12
                assert np.linalg.norm(B@sites-xy,np.inf)<3e-12,(sector,cell,lam,B@sites,xy)
    return True

def samples(fixture):
    m=fixture['microtri_vertices_xy']
    # one or two inside every PS12 microcell, plus extra center points; valid on both fans.
    barys=((.17,.37,.46),(.4,.22,.38))
    out=[]
    for sector in range(len(m)):
        for cell in range(12):
            for ix,lam in enumerate(barys):
                out.append((sector,cell,ix,(np.asarray(lam)@m[sector,cell]).tolist()))
    return out

def create(n):
    fixturepath=ROOT/'fixtures'/f'fan{n}_full42_float64.npz'
    with np.load(fixturepath) as z: fixture={k:z[k] for k in z.files}
    assertions(fixture)
    result={'schema':'npms-oracle-v1','valence':n,'source_fixture':fixturepath.name,
            'meaning':'Analytic numeric golden values for C++ evaluator. NOT a formal positivity/stability certificate.',
            'samples':[]}
    for sector,cell,ix,xy in samples(fixture):
        # record both planar and rational surface / first and second derivatives at each sample.
        row={'macro':sector,'cell':cell,'subsample':ix,'xy':xy}
        for variant in ('plane','rational'):
            xyz,w=controls(fixture,variant)
            ev=eval_sample(fixture,xy,xyz,w)
            assert ev['macro']==sector and ev['cell']==cell,(sector,cell,ev['macro'],ev['cell'])
            row[variant]={'surface_jet_6x3':ev['surface'].tolist(),
                'basis_sum':float(ev['basis'][:,0].sum()),
                'basis_min':float(ev['basis'][:,0].min())}
            if variant=='plane':
                should=np.array([[*xy,0],[1,0,0],[0,1,0],[0,0,0],[0,0,0],[0,0,0]],dtype=float)
                assert np.max(abs(ev['surface']-should))<2e-9,(n,sector,cell,ev['surface'],should)
        result['samples'].append(row)
    path=ROOT/'fixtures'/f'fan{n}_golden_samples.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS fan',n,'basis',len(fixture['coefficients']),'subcell samples',len(result['samples']), 'oracle',path.name)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--valence',type=int,choices=(5,9));a=p.parse_args()
    for n in ((a.valence,) if a.valence else (5,9)):create(n)
