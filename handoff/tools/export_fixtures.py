#!/usr/bin/env python3
"""Export authored full 42D NPMS fixed-fan demonstration basis as portable float64 fixtures.

This is an engineering data exporter, NOT an independent certificate of spline
positivity or a universal NPMS constructor. The exact source comes from the
bundled research snapshots and their rational proofs. Do not silently change
Bernstein ordering or reinterpret parameters.

Run with source archive extracted in reference_research/NPMS_Full42_Edge_Research_2026-10-30
  OPENBLAS_NUM_THREADS=1 python tools/export_fixtures.py
"""
import argparse
import json
import sys
from pathlib import Path
import hashlib
import numpy as np
import sympy as s

OUT = Path(__file__).resolve().parents[1]
SRC = OUT / 'reference_research' / 'NPMS_Full42_Edge_Research_2026-10-30'
sys.path.insert(0, str(SRC / 'base_v17'))
sys.path.insert(0, str(SRC))
from verify_v17 import run, points5, points9
import physical_nodal as X
from authentic_ps12 import V, cells

Q = s.Rational


def generate(n, points):
    res, base = run(points, f'fan{n}', True)
    reference = json.loads((SRC / 'full_positive_results.json').read_text())[f'fan{n}']
    mats = [X.physical_bernstein(points[i], points[(i + 1) % n]) for i in range(n)]
    numrows = 252*n
    grouped = [[] for _ in base]
    for desc in reference['modes']:
        kind, *idx = desc['mode']
        parts = []
        if kind == 'center_vertex':
            j = idx[0]
            parts = [(i, j, 1) for i in range(n)]
        elif kind == 'outer_vertex':
            i, j = idx
            parts = [((i-1)%n,20+j,1),(i,10+j,1)]
        elif kind == 'radial_edge':
            i, j = idx
            parts = [((i-1)%n,38+j,[-1,1,-1,-1][j]),(i,30+j,1)]
        elif kind == 'outer_edge':
            i, j = idx
            parts = [(i,34+j,1)]
        else:
            raise ValueError(desc['mode'])
        H = s.zeros(numrows, 1)
        for sector, col, sign in parts:
            H[sector*252:(sector+1)*252, 0] += sign*mats[sector][:,col]
        b = Q(desc['perturbation_beta'])
        grouped[desc['parent']].append((desc, b, H))
    parents = []
    basis = []
    grouped_sites = []
    # reconstruct exact site params from independently derived affine moments
    # The original v1.7 basis is three children of each vertex's W_i.
    A = s.zeros(n,n)
    ALPHA=Q(1,32); THETA=Q(1,5); L=Q(8)
    for i in range(n):
        A[i,i]=1-2*THETA
        A[i,(i-1)%n]=THETA
        A[i,(i+1)%n]=THETA
    P = s.Matrix([[Q(c) for c in p] for p in points]); q = A.T.inv()*P/(1-ALPHA)
    # three center sites centered at origin; other site anchors q_i
    for i in range(n+1):
        x,y=(Q(0),Q(0)) if i==0 else (q[i-1,0], q[i-1,1])
        grouped_sites.extend([(x+L,y),(x,y+L),(x-L,y-L)])
    assert len(grouped_sites)==len(base)
    control_groups=[]
    for j, parent in enumerate(base):
        group = grouped[j]
        k = len(group)
        if k:
            remainder=(parent - sum((beta*H for _,beta,H in group), s.zeros(numrows,1)))/(k+1)
            basis.append(remainder);parents.append(j);control_groups.append({'parent_index':j,'mode':'remainder','group_size':k+1})
            for desc,beta,H in group:
                basis.append((parent+beta*H)/(k+1)); parents.append(j)
                control_groups.append({'parent_index':j,'mode':desc['mode'],'group_size':k+1})
        else:
            basis.append(parent);parents.append(j)
            control_groups.append({'parent_index':j,'mode':'original','group_size':1})
    assert len(basis)==18*n+10
    assert len(parents)==len(basis)
    F=s.Matrix.hstack(*basis)
    assert min(F)>=0
    assert F*s.ones(len(basis),1)==s.ones(numrows,1)
    assert all(sum((basis[i] for i in range(len(basis)) if parents[i]==j),s.zeros(numrows,1))==base[j] for j in range(len(base)))
    sites=[grouped_sites[j] for j in parents]
    # exact affine coefficients: Bernstein representation of x and y from v1.7 sites
    for axis in range(2):
        base_affine=sum((grouped_sites[j][axis]*base[j] for j in range(len(base))),s.zeros(numrows,1))
        full_affine=F*s.Matrix([s[axis] for s in sites])
        assert full_affine==base_affine
    values=np.array(F.tolist(),dtype=np.float64)
    assert values.shape==(252*n,18*n+10)
    assert np.isfinite(values).all() and float(np.min(values)) >= 0
    coeff=values.T.reshape(len(basis),n,12,21) # axis convention: function,macro,cell,bernstein_ijk
    triples=[list(q) for q in X.ijk]
    macros=[]
    for i in range(n):
        a,b=points[i],points[(i+1)%n]
        macros.append([[0.,0.],list(map(float,a)),list(map(float,b))])
    # physical map from reference (local v1=(0,3),v2=(0,0),v3=(3,0)) to world triangle=(O,Pi,Pj)
    reference_vertices={k:[float(v) for v in pt] for k,pt in V.items()}
    microtriangles=[]
    for macro in macros:
        o,p,qv=np.array(macro,dtype=float)
        coords={key:(p + (qv-p)*np.array(ref)[0]/3-p*np.array(ref)[1]/3).tolist() for key,ref in reference_vertices.items()}
        microtriangles.append([[coords[v] for v in cell] for cell in cells])
    rounded_sites=np.array([[float(t[0]),float(t[1])] for t in sites],dtype=float)
    # return everything agent needs to support UI and fallback numeric evaluation
    npzpath=OUT/'fixtures'/f'fan{n}_full42_float64.npz'
    np.savez_compressed(npzpath,coefficients=coeff, sites_uv=rounded_sites, parent_index=np.array(parents,dtype=np.int32),macro_vertices_xy=np.array(macros,dtype=float),microtri_vertices_xy=np.array(microtriangles,dtype=float))
    manifest={
       'schema':'npms-fixed-fan-basis-v1','valence':n,'degree':5,'macro_count':n,
       'microcells_per_macro':12,'bernstein_coefficients_per_microcell':21,
       'function_count':len(basis), 'coarse_parent_count':len(base),
       'macro_vertex_order':'[center=(0,0),outer_i,outer_i_plus_1]',
       'microcell_index_order':'base_v17/authentic_ps12.py cells',
       'bernstein_order_ijk':triples,
       'basis_layout':'coefficients[function,macro,cell,bernstein_ijk]',
       'sites_layout':'sites_uv[function,2] absolute parameter coordinates',
       'parents_layout':'parent_index[function] coarse v1.7 parent ID; inherited site',
       'coordinate_mapping':'reference local v1=(0,3) maps to O; v2=(0,0) maps to outer_i; v3=(3,0) maps to outer_i_plus_1',
       'macro_vertices_xy':macros,'reference_vertices':reference_vertices,
       'reference_cells':[[str(x) for x in cell] for cell in cells],
       'controls':control_groups,
       'source':'full_positive_results.json + base_v17/verify_v17.py + base_v17/physical_nodal.py',
       'caution':'Demo ONLY, two fixed fans; complete P5 basis ill-conditioned, duplicate sites; not general topology/refinement-stable.'}
    (OUT/'fixtures'/f'fan{n}_schema.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('EXPORTED', n, 'shape', coeff.shape,'min',coeff.min(),'max',coeff.max(),'bytes',npzpath.stat().st_size,flush=True)
    return coeff,manifest

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--valence',type=int,choices=[5,9],default=None);a=parser.parse_args()
    for n,points in [(5,points5),(9,points9)]:
        if a.valence and n!=a.valence:continue
        generate(n,points)
