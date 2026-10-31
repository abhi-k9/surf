#!/usr/bin/env python3
"""Fast offline fixture checksum/shape/positive/affine/jet/seam audit; no SymPy.

Rational exact certification remains in reference_research. This is a numeric
engineering validation, not independent proof of global spline guarantees.
"""
import json, hashlib
from pathlib import Path
import numpy as np
from reference_oracle import evaluate_basis,barycentric,controls,eval_sample,assertions
ROOT=Path(__file__).resolve().parents[1]

def seams(f):
    tri=f['microtri_vertices_xy'];coef=f['coefficients']
    edges={}
    for m in range(tri.shape[0]):
        for c in range(12):
            for a,b in ((0,1),(1,2),(2,0)):
                key=tuple(sorted(tuple(np.round(v,10)) for v in (tri[m,c,a],tri[m,c,b])))
                edges.setdefault(key,[]).append((m,c))
    maxabs=np.zeros(6)
    count=0
    for key,adj in edges.items():
        if len(adj)==1:continue
        assert len(adj)==2,('nonmanifold microedge',key,adj)
        vertices=[np.array(v,dtype=float) for v in key]
        for t in (.2,.5,.8):
            xy=(1-t)*vertices[0]+t*vertices[1]
            vals=[]
            for m,c in adj:
                lam=barycentric(xy,tri[m,c]); vals.append(evaluate_basis(coef[:,m,c],tri[m,c],lam))
            maxabs=np.maximum(maxabs,abs(vals[0]-vals[1]).max(axis=0))
            count+=1
    return {'edges':sum(len(v)==2 for v in edges.values()),'samples':count,'max_abs_value_grad_hessian_per_component':maxabs.tolist()}

for n in (5,9):
    root=ROOT/'fixtures'; p=root/f'fan{n}_full42_float64.npz';b=root/f'fan{n}_full42_demo_v1.bin'
    with np.load(p) as z:f={k:z[k] for k in z.files}
    assert assertions(f)
    metadata=json.loads((root/f'fan{n}_binary_spec.json').read_text())
    assert hashlib.sha256(b.read_bytes()).hexdigest()==metadata['sha256']
    assert b.stat().st_size==metadata['file_size_bytes']
    for variant in ('plane','rational'):
        xyz,w=controls(f,variant)
        samplepoint=np.array([.3,.4,.3])@f['microtri_vertices_xy'][0,0]
        ev=eval_sample(f,samplepoint,xyz,w)
        assert np.isfinite(ev['surface']).all()
    seam=seams(f)
    print('PASS',n,'controls',len(f['coefficients']), 'seam',seam)
    assert max(seam['max_abs_value_grad_hessian_per_component'])<3e-9
