#!/usr/bin/env python3
"""Portable little-endian fixture packer, no SymPy required.

The binary is demonstrator data, not a universal tessellation format.
"""
import json,struct, hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
for n in (5,9):
    path=ROOT/'fixtures'/f'fan{n}_full42_float64.npz'
    with np.load(path) as z: a={k:z[k] for k in z.files}
    coeff=a['coefficients']
    site=a['sites_uv']; parent=a['parent_index']; nf,nm,nc,nb=coeff.shape
    assert nm==n and nc==12 and nb==21 and len(site)==nf and len(parent)==nf
    coarse_count=int(max(parent))+1
    coarse=np.zeros((coarse_count,n,12,21),dtype='<f8')
    coarse_sites=np.zeros((coarse_count,2),dtype='<f8')
    seen=set()
    for i,p in enumerate(parent):
        coarse[p]+=coeff[i]
        if int(p) not in seen:
            coarse_sites[p]=site[i]
            seen.add(int(p))
        else: assert np.array_equal(coarse_sites[p],site[i])
    assert len(seen)==coarse_count
    assert np.max(abs(coarse.sum(axis=0)-1))<1e-12
    f64=lambda x: np.ascontiguousarray(x,dtype='<f8').tobytes()
    u32=lambda x: np.ascontiguousarray(x,dtype='<u4').tobytes()
    header=struct.pack('<8s8I',b'NPMSFAN1',1,5,n,n,12,21,nf,coarse_count)
    data=b''.join([header,
       f64(a['macro_vertices_xy']),
       f64(a['microtri_vertices_xy']),
       f64(site),
       u32(parent),
       f64(coeff),
       f64(coarse_sites),
       f64(coarse)])
    dest=ROOT/'fixtures'/f'fan{n}_full42_demo_v1.bin'
    dest.write_bytes(data)
    full={'schema':'npms-fixed-fan-binary-v1','magic_ascii':'NPMSFAN1','endian':'little','float':'IEEE754 binary64','uint':'uint32',
      'header':'8-byte magic then eight little-endian uint32 values in order version,degree,valence,macro_count,microcells_per_macro,bernstein_coefficients_per_cell,fine_count,coarse_count',
      'header_size':40,
      'arrays_in_exact_order':[
          {'name':'macro_vertices_xy','type':'float64','shape':[n,3,2]},
          {'name':'microtri_vertices_xy','type':'float64','shape':[n,12,3,2]},
          {'name':'fine_sites_uv','type':'float64','shape':[nf,2]},
          {'name':'fine_parent_index','type':'uint32','shape':[nf]},
          {'name':'fine_coefficients','type':'float64','shape':[nf,n,12,21]},
          {'name':'coarse_sites_uv','type':'float64','shape':[coarse_count,2]},
          {'name':'coarse_coefficients','type':'float64','shape':[coarse_count,n,12,21]}],
      'coefficient_order':'ijk enumerated in fan{n}_schema.json: i=0..5, j=0..5-i, k=5-i-j',
      'float64_fixture_precision':'Generated from certified exact rational construction but truncated to double. Runtime checks should use numerical tolerances; original research handles exact proofs.',
      'file_size_bytes':len(data), 'sha256':hashlib.sha256(data).hexdigest(),
      'coarse_transfer':'Fine parent_index copies old control homogeneous (wP,w); coarse[pid]=sum fine coefficients of that pid at identical parametric microcell.',
      'restrictions':'Only original fixed valence-5/9 fan geometries, no generalized arbitrary-topology support.'}
    (ROOT/'fixtures'/f'fan{n}_binary_spec.json').write_text(json.dumps(full,indent=2)+'\n')
    # check direct round trip
    blob=dest.read_bytes()
    hdr=struct.unpack_from('<8s8I',blob)
    assert hdr==(b'NPMSFAN1',1,5,n,n,12,21,nf,coarse_count)
    check=np.frombuffer(blob,dtype='<f8',count=coeff.size,offset=40+a['macro_vertices_xy'].size*8+a['microtri_vertices_xy'].size*8+site.size*8+parent.size*4)
    assert np.array_equal(check.reshape(coeff.shape),coeff)
    print('PASS',dest.name,'size',len(data),'sha256',full['sha256'])
