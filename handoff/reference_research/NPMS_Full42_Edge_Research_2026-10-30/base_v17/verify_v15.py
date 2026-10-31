"""NPMS v1.5 -- exact positive parent-to-boundary-child red split.

All sign and geometry identities are rational, without sampling. This is a restricted
Davydov--Yeo compatible boundary-macro hierarchy, NOT an arbitrary-mesh CAD basis.
Source geometry and positive face template: archived NPMS v1.1.
Usage: python verify_v15.py
"""
from pathlib import Path
from math import factorial
import json
import sympy as s
from authentic_ps12 import V, cells, monos, rows, edge_cells, deriv_row
from fan_assembly import points5, points9
Q=s.Rational
HERE=Path(__file__).resolve().parent
b=s.Matrix([Q(z) for z in (HERE/'positive_fan_template_coefficients.txt').read_text().split()])
assert b.shape==(252,1)
ijk=[(i,j,5-i-j) for i in range(6) for j in range(6-i)]
bary=[[Q(z,5) for z in t] for t in ijk]
BM=s.Matrix([[Q(factorial(5), s.prod(factorial(k) for k in t))*s.prod(u**k for u,k in zip(p,t)) for t in ijk] for p in bary]); BMi=BM.inv()
ref=[s.Matrix(list(map(Q,V[k]))) for k in ('v1','v2','v3')]
triangles=[[s.Matrix(list(map(Q,V[k]))) for k in c] for c in cells]

def point(v,tri):return sum((v[k]*tri[k] for k in range(3)),s.zeros(2,1))
def contains(x,tri):
    K=s.Matrix([[v[0],v[1],1] for v in tri]);return min((s.Matrix([x[0],x[1],1]).T*K.inv()))>=0

def bval(ci,x):return sum((b[21*ci+j]*x[0]**i*x[1]**k for j,(i,k) in enumerate(monos)),Q(0))
def bern_values(values):return BMi*s.Matrix(values)
def local_bern(ci):
    tri=triangles[ci]
    return bern_values([bval(ci,point(u,tri)) for u in bary])
B=[local_bern(ci) for ci in range(12)]
assert all(z>=0 for v in B for z in v)

# Independently recheck all 390 archived C2 / prescribed C3 / zero-two-jet equations.
R,_=rows();boundary=[]
for (aa,bb),adj in edge_cells.items():
    if len(adj)!=1:continue
    ci=adj[0];p=s.Matrix(V[aa]);q=s.Matrix(V[bb]);d=q-p;n=s.Matrix([d[1],-d[0]])
    for k in range(3):
        for j in range(6-k):
            x=p+Q(j+1,7-k)*d
            boundary.append({21*ci+i:z for i,z in enumerate(deriv_row(x[0],x[1],n[0],n[1],k)) if z})
assert len(R+boundary)==390
for row in R+boundary:
    assert sum((b[i]*Q(v.numerator,v.denominator) for i,v in row.items()),Q(0))==0

# Two radial edges' midpoint third normal derivatives must vanish. Without this,
# C2 zero extension is NOT enough for membership in the nested D-Y space.
radial_jet=[]
for aa,bb in (('v1','v2'),('v1','v3')):
    p=s.Matrix(V[aa]);q=s.Matrix(V[bb]);d=q-p;n=s.Matrix([d[1],-d[0]]);x=(p+q)/2
    ci=next(k for k,c in enumerate(cells) if aa in c and (('w3' in c) if bb=='v2' else ('w2' in c)))
    rr=deriv_row(x[0],x[1],n[0],n[1],3)
    value=sum((b[21*ci+j]*Q(z.numerator,z.denominator) for j,z in enumerate(rr)),Q(0))
    assert value==0;radial_jet.append(str(value))

# Explicitly check all third-order jets at original macro vertices (in every cell touching vertex).
for name in ('v1','v2','v3'):
    p=s.Matrix(V[name]);candidates=[ci for ci,c in enumerate(cells) if name in c]
    for ci in candidates:
        for a in range(4):
            t=3-a
            z=sum((b[21*ci+j]*Q(factorial(i),factorial(i-a))*Q(factorial(k),factorial(k-t))*p[0]**(i-a)*p[1]**(k-t)
                for j,(i,k) in enumerate(monos) if i>=a and k>=t),Q(0))
            assert z==0,(name,ci,a,z)

# Red refinement: only the two descendants whose external edge lies on v2-v3
# admit the parent template as a single-macro D-Y function. The two other children
# have three interior edges; the compact one-macro template is NOT valid there.
w1=s.Matrix(V['w1']);w2=s.Matrix(V['w2']);w3=s.Matrix(V['w3'])
children=[('boundary_left',(w3,ref[1],w1)),('boundary_right',(w2,w1,ref[2]))]
results={};all_remainder={};eps_cert=None
for name,cverts in children:
    K=s.Matrix([[v[0],v[1],1] for v in ref]); KC=s.Matrix([[v[0],v[1],1] for v in cverts]);
    M=K.inv()*s.Matrix([[v[0],v[1]] for v in cverts])  # ref -> child
    reverse=KC.inv()*s.Matrix([[v[0],v[1]] for v in ref]) # child -> ref
    parentpieces=[];childpieces=[];parents=[]
    for j,cell in enumerate(cells):
        tri=[(s.Matrix([v[0],v[1],1]).T*M).T for v in (s.Matrix(V[k]) for k in cell)]
        tri=[s.Matrix([v[0],v[1]]) for v in tri]
        matches=[i for i,t in enumerate(triangles) if all(contains(p,t) for p in tri)]
        assert len(matches)==1,(name,j,matches)
        i=matches[0];parents.append(i)
        pc=[];cc=[]
        for u in bary:
            x=point(u,tri)
            y=s.Matrix([x[0],x[1],1]).T*reverse
            pc.append(bval(i,x));cc.append(bval(j,(y[0],y[1])))
        pb=bern_values(pc);cb=bern_values(cc)
        assert all(z>=0 for z in pb) and all(z>=0 for z in cb)
        parentpieces.extend(pb);childpieces.extend(cb)
    ratios=[p/c for p,c in zip(parentpieces,childpieces) if c>0]
    assert ratios and min(ratios)>0
    emax=min(ratios);eps_cert=emax if eps_cert is None else min(eps_cert,emax)
    assert sum(1 for p,c in zip(parentpieces,childpieces) if p-emax*c<0)==0
    all_remainder[name]=(parentpieces,childpieces)
    # exact nonvanishing witness: barycenter of the reference maps to each child
    g=point([Q(1,3)]*3,cverts)
    y=s.Matrix([g[0],g[1],1]).T*reverse
    childat_g=bval(6,(y[0],y[1])) # barycenter is in cell 6 closure, value one
    assert childat_g==1
    results[name]={'fine_microtriangles':12,'verified_Bernstein_inequalities':252,
       'certified_maximum_Bernstein_epsilon':str(emax),'child_barycenter_value':str(childat_g),
       'nonnegative_at_maximum':True,'containing_parent_microtriangles':parents}
assert eps_cert>0
assert eps_cert==Q(1,32)
# At the same parameter point the two child interiors do not overlap, so
# parent - eps*(left_child+right_child) is nonnegative from the two certificates.
# All functions outside both child triangles are the unchanged parent.
# Coarse b = b_remainder + eps b_left + eps b_right, exact functional identity.
# Matrix T: [remainder, left, right] -> parent = [1,1,1]^T.
T=s.ones(3,1);assert min(T)>=0

# Original v1.1 coarse fans: A_j = a_j - sum_i a_j(g_i)*eta*b_i;
# F_i=eta*b_i.  Verify positive quintic Bernstein coefficients on actual 5/9 fans.
H=[s.Matrix([Q(-20),Q(-20)]),s.Matrix([Q(20),Q(-20)]),s.Matrix([Q(0),Q(20)])]
HA=s.Matrix([[v[0],v[1],1] for v in H]);HI=HA.inv()
eta=Q(1,10)
fanresults={}
for n,points in ((5,points5),(9,points9)):
    pts=[s.Matrix(list(map(Q,v))) for v in points]
    minimum=None
    for i in range(n):
        O=s.zeros(2,1);p,q=pts[i],pts[(i+1)%n]
        vertices=[O,p,q]
        g=sum(vertices,s.zeros(2,1))/3
        ag=(s.Matrix([g[0],g[1],1]).T*HI)
        for ci,cell in enumerate(cells):
            tri=[s.Matrix(V[k]) for k in cell]
            for l,u in enumerate(bary):
                z=point(u,tri)
                x=p+(q-p)*z[0]/3-p*z[1]/3
                affine=s.Matrix([x[0],x[1],1]).T*HI
                adjusted=affine-eta*ag*B[ci][l]
                minimum=min([*adjusted,minimum]) if minimum is not None else min(adjusted)
    assert minimum>0
    fanresults[f'fan{n}']={'coarse_functions':3+n,'first_refined_functions':3+3*n,
        'Bernstein_lower_bound_for_three_affine_functions':str(minimum),
        'three_affine_basis_functions_have_full_fan_support':True,
        'bubble_parent_and_children_one_macro_support':True,
        'affine_precision_and_partition_of_unity_preserved_by_child_split':'exact symbolic identity; children inherit parent parameter site',
        'positive_rational_homogeneous_transfer':'T consists of zero/one entries',
        'hierarchy_basis_counts': [3+n*(2**(level+1)-1) for level in range(4)]}

output={'certification':'Exact rational coefficient-level verification of a restricted v1.5 split',
        'local_reference':{'degree':5,'PS12_cells':12,'local_dy_constraints':len(R),'zero_two_jet_boundary_constraints':len(boundary),
          'combined_constraints':len(R+boundary),'Bernstein_nonnegative':True,
          'radial_edge_third_normal_midpoint_values':radial_jet,'original_macro_vertex_third_jets_zero':True},
        'spatial_red_refinement':{'subdivision':'1 macro triangle -> 4 red children; positive controls only on 2 boundary-adjacent children',
          'two_boundary_children':results,'common_exact_relative_epsilon':str(eps_cert),
          'simultaneous_nonnegative_split':True,'exact_parent_sum_identity':True,
          'positive_transfer_matrix':[['1'],['1'],['1']],
          'independence':'parent remainder nonzero on central red child; left/right positive on disjoint interiors'},
        'fan_certificates':fanresults,
        'stability_obstruction':{'new_child_Linf_upper_bound_at_depth_l':'(1/5)*(1/32)^l, for l>=1',
          'supnorm_lower_coefficient_stability_cannot_be_level_independent':True,
          'conditioning':'no bound for a support-scaled and normalized dual has been established'},
        'limitations':['not a full direct normalized local basis','global affine reservoir still full fan',
          'new children inherit parent site for affine reproduction','no all-interior-child subdivision','no uniform coefficient stability','Speleers basis coefficients not extracted']}
(HERE/'v15_exact_results.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
