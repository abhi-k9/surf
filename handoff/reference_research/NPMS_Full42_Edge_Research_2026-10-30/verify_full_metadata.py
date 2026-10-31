"""Fast independent rational structure check on committed full-space positive split.
Does NOT re-run the expensive Bernstein extraction; use full_space_positive.py for that.
"""
from pathlib import Path
from collections import Counter
import json
import sympy as s
P=Path(__file__).resolve().parent
D=json.loads((P/'full_positive_results.json').read_text());Q=s.Rational
for n in (5,9):
 r=D['fan'+str(n)];modes=r['modes'];counts=Counter(x['parent'] for x in modes)
 assert r['original_positive_basis']==3*(n+1)
 assert len(modes)==15*n+7==r['added_independent_dy_nodal_modes']
 assert r['new_full_positive_basis_count']==18*n+10==r['global_full_DY_dimension']
 assert r['max_modes_from_any_parent']==max(counts.values())
 assert min(Q(x['dominance_ratio']) for x in modes)==Q(r['min_parent_mode_dominance_ratio'])
 for x in modes:
  k=counts[x['parent']]
  assert k>=1 and Q(x['dominance_ratio'])>0
  assert Q(x['perturbation_beta'])==Q(x['dominance_ratio'])/(2*k)
  assert x['sibling_count']==k+1
 assert min(Q(x['perturbation_beta']) for x in modes)==Q(r['min_perturbation_beta'])
 assert r['exact_positive_Bernstein'] and r['exact_partition_unity']
 print('PASS fast rational witness structure, fan',n,'dimensions',r['global_full_DY_dimension'],'children',r['new_full_positive_basis_count'])
