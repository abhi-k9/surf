import json
from pathlib import Path
import subprocess
import sys
ROOT = Path(__file__).resolve().parents[1]
scene = json.loads((ROOT / 'handoff/scenes/fan5_plane.json').read_text())
bad = []
for w in (0, -1, 'bad', None):
    s = json.loads(json.dumps(scene))
    s['weights'][0] = w
    bad.append(json.dumps({'scene': s, 'xy': [0, 0]}))
for change in ({'valence':9}, {'schema':'unknown'}, {'control_xyz':[]}, {'weights':[]}):
    bad.append(json.dumps({'scene':dict(scene, **change), 'xy':[0,0]}))
bad.extend(['{', '{"xy":[100,100]}', '{"xy":[NaN,0]}', '{"resolution":1000000000,"selected":0}', '{"resolution":8,"selected":999}', '{"resolution":0,"selected":0}', '{"resolution":2.5,"selected":0}'])
valid = json.dumps({'scene':scene, 'xy':[0,0]})
p = subprocess.run([sys.argv[1], str(ROOT / 'handoff'), '5', 'worker'],
                   input='\n'.join(bad+[valid])+'\n', text=True, capture_output=True, check=True, timeout=30)
rows = [json.loads(line) for line in p.stdout.splitlines()]
assert len(rows) == len(bad)+1
assert all('error' in row for row in rows[:-1])
assert rows[-1]['regular'] and abs(rows[-1]['jet'][0][0]) < 1e-9
print(f'PASS {len(bad)} malformed requests rejected; worker recovers for valid request')
