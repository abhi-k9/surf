"""Security and functional checks for the optional local HTTP bridge."""
import json
from pathlib import Path
import subprocess
import sys
from urllib.error import HTTPError
from urllib.request import Request, urlopen
ROOT=Path(__file__).resolve().parents[1]
p=subprocess.Popen([sys.executable,str(ROOT/'viewer/server.py'),'--cli',sys.argv[1],'--port','0'],stdout=subprocess.PIPE,text=True)
try:
    line=p.stdout.readline()
    port=int(line.split('port ')[1].split('.')[0])
    base=f'http://127.0.0.1:{port}'
    def request(path,body=None,headers=None):
        r=Request(base+path,data=body,headers=headers or {})
        try:
            with urlopen(r,timeout=20) as response:return response.status,response.read()
        except HTTPError as e:return e.code,e.read()
    status,body=request('/api/meta?valence=5');assert status==200
    meta=json.loads(body);assert len(meta['supports'])==100
    scene=meta['scenes']['plane']
    payload=json.dumps({'valence':5,'scene':scene,'resolution':4,'selected':0}).encode()
    headers={'Content-Type':'application/json'}
    status,body=request('/api/mesh',payload,headers);assert status==200
    m=json.loads(body);assert len(m['triangles'])==960 and m['invalid']==0
    assert request('/api/mesh',payload,dict(headers,Origin='https://untrusted.example'))[0]==403
    assert request('/',headers={'Host':'untrusted.example'})[0]==403
    assert request('/../../etc/passwd')[0]==404
    assert request('/api/mesh',b'{}',{'Content-Type':'text/plain'})[0]==400
    assert request('/api/mesh',b'{',headers)[0]==400
    assert request('/api/mesh',b'{"valence":55}',headers)[0]==400
    scene['weights'][0]=0
    bad=json.dumps({'valence':5,'scene':scene,'resolution':4,'selected':0}).encode()
    assert request('/api/mesh',bad,headers)[0]==400
    assert request('/api/mesh',payload,headers)[0]==200
    print('PASS HTTP: valid geometry, same-origin/Host restriction, path allowlist, content type, malformed data, invalid weights, recovery')
finally:
    p.terminate()
    try:p.wait(timeout=5)
    except subprocess.TimeoutExpired:p.kill();p.wait()
