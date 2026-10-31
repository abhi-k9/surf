#!/usr/bin/env python3
"""Record real Release timings, platform and process RSS (where supported)."""
import argparse
import json
import platform
from pathlib import Path
import subprocess
import time
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--cli',type=Path,default=ROOT/'build/release/npms_cli');a=p.parse_args()
info={'platform':platform.platform(),'python':platform.python_version(),'cpu':platform.processor(),'build':'Release; CMake GCC -O3 -DNDEBUG, no fast-math','compiler':subprocess.check_output(['c++','--version'],text=True).splitlines()[0]}
if Path('/proc/cpuinfo').exists():
 info['cpu']=next(x.split(':',1)[1].strip() for x in Path('/proc/cpuinfo').read_text().splitlines() if x.startswith('model name'))
(ROOT/'evidence').mkdir(exist_ok=True)
for n in (5,9):
 start=time.perf_counter();r=subprocess.run([str(a.cli),str(ROOT/'handoff'),str(n),'benchmark'],text=True,capture_output=True,check=True);result=json.loads(r.stdout);result['wall_seconds']=time.perf_counter()-start;result['hardware']=info
 try:
  import resource
  result['peak_child_rss_kib']=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
 except ImportError:
  result['peak_child_rss_kib']=None
 (ROOT/'evidence'/f'fan{n}-benchmark.json').write_text(json.dumps(result,indent=2))
 print(json.dumps(result),flush=True)
