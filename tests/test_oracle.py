"""Independent NPZ oracle versus live C++ worker, without regenerating goldens."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
HANDOFF = ROOT / 'handoff'
spec = importlib.util.spec_from_file_location('oracle', HANDOFF / 'tools/reference_oracle.py')
oracle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(oracle)

def run(cli):
    for n in (5, 9):
        with np.load(HANDOFF / 'fixtures' / f'fan{n}_full42_float64.npz', allow_pickle=False) as data:
            f = {k: data[k] for k in data.files}
        oracle.assertions(f)
        golden = json.loads((HANDOFF / 'fixtures' / f'fan{n}_golden_samples.json').read_text())
        requests, expected, basis_expected = [], [], []
        for variant in ('plane', 'rational'):
            scene = json.loads((HANDOFF / 'scenes' / f'fan{n}_{variant}.json').read_text())
            for row in golden['samples']:
                requests.append({'scene': scene, 'xy': row['xy'], 'basis': True})
                evaluation = oracle.eval_sample(f, row['xy'], scene['control_xyz'], scene['weights'])
                value = evaluation['surface']
                basis_expected.append(evaluation['basis'])
                assert np.max(np.abs(value - row[variant]['surface_jet_6x3'])) < 1e-9
                expected.append(value)
        result = subprocess.run([cli, str(HANDOFF), str(n), 'worker'],
                                input=''.join(json.dumps(r) + '\n' for r in requests),
                                text=True, capture_output=True, check=True, timeout=90)
        lines = result.stdout.splitlines()
        assert len(lines) == len(expected)
        error = max(float(np.max(np.abs(np.array(json.loads(line)['jet']) - target)))
                    for line, target in zip(lines, expected))
        basis_error = max(float(np.max(np.abs(np.array(json.loads(line)['basis']) - target))) for line,target in zip(lines,basis_expected))
        assert basis_error < 1e-9, basis_error
        assert error < 1e-9, error
        print(f'PASS independent NPZ/Python/C++ fan={n} scenes=2 samples={len(expected)} max_jet_error={error:.6g} max_all_basis_jet_error={basis_error:.6g}')
    # Check all original assets, including untouched goldens after execution.
    for line in (HANDOFF / 'MANIFEST_SHA256.txt').read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        assert hashlib.sha256((HANDOFF / name).read_bytes()).hexdigest() == digest, name
    print('PASS original manifest intact')

if __name__ == '__main__':
    run(sys.argv[1])
