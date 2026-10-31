#!/usr/bin/env python3
"""Create a deterministic source/evidence delivery ZIP; no build or network actions."""
import hashlib
from pathlib import Path
import zipfile
ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT.parent / 'NPMS_Prototype_v0.1.zip'
def files():
    return sorted(p for p in ROOT.rglob('*') if p.is_file() and
                  not any(x in {'.tools', 'build', '.git', '__pycache__'} for x in p.relative_to(ROOT).parts)
                  and p.name != 'DELIVERY_SHA256.txt')
paths = files()
manifest = ROOT / 'DELIVERY_SHA256.txt'
manifest.write_text(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT).as_posix()}\n' for p in paths))
with zipfile.ZipFile(DEST, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
    for p in sorted(paths + [manifest]):
        info = zipfile.ZipInfo('NPMS_Prototype/' + p.relative_to(ROOT).as_posix(), date_time=(2026,10,30,0,0,0))
        info.compress_type = zipfile.ZIP_DEFLATED
        info.external_attr = (0o100755 if p.name == 'npms_cli' else 0o100644) << 16
        archive.writestr(info, p.read_bytes(), compresslevel=9)
print(f'{DEST} {DEST.stat().st_size} bytes SHA256={hashlib.sha256(DEST.read_bytes()).hexdigest()}')
