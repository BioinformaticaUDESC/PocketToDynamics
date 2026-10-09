"""Fetch the exact public coordinate/dictionary snapshots used in this pilot."""
from pathlib import Path
from urllib.request import urlopen
import argparse,hashlib
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root/'raw';root.mkdir(parents=True,exist_ok=True)
SOURCES={
 '1AKE.pdb':('https://files.rcsb.org/download/1AKE.pdb','651e952f55f1317f50f4e82b7f0e99053397032cd8da5b2896f79d0af9f619a9'),
 '4AKE.pdb':('https://files.rcsb.org/download/4AKE.pdb','ff798ee8791878eb58bac1c6bed32042f51b455512ac93b50bda8fa0ff0e7f78'),
 'AP5.cif':('https://files.rcsb.org/ligands/download/AP5.cif','cb4bed433f1c4ecb298ab99a75879b841a68e1727b170fd7db68043c975a6e58')}
for name,(url,digest) in SOURCES.items():
 target=root/name
 data=target.read_bytes() if target.exists() else urlopen(url,timeout=60).read()
 if hashlib.sha256(data).hexdigest()!=digest:raise SystemExit(f'Source snapshot changed: {name}. Review explicitly; do not silently reuse.')
 if not target.exists():target.write_bytes(data)
 print(name,digest)
