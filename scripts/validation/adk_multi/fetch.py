#!/usr/bin/env python3
"""Fetch pinned RCSB PDB and CCD snapshots for the ADP redocking pilot."""
from pathlib import Path
from urllib.request import urlopen
import argparse, hashlib
SOURCES = {
    '8RJ4.pdb': ('https://files.rcsb.org/download/8RJ4.pdb', '4d1c9d763c99d84728738b99a0d4a100ea664e72e793f90be559dc360a17939d'),
    'ADP.cif': ('https://files.rcsb.org/ligands/download/ADP.cif', 'fa9eabb464c85c7c5a44146d492b907acd3151a250a4598001cc6abba8da84af'),
}
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args()
    out=a.root/'raw';out.mkdir(parents=True,exist_ok=True)
    for name,(url,digest) in SOURCES.items():
        target=out/name;data=target.read_bytes() if target.exists() else urlopen(url,timeout=60).read()
        got=hashlib.sha256(data).hexdigest()
        if got != digest: raise SystemExit(f'{name}: expected {digest}, received {got}; review the source revision')
        if not target.exists(): target.write_bytes(data)
        print(f'{name}\t{got}')
if __name__=='__main__':main()
