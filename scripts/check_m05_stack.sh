#!/usr/bin/env bash

set -euo pipefail

echo "======================================"
echo " PocketToDynamics M05 stack"
echo "======================================"

echo
echo "Python:"
python --version
command -v python

echo
echo "Vina:"
/usr/bin/vina --version
echo "/usr/bin/vina"

echo
echo "GNINA:"
/home/juan/bin/gnina --version 2>&1 | head -1
echo "/home/juan/bin/gnina"

echo
echo "Open Babel:"
/usr/bin/obabel -V
echo "/usr/bin/obabel"

echo
echo "Meeko:"
python - <<'PY'
try:
    import meeko
    print(
        getattr(
            meeko,
            "__version__",
            "installed-version-unknown",
        )
    )
    print(meeko.__file__)
except Exception as exc:
    print("NOT_FOUND")
    print(exc)
PY

echo
echo "Meeko CLIs:"
for x in \
    mk_prepare_ligand \
    mk_prepare_ligand.py \
    mk_prepare_receptor \
    mk_prepare_receptor.py
do
    printf "%-25s " "$x"
    command -v "$x" 2>/dev/null || echo "NOT_FOUND"
done
