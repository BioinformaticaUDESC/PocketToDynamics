#!/usr/bin/env bash

set -euo pipefail

ROOT="$HOME/PocketToDynamics"

cd "$ROOT"

echo "======================================"
echo " PocketToDynamics"
echo "======================================"

echo
echo "[M01] Project initialization"
python scripts/01_init_project.py

echo
echo "[M02] Input normalization"
python scripts/02_normalize_inputs.py

echo
echo "[M03] Pocket mapping"
python scripts/03_map_pockets.py

echo
echo "[M04] Job matrix builder"
python scripts/04_build_job_matrix.py

echo
echo "======================================"
echo " Pipeline preparation completed"
echo "======================================"
