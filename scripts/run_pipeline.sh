#!/usr/bin/env bash
set -euo pipefail
cd "$HOME/PocketToDynamics"
exec "${P2D_PYTHON:-python}" scripts/run_pipeline.py "$@"
