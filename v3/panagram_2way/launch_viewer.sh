#!/bin/bash
# Launch panagram viewer for the PvP01 vs PAM 2-way pangenome.
# Open http://127.0.0.1:8050/ in a browser once the server reports "Dash is running".
set -eu
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/venv/bin/activate"
export PATH="/home/anton/miniconda3/envs/bcfmod/bin:$PATH"
cd "$HERE/index"
exec panagram view . --ndebug --port 8050 --host 127.0.0.1 "$@"
