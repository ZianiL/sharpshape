#!/usr/bin/env bash
# Install sharpshape without root rights, inside this folder:
#   .venv/                    the sharpshape command (+ RNApolis 0.4.10 for stacking), pinned versions
#   tools/rnapolis-0.16.17/   RNApolis 0.16.17 annotator (base pairs, secondary structure), pinned versions
# Needs Python >= 3.12 (python3.12 on PATH); if it is missing, uv is installed in ~/.local/bin and downloads it.
# Usage: bash install.sh        (idempotent: re-running is harmless)
set -euo pipefail
cd "$(dirname "$0")"

if command -v python3.12 >/dev/null 2>&1; then
    mkvenv() { python3.12 -m venv "$1"; "$1/bin/python" -m pip install --quiet --upgrade pip; }
    install() { local venv=$1; shift; "$venv/bin/python" -m pip install --quiet "$@"; }
else
    export PATH="$HOME/.local/bin:$PATH"
    if ! command -v uv >/dev/null 2>&1; then
        echo ">>> python3.12 not found: installing uv into ~/.local/bin to download it ..."
        curl -LsSf https://astral.sh/uv/install.sh | env UV_NO_MODIFY_PATH=1 sh
    fi
    uv python install 3.12
    mkvenv() { uv venv --quiet --python 3.12 "$1"; }
    install() { local venv=$1; shift; VIRTUAL_ENV="$venv" uv pip install --quiet "$@"; }
fi

if [ ! -x .venv/bin/sharpshape ]; then
    echo ">>> Creating .venv (sharpshape) ..."
    mkvenv .venv
    install .venv -c requirements-lock.txt -e .
fi
if [ ! -x tools/rnapolis-0.16.17/bin/annotator ]; then
    echo ">>> Creating tools/rnapolis-0.16.17 ..."
    mkvenv tools/rnapolis-0.16.17
    install tools/rnapolis-0.16.17 -r tools/requirements-rnapolis.txt
fi
echo ">>> Done. Try:  .venv/bin/sharpshape tests/data/4LVV.pdb -o pockets.csv"
echo "    Check:      .venv/bin/python tests/test_against_atlas.py"
