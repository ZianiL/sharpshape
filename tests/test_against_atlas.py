"""Check an installation: run sharpshape on the 3 test structures (public PDB entries 1AJU, 3SKI, 4LVV) and
compare its 5 atlas pockets with the values of the RNA pocket atlas (tests/expected_atlas_pockets.csv):
geometry, composition and steps, base pairs, stacking, backbone links and secondary-structure motifs.
Usage: .venv/bin/python tests/test_against_atlas.py      (exit code 0 if everything matches)
"""
import sys
import tempfile
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from sharpshape.cli import describe_structure  # noqa: E402

expected = pd.read_csv(HERE / "expected_atlas_pockets.csv").set_index("PDB_Pocket")
rows = []
with tempfile.TemporaryDirectory() as tmp:
    for f in sorted((HERE / "data").glob("*.pdb")):
        rows += describe_structure(f, Path(tmp))
got = pd.DataFrame(rows).set_index("PDB_Pocket")
FLOAT = ["EV1", "EV2", "EV3", "rEV21", "rEV32", "rEV31", "R_gyration"]
n_bad = 0
for key, exp in expected.iterrows():
    if key not in got.index:
        print(f"{key}: MISSING"); n_bad += 1; continue
    g = got.loc[key]
    bad = [c for c in expected.columns
           if (abs(float(g[c]) - float(exp[c])) > 1e-3 * max(1.0, abs(float(exp[c]))) if c in FLOAT
               else int(g[c]) != int(exp[c]))]
    n_bad += bool(bad)
    print(f"{key}: {'identical to the atlas' if not bad else 'DIFFERS: ' + ', '.join(bad)}")
print(f"{len(expected) - n_bad}/{len(expected)} atlas pockets identical "
      f"({len(got)} pockets found in the test structures, {len(expected.columns)} descriptors compared)")
sys.exit(1 if n_bad else 0)
