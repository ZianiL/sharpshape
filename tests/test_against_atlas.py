"""
Check an installation: run sharpshape on the 3 test structures (PDB 1AJU, 3SKI, 4LVV) and compare their
5 atlas pockets with the values of the RNA pocket atlas (expected_atlas_pockets.csv), then the same pockets
given as pocket files.

Usage:
    .venv/bin/python tests/test_against_atlas.py      (exit code 0 if everything matches)
"""

import sys
import tempfile
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from sharpshape.cli import describe_pocket_file, describe_structure  # noqa: E402

FLOAT = ["EV1", "EV2", "EV3", "rEV21", "rEV32", "rEV31", "R_gyration"]


def differences(got, expected, columns):
    """Columns whose value differs (relative tolerance 1e-3 for the real-valued descriptors)."""
    bad = []
    for c in columns:
        if c in FLOAT:
            if abs(float(got[c]) - float(expected[c])) > 1e-3 * max(1.0, abs(float(expected[c]))):
                bad.append(c)
        elif int(got[c]) != int(expected[c]):
            bad.append(c)
    return bad


def check_structures(expected):
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for f in sorted((HERE / "data").glob("*.pdb")):
            rows += describe_structure(f, Path(tmp))
    got = pd.DataFrame(rows).set_index("PDB_Pocket")

    n_bad = 0
    for key, exp in expected.iterrows():
        if key not in got.index:
            print(f"{key}: MISSING")
            n_bad += 1
            continue
        bad = differences(got.loc[key], exp, expected.columns)
        n_bad += bool(bad)
        print(f"{key}: {'identical to the atlas' if not bad else 'DIFFERS: ' + ', '.join(bad)}")
    print(f"{len(expected) - n_bad}/{len(expected)} atlas pockets identical "
          f"({len(got)} pockets found in the test structures, {len(expected.columns)} descriptors compared)")
    return n_bad


def check_pocket_files(expected):
    columns = [c for c in expected.columns
               if c in FLOAT or c in ("Num_Residues", "Num_Atoms", "A", "U", "C", "G") or c.startswith("step_")]
    n_bad = 0
    with tempfile.TemporaryDirectory() as tmp, tempfile.TemporaryDirectory() as tmp2:
        for f in sorted((HERE / "data").glob("*.pdb")):
            describe_structure(f, Path(tmp))
        for key in expected.index:
            pdb, pocket = key.split("_", 1)
            got = describe_pocket_file(Path(tmp) / "pockets" / pdb / pocket / "pocket.pdb", Path(tmp2))
            bad = differences(got, expected.loc[key], columns)
            n_bad += bool(bad)
            print(f"{key} as a pocket file: {'identical' if not bad else 'DIFFERS: ' + ', '.join(bad)}")
    print(f"{len(expected) - n_bad}/{len(expected)} pocket files identical ({len(columns)} descriptors compared)")
    return n_bad


def main():
    expected = pd.read_csv(HERE / "expected_atlas_pockets.csv").set_index("PDB_Pocket")
    n_bad = check_structures(expected) + check_pocket_files(expected)
    sys.exit(1 if n_bad else 0)


if __name__ == "__main__":
    main()
