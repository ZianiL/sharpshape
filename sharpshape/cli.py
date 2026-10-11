"""sharpshape: extract the small-molecule binding pockets of RNA structures and compute their descriptors.

For each ligand (non-polymer, non-solvent HETATM residue not covalently linked to the RNA), the pocket is the
set of nucleotides with at least one heavy atom within the cutoff (default 6 A) of the ligand; modified
nucleotides are kept and listed. One row per pocket: geometry (eigenvalues, ratios, radius of gyration),
SHAPE archetype, composition and steps, base pairs, stacking, backbone links and secondary-structure motifs.

  sharpshape structure.pdb [more.pdb ...] -o pockets.csv [--pockets-dir DIR] [--cutoff 6] [--level residue]
"""
import argparse
import sys
import tempfile
from pathlib import Path

import pandas as pd

from . import extract
from .annotate import StructureAnnotation
from .composition import composition
from .shape import classify_shape

GEOMETRY = ["Num_Residues", "Num_Atoms", "EV1", "EV2", "EV3", "rEV21", "rEV32", "rEV31", "R_gyration"]


def describe_structure(pdb_file, pockets_dir, cutoff=6.0, level="residue"):
    """List of descriptor rows (dicts) for the pockets of one PDB file; pocket files written in pockets_dir."""
    rows = []
    kept = [r for r in extract.process_structure(Path(pdb_file), None, "server", set(), Path(pockets_dir),
                                                 cutoff, level) if r.get("status") == "kept"]
    if not kept:
        return rows
    ann = StructureAnnotation(str(pdb_file))
    for r in kept:
        d = Path(pockets_dir) / "pockets" / r["PDB"] / r["Pocket"]
        resids = [l.strip() for l in open(d / "resid.txt") if l.strip()]
        rows.append({
            "PDB_Pocket": r["PDB_Pocket"], "PDB": r["PDB"], "Pocket": r["Pocket"], "LIG": r["LIG"],
            **{k: r[k] for k in GEOMETRY},
            "archetype": classify_shape(r["rEV21"], r["rEV32"]),
            **composition(d / "pocket.pdb"),
            **ann.describe(resids),
            "modified_nt": r.get("modified_nt", ""),
            "residues": ";".join(resids),
        })
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(prog="sharpshape", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdb", nargs="+", type=Path, help="structure file(s) in PDB format")
    ap.add_argument("-o", "--output", type=Path, default=Path("pockets.csv"), help="output CSV (default pockets.csv)")
    ap.add_argument("--pockets-dir", type=Path, help="also write each pocket as DIR/pockets/<pdb>/<pocket>/pocket.pdb")
    ap.add_argument("--cutoff", type=float, default=6.0, help="distance cutoff in A (atlas: 6.0)")
    ap.add_argument("--level", choices=["residue", "atom"], default="residue",
                    help="residue: whole nucleotides (atlas); atom: only atoms within the cutoff")
    a = ap.parse_args(argv)
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        out_dir = a.pockets_dir or Path(tmp)
        for f in a.pdb:
            try:
                new = describe_structure(f, out_dir, a.cutoff, a.level)
            except Exception as e:  # keep going with the other structures
                print(f"{f}: failed ({e})", file=sys.stderr)
                continue
            print(f"{f}: {len(new)} pocket(s)", file=sys.stderr)
            rows += new
    pd.DataFrame(rows).to_csv(a.output, index=False)
    print(f"{len(rows)} pocket(s) written to {a.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
