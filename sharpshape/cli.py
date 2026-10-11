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


def check_pdb_file(pdb_file):
    """Raise ValueError if the file cannot be a PDB-format structure (so that it is reported as a failure,
    not as a structure without pockets)."""
    n_atoms, first = 0, ""
    with open(pdb_file, errors="replace") as f:
        for line in f:
            first = first or line.strip()
            if line.startswith(("ATOM", "HETATM")):
                n_atoms += 1
                if len(line) < 54:
                    raise ValueError("truncated ATOM/HETATM record (coordinates missing)")
    if n_atoms == 0:
        if first.startswith("data_"):
            raise ValueError("mmCIF file: only the PDB format is supported (convert it, e.g. with gemmi or PyMOL)")
        raise ValueError("no ATOM/HETATM records: not a PDB-format structure")


def describe_structure(pdb_file, pockets_dir, cutoff=6.0, level="residue"):
    """List of descriptor rows (dicts) for the pockets of one PDB file; pocket files written in pockets_dir.
    Raises an error for an unreadable file; returns [] for a readable structure without pockets."""
    check_pdb_file(pdb_file)
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
    rows, failed, empty = [], [], []
    with tempfile.TemporaryDirectory() as tmp:
        out_dir = a.pockets_dir or Path(tmp)
        for f in a.pdb:
            try:
                new = describe_structure(f, out_dir, a.cutoff, a.level)
            except Exception as e:  # keep going with the other structures
                print(f"{f}: ERROR - {e}", file=sys.stderr)
                failed.append(f)
                continue
            if new:
                print(f"{f}: {len(new)} pocket(s)", file=sys.stderr)
            else:
                print(f"{f}: no pocket (no ligand with >= 2 nucleotides within {a.cutoff:g} A)", file=sys.stderr)
                empty.append(f)
            rows += new
    pd.DataFrame(rows).to_csv(a.output, index=False)
    n_ok = len(a.pdb) - len(failed)
    print(f"\n{len(rows)} pocket(s) written to {a.output} | {len(a.pdb)} file(s): {n_ok} read "
          f"({len(empty)} without pocket), {len(failed)} failed", file=sys.stderr)
    if failed:
        print("Failed: " + " ".join(str(f) for f in failed), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
