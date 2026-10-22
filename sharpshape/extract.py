"""Ligand-guided RNA pocket extraction (reproducible version of the atlas protocol).

Same logic as Data_generation/.../2ndversion_LIG-RESOLV/extract_residprox_allresidtypes_noH_conf-altoptimal_lig-resolv.py
(ligand = non-polymer HETATM residue, pocket = whole residues with >= 1 heavy atom within 6 A of the
ligand, first model, hydrogens removed), with the errors found during the October 2026 check fixed:

  F1  a HETATM residue covalently linked to the RNA backbone (O3'-P <= 2.0 A with a nucleotide) is
      part of the RNA, not a ligand (5'-triphosphate G declared as a separate entity, in-chain modified
      nucleotides numbered differently from the mmCIF polymer sequence)
  F2  modified nucleotides inside the 6 A shell are detected explicitly instead of being dropped with
      the HETATM lines. Paper mode excludes the pocket, server mode keeps them in the pocket
  F3  ligand alternate locations: atoms without altloc + altloc A (or the first altloc), instead of
      altloc A only

Modes
  paper   ligands restricted to the HARIBOSS compound list (--ligand-list), OXV/Y59 removed (no SMILES),
          pockets containing a modified nucleotide excluded, >= 2 canonical nucleotides
  server  any non-polymer, non-solvent, non-covalent ligand, modified nucleotides kept in the pocket
          (flagged), >= 2 nucleotides, mmCIF optional (polymer detected from ATOM records + covalent links)

Usage
  python3 pipeline/extract_pockets.py --pdb-dir DIR --cif-dir DIR --out DIR --mode paper --ligand-list FILE
Outputs (in --out)
  pockets/<pdb>/<LIG>_<num>-<chain>/pocket.pdb   RNA residues of the pocket (one conformer per atom:
                                                 highest occupancy, first on ties, as Biopython reads it)
  pockets/<pdb>/<LIG>_<num>-<chain>/resid.txt    nucleotide IDs <chain><resname><resnum>, sorted
  pockets_index.csv   one row per ligand instance: status (kept / reason excluded), counts, geometry
"""
import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

NT = {"A", "U", "G", "C"}
SOLVENT = {"HOH", "WAT", "DOD", "NA", "K", "CL", "CA", "MG", "MN", "ZN", "CU", "FE", "CO", "NI", "CD",
           "SR", "BA", "CS", "SO4", "PO4", "NO3", "CO3", "IOD", "BR"}  # identical to the original script
NO_SMILES = {"OXV", "Y59"}
CUTOFF = 6.0
O3P_MAX = 2.0


def is_hydrogen(line):
    el = line[76:78].strip()
    return el in ("H", "D") or (not el and line[12:16].strip().lstrip("0123456789")[:1] in ("H", "D"))


def read_first_model(path):
    """{residue key: {"rec": ATOM/HETATM, "atoms": [(name, altloc, occupancy, xyz, line)]}}"""
    res = {}
    for line in open(path, errors="replace"):
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")) or is_hydrogen(line):
            continue
        key = (line[21], int(line[22:26]), line[26].strip(), line[17:20].strip())
        r = res.setdefault(key, {"rec": line[:6].strip(), "atoms": []})
        occ = float(line[54:60]) if line[54:60].strip() else 1.0
        r["atoms"].append((line[12:16].strip(), line[16], occ,
                           np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]), line))
    return res


def mmcif_polymer(cif):
    if cif is None or not Path(cif).exists():
        return set()
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    c = MMCIF2Dict(str(cif))
    out = set()
    for ch, n, ic, m in zip(c["_pdbx_poly_seq_scheme.pdb_strand_id"], c["_pdbx_poly_seq_scheme.pdb_seq_num"],
                            c["_pdbx_poly_seq_scheme.pdb_ins_code"], c["_pdbx_poly_seq_scheme.mon_id"]):
        if n in ("?", "."):
            continue
        out.add((ch.strip(), int(n), "" if ic in (".", "?", " ") else ic.strip(), m.strip()))
    return out


def atom_xyz(r, name):
    for n, alt, _, xyz, _ in r["atoms"]:
        if n == name and alt in (" ", "A", "1"):
            return xyz
    return None


def covalent_nucleotides(res):
    """Nucleotide-like residues (A/U/G/C or with a C1' atom) joined to another one by O3'-P <= 2.0 A."""
    nucl = [k for k, r in res.items() if k[3] not in SOLVENT and (k[3] in NT or atom_xyz(r, "C1'") is not None)]
    o3 = [(k, atom_xyz(res[k], "O3'")) for k in nucl]
    ph = [(k, atom_xyz(res[k], "P")) for k in nucl]
    o3 = [(k, v) for k, v in o3 if v is not None]
    ph = [(k, v) for k, v in ph if v is not None]
    linked = set()
    if o3 and ph:
        d = np.linalg.norm(np.array([v for _, v in o3])[:, None] - np.array([v for _, v in ph])[None], axis=2)
        for i, j in zip(*np.where(d <= O3P_MAX)):
            if o3[i][0] != ph[j][0]:  # a 3'-phosphate on the same residue is not a link
                linked.update((o3[i][0], ph[j][0]))
    return linked


def one_conformer(atoms):
    """One record per atom name: highest occupancy, first on ties (Biopython's behaviour)."""
    best = {}
    for a in atoms:
        if a[0] not in best or a[2] > best[a[0]][2]:
            best[a[0]] = a
    return list(best.values())


def geometry(xyz):
    ev = np.sort(np.linalg.eigvalsh(np.cov(xyz, rowvar=False)))[::-1]  # n-1 covariance, as originally
    return {"EV1": ev[0], "EV2": ev[1], "EV3": ev[2], "rEV21": ev[1] / ev[0], "rEV32": ev[2] / ev[1],
            "rEV31": ev[2] / ev[0], "R_gyration": float(np.sqrt(ev.sum())), "Num_Atoms": len(xyz)}


def process_structure(pdb_file, cif_file, mode, ligand_list, out_dir, cutoff=CUTOFF, level="residue"):
    """level='residue': whole residues (paper), level='atom': only the atoms within the cutoff
    (server option). The cutoff is 6 A for the paper, and the server allows 4-7 A."""
    pdb = pdb_file.stem.lower()
    res = read_first_model(pdb_file)
    linked = covalent_nucleotides(res)
    mmcif_poly = mmcif_polymer(cif_file)
    polymer = mmcif_poly | {k for k, r in res.items() if r["rec"] == "ATOM"} | linked
    modified = {k for k in polymer if k in res and k[3] not in NT and k[3] not in SOLVENT
                and (k in linked or atom_xyz(res[k], "C1'") is not None or res[k]["rec"] == "ATOM")}
    keys = list(res)
    all_atoms = [a for k in keys for a in res[k]["atoms"]]
    all_xyz = np.array([a[3] for a in all_atoms])
    owner = [k for k in keys for _ in res[k]["atoms"]]
    rows = []
    for k, r in res.items():
        if r["rec"] != "HETATM" or k[3] in SOLVENT:
            continue
        ch, num, ic, lig = k
        name = f"{lig}_{num}-{ch}"
        row = {"PDB_Pocket": f"{pdb}_{name}", "PDB": pdb, "Pocket": name, "LIG": lig}
        if mode == "paper" and (lig not in ligand_list or lig in NO_SMILES):
            continue  # not a HARIBOSS compound: not a candidate pocket
        if k in polymer:
            if k in linked and k not in mmcif_poly:  # recorded: named as a ligand but part of the RNA (rule F1)
                rows.append({**row, "status": "not a ligand: covalently linked to the RNA backbone"})
            continue
        alts = sorted({a[1] for a in r["atoms"]} - {" "})
        keep_alt = "A" if "A" in alts else (alts[0] if alts else " ")
        lx = np.array([a[3] for a in r["atoms"] if a[1] in (" ", keep_alt)])
        dmin = np.min(np.linalg.norm(all_xyz[:, None] - lx[None], axis=2), axis=1)
        close = np.where(dmin <= cutoff)[0]
        shell = {owner[i] for i in close} - {k}
        canon = sorted(s for s in shell if s[3] in NT and s in polymer)
        mods = sorted(s for s in shell if s in modified)
        row.update({"n_canonical_nt": len(canon), "n_modified_nt": len(mods),
                    "modified_nt": ";".join(f"{s[0]}{s[3]}{s[1]}" for s in mods)})
        members = canon + (mods if mode == "server" else [])
        if mode == "paper" and mods:
            row["status"] = "excluded: modified nucleotide within 6 A"
        elif len(members) < 2:
            row["status"] = "excluded: fewer than 2 nucleotides"
        else:
            row["status"] = "kept"
            atoms = [a for s in members for a in one_conformer(res[s]["atoms"])]
            if level == "atom":
                near = {id(all_atoms[i]) for i in close}
                atoms = [a for a in atoms if id(a) in near]
            row.update(geometry(np.array([a[3] for a in atoms])))
            row["Num_Residues"] = len(members)
            d = out_dir / "pockets" / pdb / name
            d.mkdir(parents=True, exist_ok=True)
            with open(d / "pocket.pdb", "w") as f:
                f.writelines(a[4] for a in atoms)
                f.write("END\n")
            with open(d / "resid.txt", "w") as f:
                f.writelines(f"{s[0]}{s[3]}{s[1]}\n" for s in sorted(members, key=lambda s: f"{s[0]}{s[3]}{s[1]}"))
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb-dir", required=True, type=Path)
    ap.add_argument("--cif-dir", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--mode", choices=["paper", "server"], default="paper")
    ap.add_argument("--ligand-list", type=Path, help="HARIBOSS compound IDs, one per line (paper mode)")
    ap.add_argument("--cutoff", type=float, default=CUTOFF, help="distance cutoff in A (paper: 6.0)")
    ap.add_argument("--level", choices=["residue", "atom"], default="residue")
    args = ap.parse_args()
    ligs = {l.strip() for l in open(args.ligand_list) if l.strip()} if args.ligand_list else set()
    if args.mode == "paper" and not ligs:
        raise SystemExit("paper mode needs --ligand-list")
    rows = []
    for pdb_file in sorted(args.pdb_dir.glob("*.pdb")):
        cif = args.cif_dir / f"{pdb_file.stem.lower()}.cif" if args.cif_dir else None
        rows += process_structure(pdb_file, cif, args.mode, ligs, args.out, args.cutoff, args.level)
    idx = pd.DataFrame(rows)
    idx.to_csv(args.out / "pockets_index.csv", index=False)
    print(idx.status.value_counts().to_string())


if __name__ == "__main__":
    main()
