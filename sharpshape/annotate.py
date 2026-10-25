"""
Base pairs, stacking, backbone links and secondary-structure motifs of the pocket nucleotides.

- Base pairs: RNApolis 0.16.17 annotator (--decompose-pseudoknot-free) on the first model, heavy atoms.
  Canonical = Watson-Crick and G-U wobble (Saenger XIX, XX, XXVIII), the others are non-canonical.
- Stem: nucleotide paired in a stem of at least 2 base pairs (pseudoknotted stems included).
- Loops: unpaired nucleotide of a hairpin, or of a loop typed by its number of strands (3 or more: junction,
  2: bulge if one strand has no unpaired nucleotide, otherwise internal loop). Loops with a strand that
  contains a chain end are not loops.
- ss-non-loop: unpaired nucleotide outside any loop. pseudoknotted: pseudoknot bracket in the dot-bracket.
- Stacking: RNApolis 0.4.10 annotator. Backbone links: O3'(i)-P(j) <= 2.0 A.

Base pairs, stacking and links are counted when both nucleotides are in the pocket.
RNApolis 0.16.17 runs from its own environment (tools/rnapolis-0.16.17, created by install.sh).
"""

import json
import os
import subprocess
import sys
import tempfile

from Bio.PDB import PDBParser

SS_CATEGORIES = ["stem", "bulge", "hairpin-loop", "internal-loop", "junction-loop", "ss-non-loop", "pseudoknotted"]
CANONICAL_SAENGER = {"XIX", "XX", "XXVIII"}
REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def is_hydrogen(line):
    el = line[76:78].strip()
    return el in ("H", "D") or (not el and line[12:16].strip().lstrip("0123456789")[:1] in ("H", "D"))


def first_model_copy(full_pdb, out_pdb):
    """First model, first occurrence of each atom (altloc column blanked), heavy atoms only."""
    out, seen = [], set()
    for line in open(full_pdb, errors="replace"):
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")) or is_hydrogen(line):
            continue
        key = (line[21], line[22:27], line[17:20], line[12:16])
        if key in seen:
            continue
        seen.add(key)
        out.append(line[:16] + " " + line[17:])
    with open(out_pdb, "w") as f:
        f.write("".join(out) + "END\n")
    return out_pdb


def find_rnapolis_annotator():
    """RNApolis 0.16.17 annotator: $SHARPSHAPE_RNAPOLIS, else tools/rnapolis-0.16.17/bin/annotator."""
    override = os.environ.get("SHARPSHAPE_RNAPOLIS")
    if override and os.path.isfile(override):
        return override
    default = os.path.join(REPO_DIR, "tools", "rnapolis-0.16.17", "bin", "annotator")
    if os.path.isfile(default):
        return default
    raise FileNotFoundError("RNApolis 0.16.17 annotator not found: run install.sh or set SHARPSHAPE_RNAPOLIS")


def run_rnapolis(full_pdb, workdir):
    stem = os.path.splitext(os.path.basename(full_pdb))[0]
    copy_pdb = first_model_copy(full_pdb, os.path.join(workdir, stem + "_firstmodel_heavy.pdb"))
    out_json = os.path.join(workdir, stem + "_rnapolis.json")
    subprocess.run([find_rnapolis_annotator(), "-f", "--decompose-pseudoknot-free", "-j", out_json, copy_pdb],
                   check=True, capture_output=True, text=True)
    with open(out_json) as f:
        return json.load(f)


def run_stacking(full_pdb, workdir):
    """Stacked nucleotide pairs from the RNApolis 0.4.10 annotator of the package environment."""
    annotator = os.path.join(os.path.dirname(sys.executable), "annotator")
    csv_file = os.path.join(workdir, "stacking.csv")
    subprocess.run([annotator, "--csv", csv_file, full_pdb], check=True, capture_output=True, text=True)
    pairs = []
    with open(csv_file) as f:
        f.readline()
        for line in f:
            p = line.strip().split(",")
            if len(p) >= 4 and p[2] == "stacking":
                pairs.append(("".join(p[0].split(".")), "".join(p[1].split("."))))
    return pairs


def parse_rnapolis(js):
    """Secondary-structure labels {resid: set of categories} and base pairs [(resid1, resid2, canonical)]."""
    def rid(nt):
        a = nt["auth"]
        return f"{a['chain']}{a['name']}{a['number']}"

    idx = {int(i): rid(v) for i, v in js["bpseq_index"].items()}
    lines = js["dot_bracket"].strip().split("\n")
    dbn = "".join(lines[i + 2] for i in range(0, len(lines) - 2, 3) if lines[i].startswith(">strand"))
    char = {idx[i]: dbn[i - 1] for i in idx if i <= len(dbn)}
    lab = {r: set() for r in char}

    for st in js["stems"] + js.get("pseudoknot_stems", []):
        if len(st["strand5p"]["sequence"]) >= 2:
            for s in (st["strand5p"], st["strand3p"]):
                for i in range(s["first"], s["last"] + 1):
                    if i in idx and char.get(idx[i], ".") != ".":
                        lab[idx[i]].add("stem")

    chain_of = {int(i): v["auth"]["chain"] for i, v in js["bpseq_index"].items()}

    def one_chain(s):
        return len({chain_of[i] for i in range(s["first"], s["last"] + 1) if i in chain_of}) == 1

    in_loop = set()
    for hp in js["hairpins"]:
        if not one_chain(hp["strand"]):
            continue
        for i in range(hp["strand"]["first"] + 1, hp["strand"]["last"]):
            if i in idx and char.get(idx[i]) == ".":
                lab[idx[i]].add("hairpin-loop")
                in_loop.add(idx[i])
    for lp in js["loops"]:
        strands = lp["strands"]
        if not all(one_chain(s) for s in strands):
            continue
        unpaired = [[idx[i] for i in range(s["first"], s["last"] + 1) if i in idx and char.get(idx[i]) == "."]
                    for s in strands]
        kind = ("junction-loop" if len(strands) >= 3 else
                "bulge" if min(len(u) for u in unpaired) == 0 else "internal-loop")
        for u in unpaired:
            for k in u:
                lab[k].add(kind)
                in_loop.add(k)

    for k, c in char.items():
        if c == "." and k not in in_loop:
            lab[k].add("ss-non-loop")
        if c not in "().":
            lab[k].add("pseudoknotted")

    pairs = {}
    for bp in js["base_pairs"]:
        pairs[frozenset((rid(bp["nt1"]), rid(bp["nt2"])))] = bp.get("saenger") in CANONICAL_SAENGER
    return lab, [(*sorted(p), c) for p, c in pairs.items() if len(p) == 2]


def covalent_pairs(model, pock_set):
    """Backbone links O3'(i)-P(j) <= 2.0 A between pocket residues."""
    o3, ph = [], []
    for chain in model:
        for res in chain:
            rid = f"{chain.id}{res.get_resname().strip()}{res.get_id()[1]}"
            if rid not in pock_set:
                continue
            if "O3'" in res:
                o3.append((rid, res["O3'"].coord))
            if "P" in res:
                ph.append((rid, res["P"].coord))
    return [(r1, r2) for r1, x1 in o3 for r2, x2 in ph if r1 != r2 and ((x1 - x2) ** 2).sum() ** 0.5 <= 2.0]


class StructureAnnotation:
    """Annotate one structure once, then describe any number of its pockets."""

    def __init__(self, full_pdb):
        with tempfile.TemporaryDirectory() as wd:
            self.labels, self.pairs = parse_rnapolis(run_rnapolis(full_pdb, wd))
            self.stacks = run_stacking(full_pdb, wd)
        self.model = PDBParser(QUIET=True).get_structure("x", full_pdb)[0]

    def describe(self, pocket_resids):
        ps = set(pocket_resids)

        def inside(pairs):
            return [p for p in pairs if p[0] in ps and p[1] in ps]

        bp = inside([p for p in self.pairs if p[0] != p[1]])
        return {
            "canonical.bp": sum(1 for p in bp if p[2]),
            "non.canonical.bp": sum(1 for p in bp if not p[2]),
            "stacking": len(inside(self.stacks)),
            "covalent.bonds": len(covalent_pairs(self.model, ps)),
            **{c: int(any(c in self.labels.get(r, ()) for r in pocket_resids)) for c in SS_CATEGORIES},
        }
