"""Nucleotide composition and dinucleotide / trinucleotide steps of a pocket (atlas definitions).

Steps follow real covalent backbone links, O3'(i)-P(j) <= 2.0 A between pocket nucleotides (5'->3'), not
residue numbering, so numbering gaps and chain breaks are handled. Columns: A/U/C/G counts (canonical
nucleotides), n_modified (modified nucleotides kept in the pocket), the 16 ordered steps step_XY, and the
legacy panels GU..AA (unordered dinucleotides: GU = GU + UG steps) and UGG..GUC (ordered trinucleotides).
Same code as the SharpShape web server (scripts/count_res_AUGC_combination_pocket.py).
"""
from collections import Counter, defaultdict
from itertools import product

import numpy as np

NT = {"A", "U", "G", "C"}
DI_PANEL = ["GU", "GG", "GA", "GC", "AU", "CC", "AC", "CU", "AA"]
TRI_PANEL = ["UGG", "UGC", "AUC", "AAC", "ACU", "GUC"]
STEPS = ["".join(p) for p in product("AUGC", repeat=2)]


def read_pocket(path):
    """{(chain, resnum, icode): {"name": resname, "O3'": xyz, "P": xyz}}, first conformer per atom."""
    res = {}
    for line in open(path, errors="replace"):
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[17:20].strip()
        if name in ("HOH", "WAT", "DOD"):
            continue
        k = (line[21], int(line[22:26]), line[26])
        r = res.setdefault(k, {"name": name})
        atom = line[12:16].strip()
        if atom in ("O3'", "P") and atom not in r:
            r[atom] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return res


def composition(pocket_pdb):
    res = read_pocket(pocket_pdb)
    base = {k: v["name"] for k, v in res.items()}
    links = [(i, j) for i in res for j in res
             if i != j and i[0] == j[0] and "O3'" in res[i] and "P" in res[j]
             and np.linalg.norm(res[i]["O3'"] - res[j]["P"]) <= 2.0]
    nxt = defaultdict(list)
    for i, j in links:
        nxt[i].append(j)
    canon_links = [(i, j) for i, j in links if base[i] in NT and base[j] in NT]
    steps = Counter(base[i] + base[j] for i, j in canon_links)
    tri = Counter(base[i] + base[j] + base[k] for i, j in canon_links for k in nxt[j] if base[k] in NT)
    unordered = Counter()
    for s, n in steps.items():
        unordered["".join(sorted(s))] += n
    nuc = Counter(b for b in base.values() if b in NT)
    return {
        **{n: nuc.get(n, 0) for n in "AUCG"},
        "n_modified": sum(b not in NT for b in base.values()),
        **{f"step_{s}": steps.get(s, 0) for s in STEPS},
        **{m: unordered.get("".join(sorted(m)), 0) for m in DI_PANEL},
        **{m: tri.get(m, 0) for m in TRI_PANEL},
    }
