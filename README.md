<p align="center"><img src="docs/logo.png" alt="SharpShape logo" width="320"></p>

# sharpshape

Extract the small-molecule binding pockets of RNA structures and compute their descriptors, with the
definitions of the RNA pocket atlas (Ziani et al., in preparation) and the SHAPE eigenvalue-ratio
archetypes (Ziani et al., *Comput. Struct. Biotechnol. J.* 2026, doi:10.34133/csbj.0022). This is the
descriptor code of the SharpShape web server, as a command-line tool, for complete structures or for pockets
you have already extracted.

## Install (Linux, no root needed)

Requirements: `git`, `bash`, and Python ≥ 3.12 (`python3.12`). If Python 3.12 is missing, the script
installs [uv](https://docs.astral.sh/uv/) in `~/.local/bin` (needs `curl`), which downloads it.

```bash
git clone https://github.com/ZianiL/sharpshape.git
cd sharpshape
bash install.sh                              # ~1 GB, a few minutes
.venv/bin/python tests/test_against_atlas.py # should print 5/5 atlas pockets identical
```

Everything is installed inside the folder: `.venv/` (the `sharpshape` command, pinned versions in
`requirements-lock.txt`) and `tools/rnapolis-0.16.17/` (RNApolis 0.16.17, pinned in
`tools/requirements-rnapolis.txt`). Two environments are needed because the atlas takes stacking from
RNApolis 0.4.10 and everything else from RNApolis 0.16.17.

## Usage

```bash
.venv/bin/sharpshape structure.pdb -o pockets.csv                      # one structure
.venv/bin/sharpshape structure.pdb -o pockets.csv --pockets-dir out/   # also write each pocket as a PDB file
```

### Many structures in one command

Give as many PDB files as you like: every pocket of every structure goes into the same CSV (one row per
pocket, columns `PDB` and `Pocket` identify it).

```bash
.venv/bin/sharpshape my_structures/*.pdb -o all_pockets.csv
.venv/bin/sharpshape 1aju.pdb 3ski.pdb 4lvv.pdb -o all_pockets.csv
```

For a very long list (thousands of files, more than the shell accepts on one command line), run it in
batches of 500 with `xargs`; each batch writes its own CSV (`batch_<number>.csv`):

```bash
ls my_structures/*.pdb | xargs -n 500 sh -c '.venv/bin/sharpshape "$@" -o "batch_$$.csv"' _
```

Each file is reported on the screen: `N pocket(s)`, `no pocket` (readable structure without a ligand
pocket), or `ERROR` (file that cannot be read, e.g. mmCIF or not a structure). A summary is printed at
the end; the CSV contains the pockets of every readable file, and the exit code is 1 if any file failed.
Structures are processed one after the other (about 2-3 s per structure on a laptop).

Options: `--cutoff` (Å, default 6.0 as in the atlas), `--level residue|atom` (whole nucleotides, as in the
atlas, or only the atoms within the cutoff). Input: PDB format; the first model is used and hydrogens are ignored.

### Pockets you have already extracted

```bash
.venv/bin/sharpshape --pocket-files my_pockets/*.pdb -o pockets.csv
```

Each file is treated as one pocket: its nucleotides (A/U/G/C, and modified nucleotides recognised by their C1′
atom) are kept, heavy atoms only, one conformer per atom; a ligand, water or ions left in the file are ignored.
Size, shape, archetype, composition and steps are computed. Base pairs, stacking, backbone links and the
secondary-structure motifs need the whole structure, so they are left empty: run the tool on the complete
structure to get them. Pocket name: `<file name>_pocket`.

## Pockets

A ligand is any non-polymer, non-solvent HETATM residue that is not covalently linked to the RNA backbone
(a 5′-triphosphate G or an in-chain modified nucleotide is RNA, not a ligand). Its pocket is the set of
nucleotides with at least one heavy atom within the cutoff of a ligand heavy atom (ligand alternate
locations: atoms without altloc + altloc A). Pockets with fewer than 2 nucleotides are skipped. Modified
nucleotides are kept in the pocket and listed (`modified_nt`); the atlas itself excludes such pockets.
Pocket name: `<ligand>_<residue number>-<chain>`.

## Output columns (one row per pocket)

| Columns | Definition |
|---|---|
| `Num_Residues`, `Num_Atoms` | nucleotides and heavy atoms of the pocket |
| `EV1`-`EV3`, `rEV21`, `rEV32`, `rEV31`, `R_gyration` | eigenvalues of the covariance matrix of the atom coordinates (descending), their ratios, radius of gyration = √(EV1+EV2+EV3) |
| `archetype` | SHAPE archetype: rEV21 ≥ 0.5745 → `sphere` (rEV32 ≥ 0.5433) or `disc`; otherwise `rod` (rEV32 ≥ 0.5433) or `anisotropic`. Thresholds = medians of the atlas non-redundant subset |
| `A`, `U`, `C`, `G`, `n_modified` | nucleotide counts |
| `step_XY` (16) | ordered dinucleotide steps along real backbone links (O3′–P ≤ 2.0 Å), 5′→3′ |
| `GU` … `AA`, `UGG` … `GUC` | legacy panels: unordered dinucleotides and selected trinucleotides |
| `canonical.bp`, `non.canonical.bp` | base pairs with both nucleotides in the pocket; canonical = Watson-Crick + G-U wobble (Saenger XIX, XX, XXVIII); RNApolis 0.16.17 |
| `stacking` | stacked pairs in the pocket; RNApolis 0.4.10 |
| `covalent.bonds` | backbone links O3′–P ≤ 2.0 Å in the pocket |
| `stem`, `hairpin-loop`, `bulge`, `internal-loop`, `junction-loop`, `ss-non-loop`, `pseudoknotted` | 1 if ≥ 1 pocket nucleotide has the label (RNApolis 0.16.17, `--find-gaps --decompose-pseudoknot-free`; stems of ≥ 2 canonical pairs incl. pseudoknotted stems; loops typed by number of strands) |
| `modified_nt`, `residues` | modified nucleotides kept; pocket nucleotides (`<chain><name><number>`) |

## Test

`tests/test_against_atlas.py` runs the tool on three public PDB entries (`tests/data/`: 1AJU, 3SKI, 4LVV)
and compares their five atlas pockets with the atlas values (`tests/expected_atlas_pockets.csv`), then gives the
extracted pockets back as pocket files (`--pocket-files`) and checks size, shape and composition again.

## Third-party software

Licences as declared in each package's metadata: RNApolis (MIT), Biopython (Biopython License Agreement),
NumPy and pandas (BSD 3-Clause), PuLP (BSD) with its bundled CBC solver (Eclipse Public License 1.0).
Installed by `install.sh`, not redistributed here.

## Licence

MIT licence (see `LICENSE`). The atlas data are distributed separately under CC BY 4.0 (Zenodo record of the article).
