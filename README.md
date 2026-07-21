# WATVina

WATVina is a molecular docking engine that treats **water as a first-class
citizen**. It combines a Vina-style GA/L-BFGS conformational search with an
explicit FHFT (fluid-site + GCMC network) hydration predictor, so you can
dock against a dry receptor, predict where the hydration sites are, and let
displaceable waters bias the search — all in one self-contained binary.

The engine uses watmol as its chemical source of truth (graph, force fields,
chemical features) and separates the docking search from optional
molecular-mechanics relaxation. Everything — Eigen, pugixml, nlohmann/json,
cxxopts, watmol and all force-field parameter files — is compiled into the
binary. **No RDKit, OpenMM, or external CGenFF program is required at build
time or at run time.**

## Highlights

**Docking core**

- PDB receptors and SDF/MOL ligands, all with explicit hydrogens;
- joint multi-ligand docking, flexible receptor side chains, custom
  ROOT/scaffold definitions;
- MiniWv (`miniwv-0.3`, the only production atom-scoring model) VDW +
  directional hydrogen/halogen-bond + electrostatic scoring with
  OpenFF Sage torsion potentials (proper-torsion series + 1-4 LJ) during
  GA/L-BFGS search;
- heavy-atom RMSD clustering with `num_modes` / `energy_range` filtering;
- deterministic reruns via `--seed`, parallel search via `--cpu`.

**Water (the distinguishing part)**

- `--predict_water`: two-layer hydration prediction — an FHFT scalar oxygen
  occupancy field solved on a 0.5 Å grid, followed by an explicit rigid-water
  GCMC/GCI network that resolves orientations and mutual exclusion;
- `-w waters.pdb --water_bias_weight`: a differentiable displacement bias so
  predicted (or third-party) waters participate in docking without explicit
  bridge terms;
- outputs include an OpenDX occupancy map and a PDB of ranked hydration
  sites whose occupancy/B-factor columns encode network marginal occupancy
  and a `-kBT log(p/(1-p))` free-energy representation.

**Around the docking run**

- `--mini ligand|site|site_h`: post-docking MM relaxation with AMBER14SB/Sage or
  CHARMM36/CGenFF template force fields;
- `--genph4` / `--template`: CIF pharmacophore generation from a receptor box
  or bound pose, and template-biased docking;
- a pybind11 API (`watvina_python`) built around `DockingConfig`,
  `ParseResult` and `DockingDriver` — the same code path as the CLI.

## Installation

Pre-built release files are staged in `bin/`. They are staged compiler
outputs, not portable packages — check the ABI column before use:

| Staged file | Platform | ABI requirement |
|---|---|---|
| `bin/watvina` | Linux x86-64 | ≤ `GLIBC_2.38`, `GLIBCXX_3.4.32` |
| `bin/watvina-glibc2.17` | Linux x86-64 (old distros) | ≤ `GLIBC_2.17` |
| `bin/watvina_python.cpython-312-x86_64-linux-gnu.so` | Linux x86-64 | CPython 3.12, ≤ `GLIBC_2.38` |
| `bin/watvina.exe` | Windows x86-64 | none beyond system DLLs |
| `bin/watvina_python.pyd` | Windows x86-64 | CPython 3.12 (`python312.dll`) |

Sanity-check a CLI binary:

```bash
bin/watvina --version        # watvina 2026.7.26
bin/watvina --help_water     # FHFT option group
bin/watvina --help_sasa      # SASA surface/penalty option group
```

To build any target from source (CMake ≥ 3.18, a C++20 compiler, Python 3):

```bash
make linux_release JOBS=8                      # → build/linux/release/watvina
ctest --test-dir build/linux/release --output-on-failure
```

## Quick start (CLI)

The commands below run against the bundled `demo/1a3l` example
(cytochrome-c-domain receptor + ferulic-acid-like ligand) and finish in well
under a minute on a desktop.

**1. Dock one ligand into a defined box**

```bash
build/linux/release/watvina \
  -r demo/1a3l/rec.pdb \
  -l demo/1a3l/i.sdf \
  -c demo/1a3l/vina.in \
  --exhaustiveness 8 --num_modes 10 \
  --seed 42 --cpu 8 \
  -o docked.sdf
```

`demo/1a3l/vina.in` only carries the search box; any Vina-style config file
works and every key can be overridden on the command line:

```ini
center_x = 34.352
center_y = -10.205
center_z = 21.190
size_x   = 17.864
size_y   = 15.682
size_z   = 15.930
```

**2. Read the report.** WATVina prints a per-mode score table to stdout:

```text
global_search score summary (kcal/mol)
+--+--------+-------+-------+-------+-------+-------+-------+-------+
|No|Affinity| RMSD_b| RMSD_i|    VDW|   HXBD|  Elect|  Hydro| Strain|
+--+--------+-------+-------+-------+-------+-------+-------+-------+
| 1|   -6.99|   0.00|   0.75|  -5.14|  -1.17|  -0.56|  -0.41|   1.04|
| 2|   -6.00|   3.99|   3.91|  -4.45|  -1.12|  -0.73|  -0.37|   2.08|
| 3|   -5.62|   1.35|   1.68|  -4.68|  -1.04|  -0.55|  -0.35|   2.20|
+--+--------+-------+-------+-------+-------+-------+-------+-------+
```

| Column | Meaning |
|---|---|
| `Affinity` | Final ranking score (kcal/mol); more negative is better |
| `RMSD_b` / `RMSD_i` | Heavy-atom RMSD to the best mode / to the input pose |
| `VDW` | Weighted van der Waals term (`--weight_vdw 0.40`) |
| `HXBD` | Directional H-bond (+ metal coordination) × hbond weight + halogen-bond |
| `Elect` | Electrostatics term (`--weight_electrop 0.04`) |
| `Hydro` | Apo-SASA-gated hydrophobic term (`--weight_hydrophobic 0.012`) |
| `Strain` | Current − minimum intra+torsion energy (0 in `--score_only`) |
| `PH4` | Pharmacophore-template matching term (shown only with `--template`) |
| `Water` | Discrete-water displacement bias (shown only with `--water`) |

**3. Read the poses.** SDF output carries per-mode property tags:

| SDF tag | Meaning |
|---|---|
| `WV_ENERGY` | Reported ranking score (same as `Affinity`) |
| `WATVINA_AFFINITY` | (raw rescore − pharmacophore − water) ÷ torsion normalization |
| `WATVINA_VDW` | vdW × `weight_vdw` ÷ normalization |
| `WATVINA_HXBD` | (hbond + metal) × `weight_hbond` + halogen-bond, ÷ normalization |
| `WATVINA_ELECT` | electro × `weight_electrop` ÷ normalization |
| `WATVINA_HYDRO` | hydrophobic × `weight_hydrophobic` ÷ normalization |
| `WATVINA_STRAIN` | strain ÷ normalization |
| `WATVINA_PH4` | Pharmacophore term (shown only with `--template`) |
| `WATVINA_WATER` | Water-displacement bias (shown only with `--water`) |
| `WV_SASA_PENALTY` | SASA environment penalty (only when a SASA penalty applies) |
| `WV_MODE` / `WV_LIGAND` | Mode rank / ligand index in joint runs, 1-based |

SDF output preserves atom/residue identity and the molecule's actual partial
charges. PDBQT input/output is not supported.

## Preparing inputs

- **All hydrogens are mandatory**, including non-polar ones, on receptors
  and ligands. WATVina uses them for the directional H-bond model; a
  stripped-hydrogen input silently degrades H-bond scoring.
- **Receptor**: PDB. Flexible receptor side chains are expressed through
  PDB REMARK FLEX records inside the receptor file. Partial charges already
  present are kept; missing/zero charges fall back to the embedded templates
  (formal charge or common oxidation state for metals).
- **Ligand**: SDF or MOL. Use `--multiligs_sdf` to dock every `$$$$`-separated
  record independently (batch mode writes to `--out_dir`, default `OUTDIR`).
- **Box**: center it on the binding site, typically the centroid of a
  co-crystallized ligand, and pad each axis by ~8–10 Å beyond the ligand
  extent. Box size defaults to 25 Å per axis if unset.

## Run modes at a glance

| Mode | Command fragment | What it does |
|---|---|---|
| Global docking | *(default)* | GA global search + L-BFGS refinement |
| Score only | `--score_only` | Evaluate the input pose, no search |
| Local optimize | `--local_only [--local_steps N]` | L-BFGS from the input pose |
| Random placement | `--randomize [--randomize_steps N]` | Sample N random box poses, keep lowest objective |
| Batch | `--multiligs_sdf lib.sdf --out_dir OUTDIR` | Dock each record independently |
| MM relaxation | `--mini ligand\|site\|site_h --mini_forcefield charmm36` | Relax output poses (see `--help_relax`) |
| Pharmacophore | `--genph4` / `--template tpl.cif` | Generate/apply CIF pharmacophores; receptor SHP cutoff is tunable with `--ph4_shp_cutoff` (see `--help_ph4`) |
| Water prediction | `--predict_water` | FHFT + GCMC hydration sites (see `--help_water`) |

Every option group has its own help screen: `--help` (docking),
`--help_water`, `--help_sasa`, `--help_ph4`, `--help_relax`. The defaults printed there are
authoritative; commonly tuned ones:

| Option | Default | Guidance |
|---|---|---|
| `--exhaustiveness` | 8 | Global-search effort; raise for large/torsion-rich ligands |
| `--population` / `--ga_search` | 8 / 4 | GA population and search-size multiplier |
| `--num_modes` | 10 | Modes kept after RMSD clustering (`--rmsd 1.5`) |
| `--energy_range` | 3.0 | kcal/mol window above the best mode |
| `--cpu` | 0 | 0 = use all hardware threads |
| `--seed` | 0 | 0 = non-deterministic (entropy ^ pid ^ time); set for reproducible runs |
| `--verbosity` | 1 | 0 quiet … 3 debug |
| `--accessibility_policy` | `legacy_apo` | Named apo accessibility treatment; see `--help_sasa` |
| `--weight_hydrophobic` | 0.012 | Independent hydrophobic-channel weight |

## Water-assisted docking

Predict scalar FHFT oxygen occupancy and NMS hydration sites:

```bash
build/linux/release/watvina -r demo/1a3l/rec.pdb -c demo/1a3l/vina.in --predict_water
# → rec_fhft_occupancy.dx  (OpenDX occupancy map)
# → rec_fhft_waters.pdb    (ranked hydration sites)
```

FHFT first solves a scalar oxygen occupancy field with a radial cooperative
kernel. The protein–water external field uses a rigid TIP3P/TIP4P SO(3)
orientation integral with the same two-angle directional hydrogen-bond model
used by docking. Candidate peaks then enter an explicit rigid-water GCMC/GCI
network calculation. With the defaults, protein–water electrostatics use the
GRID-style `4r` dielectric and the C2-reduced orientation quadrature is
`8 × 12 × 6` (576 independent orientations).
The CLI reports the watmol parameterization provenance and coverage used for
the receptor, including any element-level fallback or unresolved atom count.

In the sites PDB, the occupancy column is the **network marginal occupancy**;
the B-factor column is the corresponding network free-energy representation
`-kBT log(p/(1-p))`, **not** an absolute binding free energy. Closely spaced
(< 2.3 Å) high-occupancy peaks are mutually exclusive alternative positions
for the same high-density region.

Metal receptors receive FHFT steric exclusion from common-oxidation-state
ionic radii; zero-charge metals fall back to formal charge or the common
oxidation state.

Feed predicted (or third-party, fixed-column) waters back into docking as a
single differentiable displacement bias:

```bash
build/linux/release/watvina -r demo/1a3l/rec.pdb -l demo/1a3l/i.sdf \
  -c demo/1a3l/vina.in \
  -w rec_fhft_waters.pdb --water_bias_weight 0.1 -o hydrated.sdf
```

The coupling has no explicit water-bridge term and does not rerun FHFT per
pose. `WaterBias` stays separate from the reported `Affinity`;
`--help_water` lists the authoritative defaults.

## Python quick start

The Python module is a native extension (`bin/watvina_python.*.so` /
`.pyd`), ABI-locked to CPython 3.12 on its platform. With the file on
`PYTHONPATH`:

```python
import watvina_python as wv

# parse_argv() shares the CLI's cxxopts/config merge logic — nothing drifts.
parsed = wv.parse_argv([
    "-c", "demo/1a3l/vina.in",
    "-r", "demo/1a3l/rec.pdb",
    "-l", "demo/1a3l/i.sdf",
    "--exhaustiveness", "8", "--seed", "42", "--cpu", "8",
])
if not parsed.success:
    raise RuntimeError(parsed.error_message)

driver = wv.DockingDriver.from_parse(parsed)
result = driver.run()
if not result.ok:
    raise RuntimeError(result.error_message)

for rank, pose in enumerate(result.poses, start=1):
    print(rank, round(pose.energy, 3), pose.pose.num_torsions)

with open("docked.sdf", "w") as handle:
    handle.write(driver.poses_to_string(result, "sdf", num_modes=10, energy_range=3.0))
```

Water prediction is one call away and returns structured results
(sites, occupancy grid, diagnostics, timings):

```python
parsed = wv.parse_argv(["-r", "demo/1a3l/rec.pdb", "-c", "demo/1a3l/vina.in",
                        "--predict_water"])
driver = wv.DockingDriver.from_parse(parsed)
prediction = driver.predict_water()
print(len(prediction.sites), prediction.diagnostics.converged)
print(prediction.parameterization.provenance.source)
print(prediction.parameterization.degraded,
      prediction.parameterization.atoms_unresolved)

open("waters.pdb", "w").write(wv.format_water_pdb(prediction))
open("occupancy.dx", "w").write(wv.format_water_open_dx(prediction))
```

Long-running calls (`run()`, `score_only()`, `predict_water()`,
`eval_at_pose()`) release the GIL, so multiple drivers can work from
`concurrent.futures.ThreadPoolExecutor`.

## Repository layout

```text
src/watmol/       chemistry: graph, force fields, features, embedded parameter data
src/lib/          docking engine: search, scoring, FHFT water, pharmacophore, relax, io
src/python/       pybind11 bindings (watvina_python)
src/thirdparty/   eigen5 (single shared Eigen), pugixml, nlohmann/json, cxxopts
cmake/toolchains/ MinGW cross toolchain
scripts/          force-field data builders and research utilities
demo/             runnable receptor/ligand examples (1a3l, 3ry2_biotin, akt1, dna, ...)
build/            out-of-source build trees (one per release target)
bin/              staged release artifacts
```

## Performance notes

- Search scales with `--exhaustiveness × --population × --ga_search`;
  doubling exhaustiveness roughly doubles runtime.
- `--cpu` parallelizes across search tasks, not within one task — a single
  easy ligand will not use every core.
- FHFT runtime is dominated by the external-field build (orientation
  quadrature). On the `demo/1a3l` box (~18 Å edge, 0.5 Å grid) expect tens of
  seconds; it grows with box volume and quadrature order.
- Set `--seed` for reproducible rankings; leave it 0 for production campaigns.

## Current limits

- `--mini full` is not enabled; batch input does not support `--mini`.
- TIP4P network prediction requires an explicit calibrated
  `--water_network_mu`; only the current TIP3P network has an `auto` value.
- The docking water coupling models displacement only — no explicit bridges,
  no per-pose FHFT reranking.
- The Python artifacts are ABI-specific extension modules; wheels and a
  stable cross-version Python ABI are not provided.

## License

Apache License 2.0. The original AutoDock Vina/Scripps copyright notices are
retained. WATVina is maintained by Dr. XU Ximing.
