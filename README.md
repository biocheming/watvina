# WATVina 2026.7.20 Release Guide

![](https://github.com/biocheming/watvina/blob/main/watvina_logo.png)

This directory is a self-contained x86-64 binary release. WATVina embeds its chemical-perception and force-field data. Receptor, ligand and configuration files are supplied by the user.

## Choose The Correct File

| File | Platform and ABI |
|---|---|
| `watvina` | Linux x86-64 CLI; maximum `GLIBC_2.38`, `GLIBCXX_3.4.32` |
| `watvina-glibc2.17` | Linux x86-64 CLI; verified maximum `GLIBC_2.17` |
| `watvina.exe` | Windows x86-64 CLI |
| `watvina_python.cpython-312-x86_64-linux-gnu.so` | Linux x86-64 CPython 3.12 extension |
| `watvina_python.pyd` | Windows x86-64 CPython 3.12 extension; requires `python312.dll` |
| `watvina-python-guide.md` | Complete Python API and workflow guide |

The normal Linux binary is preferred on a current distribution. Use `watvina-glibc2.17` when the normal binary reports that a `GLIBC` or `GLIBCXX` version is unavailable. The Python extensions are ABI-specific: the `.so` and `.pyd` are not interchangeable and cannot be imported by Python 3.11 or 3.13.

Check the executable before starting:

```bash
chmod +x watvina watvina-glibc2.17
./watvina --version
./watvina --help
```

On Windows PowerShell:

```powershell
.\watvina.exe --version
.\watvina.exe --help
```

## Prepare Inputs

WATVina accepts a receptor in PDB or PDBQT and a ligand in SDF or PDBQT. Both receptor and ligand must contain all explicit hydrogens, including non-polar hydrogens. Bond orders, formal charges and protonation states must be chemically correct before docking; the program does not silently repair an ambiguous chemical identity.

Define the search box in a Vina-style text file such as `vina.in`:

```ini
center_x = 17.3635
center_y = 6.1955
center_z = 12.6360
size_x   = 15.393
size_y   = 18.589
size_z   = 18.688
```

The center is normally the centroid of a known ligand or the pocket center. The box should include the ligand plus enough room to translate and rotate; padding each axis by roughly 8-10 A beyond the bound ligand is a useful starting point.

## Dock A Ligand

```bash
./watvina \
  -c vina.in \
  -r receptor.pdb \
  -l ligand.sdf \
  --exhaustiveness 8 \
  --num_modes 10 \
  --seed 42 \
  --cpu 8 \
  -o docked.sdf
```

`--seed 42` makes repeated runs reproducible. `--cpu 0` uses all hardware threads. Increase `--exhaustiveness` for large or torsion-rich ligands.

The terminal table reports the final ranking score (`Affinity`) and its VDW, directional hydrogen/halogen bond (`HXBD`), electrostatic, intramolecular, torsion, pharmacophore and water components. More negative `Affinity` is better within the same scoring setup; it is not an experimental binding free energy.

The SDF contains one record per retained mode and includes:

| Property | Meaning |
|---|---|
| `WV_ENERGY` | Reported ranking score, equal to `Affinity` |
| `WV_RAW_ENERGY` | Raw scoring-kernel total before final post-processing |
| `WV_WATER_BIAS` | Discrete-water displacement bias |
| `WV_MODE` | One-based mode rank |
| `WV_LIGAND` | One-based ligand index in a joint multi-ligand run |

## Score Or Refine An Existing Pose

Score the coordinates exactly as supplied, without searching:

```bash
./watvina -c vina.in -r receptor.pdb -l pose.sdf --score_only
```

Optimize locally from the supplied pose:

```bash
./watvina -c vina.in -r receptor.pdb -l pose.sdf \
  --local_only --local_steps 80 -o local.sdf
```

## Post-Docking MM Relaxation

Relax only the ligand after docking:

```bash
./watvina -c vina.in -r receptor.pdb -l ligand.sdf \
  --mini ligand --mini_forcefield charmm \
  --mini_steps 80 -o docked_relaxed.sdf
```

Use `--mini site` to relax the ligand and supported flexible pocket atoms. The CHARMM route applies biopolymer templates first, two-hop CGenFF chemical perception to non-standard chemistry, and miniCGenFF as the final organic small-molecule fallback. Metal complexes are handled separately. Read the printed parameterization table: unresolved atoms or bonded terms are fatal; fallback and pruned-improper counts remain visible diagnostics. Run `./watvina --help_relax` for optimizer and tolerance controls.

## Batch Docking

Dock every SDF record independently:

```bash
./watvina -c vina.in -r receptor.pdb \
  --multiligs_sdf library.sdf --out_dir results \
  --exhaustiveness 8 --cpu 8
```

Batch mode creates one result per input molecule. Post-docking `--mini` is not currently supported in batch mode.

## FHFT Water Prediction

```bash
./watvina -r receptor.pdb -c vina.in --predict_water \
  --water_dx_out occupancy.dx \
  --water_sites_out waters.pdb
```

The box in `vina.in` is the prediction domain. `occupancy.dx` is the scalar oxygen occupancy field. In `waters.pdb`, the occupancy column is the explicit water-network marginal occupancy and the B-factor is `-kBT log(p/(1-p))`; it is not an absolute binding free energy. Defaults use TIP3P, the `4r` protein-water dielectric, `8 x 12 x 6` C2-reduced orientation quadrature and the calibrated TIP3P network chemical potential. TIP4P requires an explicitly calibrated `--water_network_mu`.

Run `./watvina --help_water` for field, dielectric and site controls.

The predicted `waters.pdb` can be used as a discrete-water displacement bias:

```bash
./watvina -c vina.in -r receptor.pdb -l ligand.sdf \
  -w waters.pdb --water_bias_weight 0.1 -o docked_with_water.sdf
```

This bias does not add explicit water bridges or re-solve the FHFT network for each pose.

## Pharmacophore Workflows

Generate a receptor-pocket pharmacophore in the configured box:

```bash
./watvina -c vina.in -r receptor.pdb --genph4 -o ph4.cif
```

Apply a CIF pharmacophore during docking:

```bash
./watvina -c vina.in -r receptor.pdb -l ligand.sdf \
  --template ph4.cif -o docked_ph4.sdf
```

Run `./watvina --help_ph4` for the pharmacophore-specific controls.

## Python

Linux CPython 3.12:

```bash
PYTHONPATH=/path/to/this/directory python3.12 -c \
  'import watvina_python as wv; print(wv.__version__)'
```

Windows PowerShell with x86-64 CPython 3.12:

```powershell
$env:PYTHONPATH = "C:\path\to\this\directory"
python -c "import watvina_python as wv; print(wv.__version__)"
```

See `watvina-python-guide.md` for complete examples covering file and in-memory inputs, docking, score-only evaluation, structured results, pose serialization, mini relaxation, FHFT fields and water-site validation.

## Common Failures

| Symptom | Check |
|---|---|
| Missing `GLIBC_*`/`GLIBCXX_*` | Use `watvina-glibc2.17` on older Linux |
| Python import fails | Match OS, x86-64 architecture and CPython 3.12 exactly |
| Weak or missing H-bonds | Confirm all hydrogens and protonation states are explicit |
| Parse or typing error | Confirm SDF bond orders/formal charges and PDB residue identity |
| Empty or implausible search | Check box center/size and ligand coordinates |
| Non-reproducible ranks | Set a non-zero `--seed` and keep CPU/search settings fixed |

The four help screens are references for exact option defaults; this guide describes when and why to use the main workflows:

```bash
./watvina --help
./watvina --help_relax
./watvina --help_water
./watvina --help_ph4
```
