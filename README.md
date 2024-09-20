# Water Model supported protein-ligand docking with Autodock Vina engine. 

![Watvina](https://github.com/biocheming/watvina/blob/main/watvina_logo.png)


Watvina facilitates drug design with support for explicit or implicit waters, pharmacophore or position-constrained docking, and external torsion parameters (akin to amber/gaff/charmm force fields).
 
## 1. IMPORTANT NOTES
### 1.1 Receptor and ligand(s) should be carefully prepared. KEEP ALL THE NON-POLAR AND POLAR HYDROGENS. 
### 1.2 Minimize/Optimize the ligand carefully. Watvina calculates some torsion penalty.
### 1.3 Hydroxy hydrogens which do not formed intra hydrogen bond in the receptor are suggested to be flexible. Another method is to change the OA and HD atom types to OW and HW respectively.
### 1.4 The ```pywatvina``` surpports ```pdb``` and ```pdbqt``` formats for receptor, ```sdf``` and ```pdbqt``` for ligand


## 2 USAGE
### 2.1 USAGE for explicit water model

Opendx format was not supported now, but later on it will be online. In addition, the desolvation weight has to be adjusted mannually due to the energy calculated from different methods.

```
watvina --config vina.conf --water water.pdb
```
In the ```water.pdb``` file with energy value (calculated from GIST, WATSite, Watermap, lesite etc) in the beta factor column. Keep the Oxygen atoms only, and with a resname ```HOH```. 

### 2.2 USAGE for 'implicit water' model

Desolvation was calculated from water probe generated energy map

```
watvina --config vina.conf --implicitsol
```

### 2.3 USAGE for template based docking

```
watvina --config vina.conf --template a_pseudo_pharmacophore_pdb_file.pdb
```
the format of ```a_pseudo_pharmacophore_pdb_file.pdb```is in ```pdb```format

```
ATOM      1 CH   HVY P   2      -0.605  29.778  21.561  2.00  0.20           C
ATOM      2 CH   HVY P  12       1.270  28.887  19.820  2.00  0.20           C
ATOM      3 CH   HVY P  13       2.006  27.868  19.172  2.00  0.20           C
ATOM      4 CA   ARO P  13       2.006  27.868  19.172  0.70  0.27           C
ATOM      5 CH   HVY P  17       2.939  28.100  18.157  2.00  0.20           C
ATOM      6 ND   DON P  17       2.939  28.100  18.157  1.00  0.41           N
ATOM      7 CH   HVY P  14       1.727  26.553  19.510  2.00  0.20           C
ATOM      8 CH   HVY P   3      -0.477  29.435  22.941  2.00  0.20           C
ATOM      9 OA   ACC P   3      -0.477  29.435  22.941  1.00  0.48           O
ATOM     10 CH   HVY P   7      -3.861  31.122  17.827  2.00  0.20           C
ATOM     11 CH   HVY P  18       2.859  29.108  17.251  2.00  0.20           C
ATOM     12 CH   HVY P  20       2.164  30.901  16.325  2.00  0.20           C
ATOM     13 CH   HVY P  23       1.204  31.967  16.050  2.00  0.20           C
ATOM     14 CA   ARO P  23       1.204  31.967  16.050  0.70  0.26           C
ATOM     15 CH   HVY P  21       3.350  30.513  15.703  2.00  0.20           C
ATOM     16 CA   ARO P  21       3.350  30.513  15.703  0.70  0.43           C
ATOM     17 CH   HVY P  22       3.811  29.380  16.306  2.00  0.20           C
ATOM     18 OA   ACC P  22       3.811  29.380  16.306  1.00  0.70           O
ATOM     19 CH   HVY P  24       1.152  32.577  14.797  2.00  0.20           C
ATOM     20 CA   ARO P  24       1.152  32.577  14.797  0.70  0.30           C

```
where ```resname``` for pharmacophore type; ```occupancy```for cutoff distance and ```b-factor``` for contribution weight.

If only keep the ```HVY```, the heavy atoms, the model is quite similar to a molecular shape,
while other pharmacophores for the colors in shape.
watvina can generate a initial template ```ph4.pdb``` file based on the protein and ligand interaction.
```
watvina --config vina.conf --score_only --genph4
```

other pharmacophores: 

```DON```: hbond donor;

```ARO```: for aromatic carbons;

```PCG```: for positively charged nitroge or guanidine carbon;

```NCG```: for negatively charged center;

```SGM```: Cl, Br, I, S;

```HYP```: hydrophobic atoms(hydrophobic carbon, Cl, Br, I)

```HVY```: any heavy atoms(not hydrogen)

### 2.4 position[scaffold, core] constrained docking.

position constrained docking is useful for FEP, enzymatic pre-active pose prediction...

```
watvina --config vina.conf --tramplitude 0
```

```--tramplitude 0``` will freeze the ```ROOT```  of the ligand. 

### 2.6 rigid docking
```
watvina --config vina.conf --toramplitude 0
```
```--toramplitude 0``` will frozen the torsions of the ligand

### 2.7 The pdbqt files for receptors and ligands are prepared from their pdb file by mgltools, or from rdkit2pdbqt.py

```
rdkit2pdbqt.py -l lig.sdf
rdkit2pdbqt.py -r rec.pdb
```
or setup the scaffold as ROOT of the pdbqt
```
rdkit2pdbqt.py -l lig.sdf core.sdf
```
### 2.8 External torsion parameter

External torsion parameters in the header of ligand file

the format is
```
REMARK TORSION INDEX  i  j  k  l  V/2   theta_0  n
```
for examle:
```
REMARK TORSION INDEX  18  17  16  21  0.16   0  3
```
### 2.9 Help information
```
Input options:
  -r, --receptor arg         Rigid part of the receptor (PDBQT)
      --flex arg             Flexible part of the receptor (pdbqt)
      --template arg         Template ph4 (pdb)
  -w, --water arg            Water file (O coordinates file with resname HOH, 
                             energy in the beta column)
      --pharma arg           Pharmacophore[ph4] constrained file
      --ligands_dir arg      Directory for ligands
      --multiligs_pdbqt arg  PDBQT file containing multi ligands
  -l, --ligand arg           Ligand.pdbqt(ligand file, PDBQT/SDF)

 Output options:
  -o, --out arg             Output models (PDBQT/SDF)
      --out_dir arg         Output directory for batch mode (default: OUTDIR)
      --score_cutoff arg    The cutoff score for output (default: 10000)
      --ph4_cutoff arg      The cutoff of ph4 for output (default: -1.0)
      --timeout_cutoff arg  The timeout for each thread (default: 600)
      --log arg             Optionally, write log file

 Search Space options:
      --center_x arg  X coordinate of the center
      --center_y arg  Y coordinate of the center
      --center_z arg  Z coordinate of the center
      --size_x arg    Size in the X dimension (Angstroms) (default: 25)
      --size_y arg    Size in the Y dimension (Angstroms) (default: 25)
      --size_z arg    Size in the Z dimension (Angstroms) (default: 25)

 Misc options:
      --cpu arg             The number of CPUs to use
      --seed arg            Explicit random seed
      --exhaustiveness arg  Exhaustiveness of the global search (default: 8)
      --population arg      Population size for genetic algorithm (default: 8)
      --ga_search arg       Amplitude for ga searching loop size (default: 4)
      --num_modes arg       Maximum number of binding modes to generate 
                            (default: 10)
      --rmsd arg            Modes clustering cutoff (default: 1.5)
      --energy_range arg    Maximum energy difference between the best binding 
                            mode and the worst one displayed (kcal/mol) 
                            (default: 3.0)

 Config File options:
  -c, --config arg  The above options can be put here

 Advanced options:
      --score_only           Score only - search space can be omitted
      --genph4               Generate pseudo pharmacophore model in pdb format
      --local_only           Do local search only
      --implicitsol          Implicit solvation model result in an 
                             wgrid.pdb
      --grid_space arg       Grid space (default: 0.375)
      --weight_vdw arg       Vdw weight (default: 0.184)
      --weight_hbond arg     Hydrogen bond weight (default: 1.000)
      --weight_electrop arg  Polar repulsion or hydrophobic attraction 
                             (default: 0.250)
      --weight_desol arg     Desolvation weight [depends on water model used] 
                             (default: -0.500)
      --wclash_dist arg      Clash distance with water [depends on water model 
                             used] (default: 0.500)
      --weight_torsion arg   External torsion weight [depends on forcefield or 
                             unit ] (default: 0.300)
      --relax_only           Do relax only without BFGS refinement for local 
                             searching
      --local_steps arg      Local relax steps (default: 2000)
      --tramplitude arg      Amplitude for translation/rotation in local relax 
                             (default: 1.000)
      --toramplitude arg     Amplitude for torsion in local relax (default: 
                             1.000)

 Information options:
  -h, --help           Display usage summary
      --help_advanced  Display usage summary with advanced options
      --verbosity arg  Display IO information (default: 2)
  -v, --version        Display program version
```
