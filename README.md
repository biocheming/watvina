# Water Model supported protein-ligand docking with Autodock Vina engine. 

![Watvina](https://github.com/biocheming/watvina/blob/main/watvina_logo.png)

For drug design purpose, explicit or implicit waters, pharmacophore or position constrained docking, external torsion parameter (in amber/gaff/charmm force field like parameters) are supported in watvina. 

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
ATOM      1 CH   HVY P   1      -9.099 208.885 112.556  1.00  0.20           C
ATOM      2 CA   ARO P   2      -9.099 208.885 112.556  0.70  0.36           C
ATOM      3 CH   HVY P   3     -10.107 208.439 111.792  1.00  0.20           C
ATOM      4 CH   HVY P   4      -9.242 209.208 114.013  1.00  0.20           C
ATOM      5 CA   ARO P   5      -9.242 209.208 114.013  0.70  0.25           C
ATOM      6 CH   HVY P   6      -8.067 209.594 114.784  1.00  0.20           C
ATOM      7 CH   HVY P   7      -6.765 209.768 114.098  1.00  0.20           C
ATOM      8 CA   ARO P   8      -6.765 209.768 114.098  0.70  0.27           C
ATOM      9 CH   HVY P   9      -5.422 210.380 115.054  1.00  0.20           C
ATOM     10 CH   HVY P  10      -7.783 209.022 111.886  1.00  0.20           C
```
where ```resname``` for pharmacophore type; ```occupancy```for cutoff distance and ```b-factor``` for contribution weight.

If only keep the ```HVY```, the heavy atoms, the model is quite similar to a molecular shape,
while other pharmacophores for the colors in shape.

other pharmacophores: 

```DON```: hbond donor;

```ARO```: for aromatic carbons;

```PCG```: for positively charged nitroge or guanidine carbon;

```NCG```: for negatively charged center;

```SGM```: Cl, Br, I, S;

```HYP```: hydrophobic atoms(hydrophobic carbon, Cl, Br, I)

```HVY```: any heavy atoms(not hydrogen)

### 2.4 position constrained docking.

position constrained docking is useful for FEP, enzymatic pre-active pose prediction...

```
watvina --config vina.conf --tramplitude 0
```

```--tramplitude 0``` will frozen the ```ROOT```  of the ligand. 

### 2.6 rigid docking
```
watvina --config vina.conf --toramplitude 0
```
```--toramplitude 0``` will frozen the torsions of the ligand

### 2.7 The pdbqt files for receptors and ligands are prepared from their pdb file by mgltools, or from rdkit2pdbqt.py(using opencadd,)

```
rdkit2pdbqt -l lig.sdf
rdkit2pdbqt -r rec.pdb
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
Input:
  -r [ --receptor ] arg                rigid part of the receptor (PDBQT)
  --flex arg                           flexible part of the receptor (pdbqt)
  --template arg                       template p_ph4 (pdb)
  -w [ --water ] arg                   water file (O coordicates file with 
                                       resname HOH, energy in the beta column)
  --pharma arg                         pharmacophore restraints
  -l [ --ligand ] arg                  ligand.pdbqt(ligand file, PDBQT)
  --ligands_dir arg                    directory for ligands 

Search space (required):
  --center_x arg                       X coordinate of the center
  --center_y arg                       Y coordinate of the center
  --center_z arg                       Z coordinate of the center
  --size_x arg                         size in the X dimension (Angstroms)
  --size_y arg                         size in the Y dimension (Angstroms)
  --size_z arg                         size in the Z dimension (Angstroms)

Output (optional):
  -o [ --out ] arg                     output models (PDBQT), the default is 
                                       chosen based on the ligand file name
  --out_dir arg (=OUTPUT)              output directory for batch mode
  --log arg                            optionally, write log file
  --genph4                             generate pseudo pharmacophore model in 
                                       pdb format

Advanced options (see the manual):
  --score_only                         score only - search space can be omitted
  --local_only                         do local search only
  --implicitsol                        implicit solvation model result in a 
                                       implicitsol.pdb
  --grid_space arg (=0.375)            grid space
  --weight_vdw arg (=0.184)            vdw weight
  --weight_hbond arg (=1)              Hydrogen bond weight
  --weight_electrop arg (=0.25)        polar repulsion or hydrophobic 
                                       attraction
  --weight_desol arg (=-0.5)           desolvation weight[depends on water 
                                       model used]
  --wclash_dist arg (=0.5)             clash distance with water[depends on 
                                       water model used]
  --weight_torsion arg (=0.300)       external torsion weight[depends on 
                                       forcefield or unit in kj/mol or 
                                       kcal/mol]
  --relax_only                         do relax only without BFGS refinement 
                                       for local searching
  --local_steps arg (=2000)            local relax steps
  --tramplitude arg (=1)               amplitude for translation/rotation in 
                                       local relax
  --toramplitude arg (=1)              amplitude for torsion in local relax

Misc (optional):
  --cpu arg                            the number of CPUs to use (the default 
                                       is to try to detect the number of CPUs 
                                       or, failing that, use 1)
  --seed arg                           explicit random seed
  --exhaustiveness arg (=8)            exhaustiveness of the global search 
                                       (roughly proportional to time): 1+
  --population arg (=8)                population size for genetic algorithm
  --ga_search arg (=4)                 amplitude for ga searching loop size: 0,
                                       1, 2...
  --num_modes arg (=10)                maximum number of binding modes to 
                                       generate
  --rmsd arg (=1.5)                    modes clustering cutoff
  --energy_range arg (=3)              maximum energy difference between the 
                                       best binding mode and the worst one 
                                       displayed (kcal/mol)

Configuration file (optional):
  -c [ --config ] arg                  the above options can be put here

Information (optional):
  --help                               display usage summary
  --help_advanced                      display usage summary with advanced 
                                       options
  --verbosity arg                      display IO information
  --version                            display program version
```
