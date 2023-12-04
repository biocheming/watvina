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
ATOM     20 CH   HVY P  20      17.244  18.241  50.736  2.00  0.20           C
ATOM     20 ND   DON P  20      17.244  18.241  50.736  1.00  0.58           N
ATOM     18 CH   HVY P  18      18.452  17.752  51.268  2.00  0.20           C
ATOM     18 CA   ARO P  18      18.452  17.752  51.268  0.70  0.28           C
ATOM     21 CH   HVY P  21      16.957  19.505  50.313  2.00  0.20           C
ATOM     22 CH   HVY P  22      15.713  19.683  49.527  2.00  0.20           C
ATOM     22 CA   ARO P  22      15.713  19.683  49.527  0.70  0.28           C
ATOM     36 CH   HVY P  36      17.643  20.506  50.528  2.00  0.20           C
ATOM     36 OA   ACC P  36      17.643  20.506  50.528  1.00  0.56           O
ATOM     14 CH   HVY P  14      20.800  18.082  51.969  2.00  0.20           C
ATOM     14 CA   ARO P  14      20.800  18.082  51.969  0.70  0.28           C
ATOM     13 CH   HVY P  13      21.893  18.937  52.110  2.00  0.20           C
ATOM     13 ND   DON P  13      21.893  18.937  52.110  1.00  0.35           N
ATOM     15 CH   HVY P  15      20.895  16.683  52.364  2.00  0.20           C
ATOM     15 CA   ARO P  15      20.895  16.683  52.364  0.70  0.28           C
ATOM     37 CH   HVY P  37      22.181  16.074  52.953  2.00  0.20           C
ATOM     16 CH   HVY P  16      19.741  15.840  52.195  2.00  0.20           C
ATOM     16 CA   ARO P  16      19.741  15.840  52.195  0.70  0.36           C
ATOM     17 CH   HVY P  17      18.520  16.379  51.645  2.00  0.20           C
ATOM     17 CA   ARO P  17      18.520  16.379  51.645  0.70  0.28           C
ATOM     19 CH   HVY P  19      19.616  18.588  51.441  2.00  0.20           C

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

### 2.4 position constrained docking.

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
