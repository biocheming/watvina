# Water Model supported protein-ligand docking with Autodock Vina engine. 

For drug design purpose, explicit or implicit waters, pharmacophore or position constrained docking, external torsion parameter (in amber/gaff/charmm force field like parameters) are supported in watvina. 

## 1. IMPORTANT NOTES
### 1.1 Receptor and ligand(s) should be carefully prepared. KEEP ALL THE NON-POLAR AND POLAR HYDROGENS. 
### 1.2 Minimize/Optimize the ligand carefully. Watvina calculates some torsion penalty.
### 1.3 Hydroxy hydrogens which do not formed intra hydrogen bond in the receptor are suggested to be flexible. Another method is to change the OA and HD atom types to OW and HW respectively.
### 1.4 Powershell is not recommended for ```rdkit2pdbqt.py```


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

### 2.3 USAGE for pharmacophore constrained docking

```
watvina --config vina.conf --pharma a_txt_file.txt
```
the format of ```a_txt_file.txt```is:

```name x y z cut_off_distance weight```

for example:

```ACC 0.0 0.0 0.0 0.3 0.5```

```ACC```: the hbond acceptor

```0.0 0.0 0.0```: the coordinate of hbond acceptor 

```0.3```: the cut off distance is 0.3

```0.5```: the award weight is 0.5, finally the score is (1+0.5)*watvina_score

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
watvina --config vina.conf --pharma a_txt_file.txt
```

similar to pharmacophore constrained docking, the position constrained pharmacophore name is:

```BFAX```

The ```X``` is the atom serial number in the ligand pdbqt file.

for example, ```BFA4```, the atom with serial number ```4``` will be fixed to a position.

### 2.5 The pdbqt files for receptors and ligands are prepared from their pdb file by mgltools, or from rdkit2pdbqt.py(using opencadd,)

```
rdkit2pdbqt -l lig.sdf
rdkit2pdbqt -r rec.pdb
```

### 2.6 External torsion parameter

External torsion parameters in the header of ligand file

the format is
```
REMARK TORSION INDEX  i  j  k  l  V/2   theta_0  n
```
for examle:

```
REMARK TORSION INDEX  18  17  16  21  0.16   0  3
```

### 2.7 Help information
```
Input:
  --receptor arg                        rigid part of the receptor (PDBQT)
  --flex arg                            flexible part of the receptor (pdbqt)
  --water arg                           water file (O coordicates file with 
                                        resname HOH, energy in the beta column)
  --pharma arg                          pharmacophore restraints
  --ligand arg                          ligand.pdbqt(ligand file, PDBQT)

Search space (required):
  --center_x arg                        X coordinate of the center
  --center_y arg                        Y coordinate of the center
  --center_z arg                        Z coordinate of the center
  --size_x arg                          size in the X dimension (Angstroms)
  --size_y arg                          size in the Y dimension (Angstroms)
  --size_z arg                          size in the Z dimension (Angstroms)

Output (optional):
  --out arg                             output models (PDBQT), the default is 
                                        chosen based on the ligand file name
  --log arg                             optionally, write log file

Advanced options (see the manual):
  --score_only                          score only - search space can be 
                                        omitted
  --local_only                          do local search only
  --implicitsol                         implicit solvation model result in a 
                                        implicitsol.pdb
  --grid_space arg (=0.375)             grid space
  --randomize_only                      randomize input, attempting to avoid 
                                        clashes
  --weight_vdw arg (=0.193000004)       vdw weight
  --weight_hbond arg (=0.600000024)     Hydrogen bond weight
  --weight_electrop arg (=0.140000001)  polar repulsion or hydrophobic 
                                        attraction
  --weight_desol arg (=-1)              desolvation weight[depends on water 
                                        model used]
  --wclash_dist arg (=1.39999998)       clash distance with water[depends on 
                                        water model used]
  --weight_torsion arg (=0.300000012)   external torsion weight[depends on 
                                        forcefield or unit in kj/mol or 
                                        kcal/mol]

Misc (optional):
  --cpu arg                             the number of CPUs to use (the default 
                                        is to try to detect the number of CPUs 
                                        or, failing that, use 1)
  --seed arg                            explicit random seed
  --exhaustiveness arg (=1)             exhaustiveness of the global search 
                                        (roughly proportional to time): 1+
  --population arg (=10)                 population size for genetic algorithm
  --ga_search arg (=1)                  amplitude for ga searching loop size:1,2,3...
  --num_modes arg (=20)                 maximum number of binding modes to 
                                        generate
  --rmsd arg (=1.5)                     modes clustering cutoff
  --energy_range arg (=3)               maximum energy difference between the 
                                        best binding mode and the worst one 
                                        displayed (kcal/mol)

Configuration file (optional):
  --config arg                          the above options can be put here

Information (optional):
  --help                                display usage summary
  --help_advanced                       display usage summary with advanced 
                                        options
  --version                             display program version

```
