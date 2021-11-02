# Water Model supported protein-ligand docking with Autodock Vina engine. 

For drug design purpose, explicit or implicit waters, pharmacophore or position constrained docking are supported in watvina. 

## 1. IMPORTANT NOTES
### 1.1 Receptor and ligand(s) should be carefully prepared. KEEP ALL THE NON-POLAR AND POLAR HYDROGENS. 
### 1.2 Minimize/Optimize the ligand carefully. Watvina calculates some torsion penalty.
### 1.3 Hydroxy hydrogens are suggested to be flexible in docking.


## 2 USAGE
### 2.1 USAGE for explicit water model

```
watvina --config vina.conf --water water.pdb
```
In the ```water.pdb``` file with energy value (calculated from GIST, WATSite, Watermap, lesite etc) in the beta factor column. Keep the Oxygen atoms only, and with a resname ```HOH```. 

### 2.2 USAGE for implicit water model

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

```ARO```: for aromatic carbons.

### 2.4 position constrained docking.

position constrained docking is useful for FEP, enzymatic pre-active pose prediction...

```
watvina --config vina.conf --pharma a_txt_file.txt
```

similar to pharmacophore constrained docking, the constrained pharmacophore name is:

```BFA```

the atom in the ligand pdbqt file with beta factor value 100.0 is required. 


### 2.4 The pdbqt files for receptors and ligands are prepared from their pdb file by mgltools
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
  --randomize_only                      randomize input, attempting to avoid 
                                        clashes
  --weight_gauss arg (=-0.360000014)    gauss weight
  --weight_repulsion arg (=0.699999988) repulsion weight
  --weight_hydrogen arg (=0.600000024)  Hydrogen bond weight
  --weight_electrop arg (=0.400000006)  Electro polar weight
  --weight_desol arg (=0.0700000003)    desolvation weight
  --wclash_dist arg (=1.39999998)       clash distance with water

Misc (optional):
  --cpu arg                             the number of CPUs to use (the default 
                                        is to try to detect the number of CPUs 
                                        or, failing that, use 1)
  --seed arg                            explicit random seed
  --exhaustiveness arg (=8)             exhaustiveness of the global search 
                                        (roughly proportional to time): 1+
  --population arg (=5)                 population size for genetic algorithm
  --ga_search arg (=5)                  amplitude for ga searching loop size
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
