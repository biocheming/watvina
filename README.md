# watvina
implicit or explicit water model in protein-ligand docking with vina engine. supporting pharmacophore /position constrained docking
IMPORTANT NOTES:
1. receptor and ligand should be carefully prepared. Keep all the non-polar and polar hydrogens.
2. minimize/optimize the ligand carefully. watvina calculates some torsion penalty.

USAGE: 
1. explicit water model
```watvina --config vina.conf --water water.pdb```
In the water.pdb file, the energy value in beta factor column, which can be from gist, watsite, 3d-rism, lesite etc.
2. implicit water model
```watvina --config vina.conf --implicitsol```
3. pharmacophore constrained docking
```watvina --config vina.conf --pharma a_txt_file.txt```
format in the a_txt)fle.txt
```name x y z cut_off_distance weight```
for example
```ACC 0.0 0.0 0.0 0.3 0.5```
ACC: hbond acceptor
0.0 0.0 0.0: the coordinate of the center of ACC
0.3: within a cut off distance 0.3 A of the ACC center
0.5: (1+0.5)* watvina_score will be calculated

other pharmacophores,  DON for hbond donor, ARO for aromatic ring. 

4. position constrained docking
```watvina --config vina.conf --pharma a_txt_file.txt```
like the pharmacophore constrained docking, in the ```a_txt_file.txt```file,
```name x y z cut_off_distance weight```
for example
```BFA 0.0 0.0 0.0 0.3 0.5```
the atom with beta factor in the pdbqt file should be set to 100.0 for docking.
If this atom falls within the 0.3 A of BFA center, new score will be (1+0.5) * watvina_score.

