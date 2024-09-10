import watvina_python as wvpy
from rdkit2pdbqt import *

wv = wvpy.WATVina()
wv.__init__(cpu=8, seed=0, verbosity=1)

#从pdb文件读取受体
receptor_mol=Chem.MolFromPDBFile("rec.pdb", removeHs=False)
receptor_lines=MolToPDBQTBlock(receptor_mol, False, False, True)
wv.set_receptor_from_string(receptor_lines)

# Step 2.x:  读取水，或者药效团设置
#wv.set_pharmacophore_from_file("pharm.txt")
#wv.set_water("water.pdb",implicitsol=False)

# Step 3: 设置打分权重和格点
wv.set_watvina_weights(weight_vdw=0.184,weight_hb=1.00,weight_elep=0.25)
wv.set_extra_constraints(weight_desol=-0.500, wclash_dist=0.500, weight_torsion=0.300)
wv.set_grid_dims(center_x=0.55,center_y=30,center_z=16.5,size_x=19,size_y=23,size_z=23,granularity=0.375)

#Step 4: 根据打分方程预先计算不同距离的打分, 这个在global_search，score, optimize, relax时候都需要
#不带参数则只储存配体原子类型的格点能量
#wv.set_precalculate_sf()
#watvina的格点能量计算需要点时间, 所以我们如果针对不同的配体最好预先计算所有原子类型
wv.set_precalculate_sf(prec_full_atomtypes=True)


#计算所有原子类型的格点
#后面所有的配体都可以用这一个格点能量，不需要重复计算
wv.compute_watvina_maps()

#Step 5 : 从sdf文件中读取配体
ligand_mol=Chem.MolFromMolFile("i.sdf",removeHs=False)
ligand_lines=MolToPDBQTBlock(ligand_mol, True, False, True)
wv.set_ligand_from_string(ligand_lines)
wv.pose_atomids=[x-1 for x in wv.pose_atomids]
#wv.pose_atomids
wv.score()
wv.optimize()

wv.global_search(exhaustiveness=8,n_poses=5, min_rmsd=1.5,energy_range=3,population_size=8,
                 ga_searching=4,refinement=True, tramplitude=1.00)

print(wv.poses_score)
