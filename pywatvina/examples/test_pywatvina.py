#!/usr/bin/env python
# coding: utf-8

# # pyWATVina (Python API for watvina) 教程
#-------------------------------------
# ## pyWATVina的安装流程：
# =======================
# **pyWATVina的安装较为复杂，请按照以下流程安装** 
# -------------------
# 1. **解压 pywatvina.zip**
# `unzip pywatvina.zip`
# 
# 2. **解压 pywatvina的静态库**
# `cd watvina`
# `ar -x libwatvina.a` 
# 该步骤会生成若干的 `.o` 格式的文件
#
# 3. **重新编译成动态库，以便python调用**
# **注意要在自己python环境下的boost版本 **
# `g++ -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro -g -fwrapv -O2 -Wl,-Bsymbolic-functions -Wl,-z,relro -g -fwrapv -O2 -g -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 *.o -o _watvina_wrapper.so -lboost_thread -lboost_serialization -lboost_filesystem -lboost_program_options`
# 该步骤会生成 `_watvina_wrapper.so` 动态库
#
# 4. **删除不必要的文件** 
# `$rm *.a *.o`
# 
# 5. **将本文件夹`watvina`拷贝到python环境** 
# 一般是 `${你的python环境下lib/python3.x/dist-pakcages/}`

# Step 1: 加载watvina模块和rdkit2pdbqt.py的函数
from watvina.watvina_wrapper import WATVina
from watvina.rdkit2pdbqt import *

# Step 2: 建立一个对接任务，采用8个cpu核
wv=WATVina(8)

#从pdb文件读取受体
receptor_mol=Chem.MolFromPDBFile("rec.pdb", removeHs=False)
receptor_lines=MolToPDBQTBlock(receptor_mol, False, False, True)
wv.set_receptor_from_string(receptor_lines)

# Step 2.x:  读取水，或者药效团设置
#wv.set_pharmacophore_from_file("pharm.txt")
#wv.set_water("water.pdb",implicitsol=False)

# Step 3: 设置打分权重和格点
wv.set_watvina_weights(weight_vdw=0.193,weight_hb=0.6,weight_elep=0.15)
wv.set_extra_constraints(weight_desol=-0.500, wclash_dist=0.500, weight_torsion=0.300)
wv.set_grid_dims(center_x=0.55,center_y=30,center_z=16.5,size_x=19,size_y=23,size_z=23,granularity=0.35)

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


#Step 5.x1 只给初始构象打分
#wv.score()

#Step 5.x2 只对初始构象进行松弛(relax)
#wv.relax_structure(relax_steps=10000, tramplitude=0.001, rotamplitude=1.00)

#Step 5.x3 只对初始构象进行局部优化
#wv.optimize(max_steps=100)

#全局搜索构象
wv.global_search(exhaustiveness=12,n_poses=5,min_rmsd=1.5,energy_range=3,population_size=8,
                 ga_searching=4,refinement=True)

print(f'Watvina genetated {len(wv.poses_coords)} conformation.')
print(f'The score of the 1st conformation: {wv.poses_score[0]}')

#将对接后的坐标，能量等写入新的构象，并输出文件
#在后面的循环当中，如果只设置conf_w的属性，并不能将属性输出到i_out.sdf文件中，
#如果只设置ligand_mol的属性，最后一个构象缺失属性
from rdkit import Chem
from rdkit.Geometry import Point3D

#print(f'default conformers num: {ligand_mol.GetNumConformers()}')
conf = ligand_mol.GetConformer(0)

writer = Chem.SDWriter('i_out.sdf')
for pose_id in range(len(wv.poses_coords)):
    for atomid,coord in zip(wv.pose_atomids,wv.poses_coords[pose_id]):
        conf.SetAtomPosition(atomid,(coord[0],coord[1],coord[2]))
    #第0个构象是原来的，因此pose_id是0的时候，conf_id是1
    ligand_mol.AddConformer(conf,assignId=True)
    #处理第pose_id + 1的构象
    
    conf_id   = pose_id + 1
    conf_w = ligand_mol.GetConformer(conf_id)
    conf_score='%.2f' % wv.poses_score[pose_id]
    conf_vdw='%.2f' % wv.poses_vdw[pose_id]
    conf_hbond='%.2f' % wv.poses_hbond[pose_id]
    conf_electrop='%.2f' % wv.poses_electrop[pose_id]

    ligand_mol.SetProp('ConfID', f'{conf_id}')
    ligand_mol.SetProp("Score",conf_score)
    ligand_mol.SetProp("VDW",conf_vdw)
    ligand_mol.SetProp("Hbond",conf_hbond)
    ligand_mol.SetProp("Electrop",conf_electrop)
    
    conf_w.SetProp('ConfID', f'{conf_id}')
    conf_w.SetProp("Score",conf_score)
    conf_w.SetProp("VDW",conf_vdw)
    conf_w.SetProp("Hbond",conf_hbond)
    conf_w.SetProp("Electrop",conf_electrop)    
    writer.write(ligand_mol,confId=conf_id)
    #print(f'Processing pose_id {pose_id} and conf_id {conf_id}, and conf real id is {conf_w.GetId()}')

ligand_mol.GetNumConformers()

