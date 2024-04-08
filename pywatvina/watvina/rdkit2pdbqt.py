#!/usr/bin/env python3
'''
Original code from Open Drug Discovery Toolkit (ODDT), J Cheminform 7, 26 (2015). 
Modified by Ximing XU, Rongfeng ZOU, Hongrui LIN
xuximing@ouc.edu.cn
version 2023-12-07 tested with RDkit 2024.03.1pre
'''
from __future__ import absolute_import, print_function
from math import isnan, isinf
from itertools import combinations

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondType
import sys
from rdkit.Chem import rdFMCS
from rdkit.Chem.rdMolAlign import AlignMol
from rdkit.Chem import ChemicalForceFields

#### By Lin Hongrui
def get_core_alignment_for_template_docking(reference_mol, query_mol, core_atom_mapping_dict):
    # the initial position of query_mol is random, so align to the reference_mol firstly
    _ = AlignMol(query_mol, reference_mol, atomMap=list(core_atom_mapping_dict.items()))

    # assign template positions from reference mol to query mol
    core_fixed_query_conformer = Chem.Conformer(query_mol.GetNumAtoms())
    reference_conformer = reference_mol.GetConformer()
    query_conformer = query_mol.GetConformer()

    for query_atom_idx in range(query_mol.GetNumAtoms()):
        if query_atom_idx in core_atom_mapping_dict:
            reference_atom_idx = core_atom_mapping_dict[query_atom_idx]
            atom_position = reference_conformer.GetAtomPosition(reference_atom_idx)
            core_fixed_query_conformer.SetAtomPosition(query_atom_idx, atom_position)
        else:
            atom_position = query_conformer.GetAtomPosition(query_atom_idx)
            core_fixed_query_conformer.SetAtomPosition(query_atom_idx, atom_position)

    query_mol.RemoveAllConformers()
    query_mol.AddConformer(core_fixed_query_conformer)

    # optimize conformer using chemical forcefield
    ff_property = ChemicalForceFields.MMFFGetMoleculeProperties(query_mol, 'MMFF94s')
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(query_mol, ff_property, confId=0)

    for query_atom_idx in core_atom_mapping_dict.keys():
        reference_atom_idx = core_atom_mapping_dict[query_atom_idx]
        core_atom_position = reference_conformer.GetAtomPosition(reference_atom_idx)
        virtual_site_atom_idx = ff.AddExtraPoint(core_atom_position.x, core_atom_position.y, core_atom_position.z, fixed=True) - 1
        ff.AddDistanceConstraint(virtual_site_atom_idx, query_atom_idx, 0, 0, 20.0) #100 is too high, By Rongfeng

    ff.Initialize()

    max_minimize_iteration = 5
    for _ in range(max_minimize_iteration):
        minimize_seed = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
        if minimize_seed == 0:
            break

    query_mol.SetProp('aligned_conformer_energy', str(ff.CalcEnergy()))

    core_atom_idx_list = list(core_atom_mapping_dict.keys())

    return core_atom_idx_list
####
##From ODDT, modified by Ximing
def PDBQTAtomLines(mol, donors, acceptors):
    atom_lines = [line.replace('HETATM', 'ATOM  ')
                  for line in Chem.MolToPDBBlock(mol).split('\n')
                  if line.startswith('HETATM') or line.startswith('ATOM')]

    pdbqt_lines = []
    for idx, atom in enumerate(mol.GetAtoms()):
        pdbqt_line = atom_lines[idx][:56]

        pdbqt_line += '0.00  0.00    '  # append empty vdW and ele
        # Get charge
        charge = 0.
        fields = ['_MMFF94Charge', '_GasteigerCharge', '_TriposPartialCharge']
        for f in fields:
            if atom.HasProp(f):
                charge = atom.GetDoubleProp(f)
                break
        # FIXME: this should not happen, blame RDKit
        if isnan(charge) or isinf(charge):
            charge = 0.
        pdbqt_line += ('%.3f' % charge).rjust(6)

        # Get atom type
        pdbqt_line += ' '
        atomicnum = atom.GetAtomicNum()
        atomhybridization = atom.GetHybridization()
        atombondsnum = atom.GetDegree()
        if atomicnum == 6 and atom.GetIsAromatic():
            pdbqt_line += 'A '
        elif atomicnum == 7 and idx in acceptors:
            pdbqt_line += 'NA'
        elif atomicnum == 8 and idx in acceptors:
            pdbqt_line += 'OA'
        elif atomicnum == 1 and atom.GetNeighbors()[0].GetIdx() in donors:
            pdbqt_line += 'HD'
        elif atomicnum == 1 and atom.GetNeighbors()[0].GetIdx() not in donors:
            pdbqt_line += 'H '
        elif atomicnum == 16 and ( (atomhybridization == Chem.HybridizationType.SP3 and atombondsnum != 4) or atomhybridization == Chem.HybridizationType.SP2 ):
            pdbqt_line += 'SA'
        else:
            if len(atom.GetSymbol()) >1:
                pdbqt_line += atom.GetSymbol()    
            else:
                pdbqt_line += (atom.GetSymbol() + ' ') 			
        pdbqt_lines.append(pdbqt_line)
    return pdbqt_lines

def MolCoreToPDBQTBlock(mol, core_atom_idx_list, flexible=True, addHs=False, computeCharges=False):
    # make a copy of molecule
    mol = Chem.Mol(mol)

    # if flexible molecule contains multiple fragments write them separately
    if flexible and len(Chem.GetMolFrags(mol)) > 1:
        return ''.join(MolToPDBQTBlock(frag, flexible=flexible, addHs=addHs, computeCharges=computeCharges)
                       for frag in Chem.GetMolFrags(mol, asMols=True))

    # Identify donors and acceptors for atom typing
    # Acceptors
    patt = Chem.MolFromSmarts('[$([O;H1;v2]),'
                              '$([O;H0;v2;!$(O=N-*),'
                              '$([O;-;!$(*-N=O)]),'
                              '$([o;+0])]),'
                              '$([n;+0;!X3;!$([n;H1](cc)cc),'
                              '$([$([N;H0]#[C&v4])]),'
                              '$([N&v3;H0;$(Nc)])]),'
                              '$([F;$(F-[#6]);!$(FC[F,Cl,Br,I])])]')
    acceptors = list(map(lambda x: x[0], mol.GetSubstructMatches(patt, maxMatches=mol.GetNumAtoms())))
    # Donors
    patt = Chem.MolFromSmarts('[$([N&!H0&v3,N&!H0&+1&v4,n&H1&+0,$([$([Nv3](-C)(-C)-C)]),'
                              '$([$(n[n;H1]),'
                              '$(nc[n;H1])])]),'
                              # Guanidine can be tautormeic - e.g. Arginine
                              '$([NX3,NX2]([!O,!S])!@C(!@[NX3,NX2]([!O,!S]))!@[NX3,NX2]([!O,!S])),'
                              '$([O,S;H1;+0])]')
    donors = list(map(lambda x: x[0], mol.GetSubstructMatches(patt, maxMatches=mol.GetNumAtoms())))
    if addHs:
        mol = Chem.AddHs(mol, addCoords=True, onlyOnAtoms=donors, )
    if addHs or computeCharges:
        AllChem.ComputeGasteigerCharges(mol)

    atom_lines = PDBQTAtomLines(mol, donors, acceptors)
    assert len(atom_lines) == mol.GetNumAtoms()

    pdbqt_lines = []

    pdbqt_lines.append('REMARK  Name = ' + (mol.GetProp('_Name') if mol.HasProp('_Name') else ''))
    if flexible:
        # Find rotatable bonds
        #rot_bond = Chem.MolFromSmarts('[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]') # From Chemaxon
        rot_bond  = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]') #single and not ring, really not in ring?
        bond_atoms = list(mol.GetSubstructMatches(rot_bond))
        
        exclude_list = ['[NX3]-[CX3]=[O,N]', 'C(=[O;!R])N', '[!#1]-[C;H3,Cl3,F3,Br3]']
        excluded_atoms = []
        for excluded_smarts in exclude_list:
            excluded_bond = Chem.MolFromSmarts(excluded_smarts)
            excluded_atoms += [(x[0],x[1]) for x in list(mol.GetSubstructMatches(excluded_bond))]

        for excld_atom in excluded_atoms:
            excld_atom_reverse=(excld_atom[1],excld_atom[0])
            if excld_atom in bond_atoms:
                bond_atoms.remove(excld_atom)
            elif excld_atom_reverse in bond_atoms:
                bond_atoms.remove(excld_atom_reverse)

        tertiary_amide_bonds = Chem.MolFromSmarts('[NX3]([!#1])([!#1])-[CX3]=[O,N]')
        tertiary_amide_bond_atoms=[(x[0],x[3]) for x in list(mol.GetSubstructMatches(tertiary_amide_bonds))]
        if len(tertiary_amide_bond_atoms) > 0:
            for tertiary_amide_bond_atom in tertiary_amide_bond_atoms:
                bond_atoms.append(tertiary_amide_bond_atom)
###By Xu Ximing
        to_delete_core_atoms = []
        for i in range(len(bond_atoms)):
            reverse_bond_atoms_i=(bond_atoms[i][1],bond_atoms[i][0])
            for j in core_atom_idx_list:
                atom_j = mol.GetAtomWithIdx(j)
                for neighbor_atom in atom_j.GetNeighbors():
                    neighbor_atom_index = neighbor_atom.GetIdx()
                    if neighbor_atom_index in core_atom_idx_list:
                        bond_atoms_in_core = (j, neighbor_atom_index)
                        if bond_atoms_in_core == bond_atoms[i] or bond_atoms_in_core == reverse_bond_atoms_i:
                            to_delete_core_atoms.append(bond_atoms[i])

        to_delete_core_atoms = list(dict.fromkeys(to_delete_core_atoms))
        for i in to_delete_core_atoms:
            bond_atoms.remove(i)


        atom_lines = PDBQTAtomLines(mol, donors, acceptors) # update coordinate
####
        # Fragment molecule on bonds to ge rigid fragments
        bond_ids = [mol.GetBondBetweenAtoms(a1, a2).GetIdx()
                    for a1, a2 in bond_atoms]
        
        if bond_ids:
            for i, b_index in enumerate(bond_ids):
                tmp_frags= Chem.FragmentOnBonds(mol, [b_index], addDummies=False)
                tmp_frags_list=list(Chem.GetMolFrags(tmp_frags))
                #tmp_bigger=0
                if len(tmp_frags_list) == 1:
                    del bond_ids[i]
                    del bond_atoms[i]

            mol_rigid_frags = Chem.FragmentOnBonds(mol, bond_ids, addDummies=False)
        else:
            mol_rigid_frags = mol

        num_torsions = len(bond_atoms)
        # Active torsions header
        pdbqt_lines.append('REMARK  %i active torsions:' % num_torsions)
        pdbqt_lines.append('REMARK  status: (\'A\' for Active; \'I\' for Inactive)')
        for i, (a1, a2) in enumerate(bond_atoms):
            pdbqt_lines.append('REMARK%5.0i  A    between atoms: _%i  and  _%i'
                               % (i + 1, a1 + 1, a2 + 1))

        frags = list(Chem.GetMolFrags(mol_rigid_frags))

        #list frag  from which bonds ?
        fg_bonds=[]
        fg_num_rotbonds={}
        for fg in frags:
            tmp_bonds=[]
            for a1,a2 in bond_atoms:
                if a1 in fg or a2 in fg:
                    tmp_bonds.append(mol.GetBondBetweenAtoms(a1, a2).GetIdx())
            if tmp_bonds:
                fg_bonds.append(tmp_bonds)
            else:
                fg_bonds.append(None)
            fg_num_rotbonds[fg] = len(tmp_bonds)

        # frag with long branch ?
        fg_bigbranch={}
        for i, fg_bond in enumerate(fg_bonds):
            tmp_bigger=0
            frag_i_mol=frags[i]
            if fg_bond != None: # for rigid mol
                tmp_frags= Chem.FragmentOnBonds(mol, fg_bond, addDummies=False)
                tmp_frags_list=list(Chem.GetMolFrags(tmp_frags))
                for tmp_frag_j in tmp_frags_list:
                    len_tmp_fg_j=len(tmp_frag_j)
                    if frag_i_mol == tmp_frag_j:
                        pass
                    else:
                        if len_tmp_fg_j > tmp_bigger:
                            tmp_bigger=len_tmp_fg_j
                #print(f'REMARK FRAG: {i} : {len(frags[i])} : {tmp_bigger} ')
            fg_bigbranch[frags[i]] = tmp_bigger

        def weigh_frags(frag):
            return fg_bigbranch[frag], -fg_num_rotbonds[frag],   # bond_weight
        frags = sorted(frags, key=weigh_frags)
####By Zou Rongfeng, Xu Ximing
        def count_match(a, b):
            n = 0
            for i in a:
                for j in b:
                    if int(i) == int(j):
                        n = n + 1
            return n
        match_num = []
        for i in range(len(frags)):
            x = frags[i]
            match_num.append(count_match(x, core_atom_idx_list))
        pop_num = match_num.index(max(match_num))

####By Zou Rongfeng 
        # Start writting the lines with ROOT
        pdbqt_lines.append('ROOT')
        frag = frags.pop(pop_num)

        for idx in frag:
            pdbqt_lines.append(atom_lines[idx])
        pdbqt_lines.append('ENDROOT')

        # Now build the tree of torsions usign DFS algorithm. Keep track of last
        # route with following variables to move down the tree and close branches
        branch_queue = []
        current_root = frag
        old_roots = [frag]

        visited_frags = []
        visited_bonds = []
        while len(frags) > len(visited_frags):
            end_branch = True
            for frag_num, frag in enumerate(frags):
                for bond_num, (a1, a2) in enumerate(bond_atoms):
                    if (frag_num not in visited_frags and
                        bond_num not in visited_bonds and
                        (a1 in current_root and a2 in frag or
                         a2 in current_root and a1 in frag)):
                        # direction of bonds is important
                        if a1 in current_root:
                            bond_dir = '%i %i' % (a1 + 1, a2 + 1)
                        else:
                            bond_dir = '%i %i' % (a2 + 1, a1 + 1)
                        pdbqt_lines.append('BRANCH %s' % bond_dir)
                        for idx in frag:
                            pdbqt_lines.append(atom_lines[idx])
                        branch_queue.append('ENDBRANCH %s' % bond_dir)

                        # Overwrite current root and stash previous one in queue
                        old_roots.append(current_root)
                        current_root = frag

                        # remove used elements from stack
                        visited_frags.append(frag_num)
                        visited_bonds.append(bond_num)

                        # mark that we dont want to end branch yet
                        end_branch = False
                        break
                    else:
                        continue
                    break  # break the outer loop as well

            if end_branch:
                pdbqt_lines.append(branch_queue.pop())
                if old_roots:
                    current_root = old_roots.pop()
        # close opened branches if any is open
        while len(branch_queue):
            pdbqt_lines.append(branch_queue.pop())
        pdbqt_lines.append('TORSDOF %i' % num_torsions)
    else:
        pdbqt_lines.extend(atom_lines)

    return '\n'.join(pdbqt_lines)

def MolToPDBQTBlock(mol, flexible=True, addHs=False, computeCharges=False):
    core_atom_idx_list = []
    return MolCoreToPDBQTBlock(mol, core_atom_idx_list, flexible, addHs, computeCharges)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog='rdkit2pdbqt.py', description='For example: rdkit2pdbqt.py -l ligand.sdf scaffold.sdf')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-r', '--receptor', help='receptor file in pdb format')
    group.add_argument('-l', '--ligand', help='ligand file in sdf format, appending a core sdf file if available', nargs='+')
    args = parser.parse_args()

    if args.receptor:
        receptor_mol=Chem.MolFromPDBFile(sys.argv[2], removeHs=False)
        receptor_lines=MolToPDBQTBlock(receptor_mol, False, False, True)
        print(receptor_lines)
        exit()
    elif args.ligand:
        lig_and_core = args.ligand
        mol=Chem.MolFromMolFile(lig_and_core[0],removeHs=False, sanitize = False)
        mol_problems = Chem.DetectChemistryProblems(mol)
        if len(mol_problems) > 0:
            for problem in mol_problems:
                if "N, 4" in problem.Message():
                    at_idx = problem.GetAtomIdx()
                    atom = mol.GetAtomWithIdx(at_idx)
                    chg = atom.GetFormalCharge()
                    print(f'REMARK    N {at_idx} with formal charge {chg}')
                    atom.SetFormalCharge(1)
                    atom.UpdatePropertyCache()
                else:
                    print(problem.Message())
                    exit()
        Chem.SanitizeMol(mol)
        core_atom_idx_list = []
        if len(lig_and_core) == 2:
            #scaffold=Chem.MolFromMolFile(lig_and_core[1],removeHs=False)
            #For better solution: https://github.com/proteneer/timemachine/blob/master/timemachine/fe/atom_mapping.py
            #mcs = Chem.MolFromSmarts(rdFMCS.FindMCS([mol, scaffold]).smartsString)
            #scaffold_atoms = list(mol.GetSubstructMatch(mcs))
            #ref_atom_mapping_list = list(scaffold.GetSubstructMatch(mcs))
            
            #We treat the lig_and_core as the scaffold to reduce the MCS searching time.
            scaffold=Chem.MolFromMolFile(lig_and_core[1],removeHs=False)
            mcs= Chem.RemoveHs(scaffold)
            scaffold_atoms = list(mol.GetSubstructMatch(mcs))
            ref_atom_mapping_list = list(scaffold.GetSubstructMatch(mcs))
            core_atom_mapping_dict = dict(zip(scaffold_atoms, ref_atom_mapping_list))
            core_atom_idx_list = get_core_alignment_for_template_docking(scaffold, mol, core_atom_mapping_dict)		

        pdbqtlines=MolCoreToPDBQTBlock(mol, core_atom_idx_list, True, False, True)
        print(pdbqtlines)
        exit()            
