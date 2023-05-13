#!/usr/bin/env python3
'''
Copyright 2015, Maciej Wojcikowski Revision caf5d84a.
Wojcikowski, M., Zielenkiewicz, P. & Siedlecki, P. 
Open Drug Discovery Toolkit (ODDT): a new open-source player in the drug discovery field. 
J Cheminform 7, 26 (2015). https://doi.org/10.1186/s13321-015-0078-2

Modified by XU Ximing
xuximing@ouc.edu.cn
2023.04.22
'''
from __future__ import absolute_import, print_function
from math import isnan, isinf
from itertools import combinations

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondType
import sys

METALS = (3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
          37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57,
          58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74,
          75, 76, 77, 78, 79, 80, 81, 82, 83, 87, 88, 89, 90, 91, 92, 93, 94,
          95, 96, 97, 98, 99, 100, 101, 102, 103)


def PathFromAtomList(mol, amap):
    out = []
    for i, j in combinations(amap, 2):
        bond = mol.GetBondBetweenAtoms(i, j)
        if bond:
            out.append(bond.GetIdx())
    return out



def AtomListToSubMol(mol, amap, includeConformer=False):
    """
    Parameters
    ----------
        mol: rdkit.Chem.rdchem.Mol
            Molecule
        amap: array-like
            List of atom indices (zero-based)
        includeConformer: bool (default=True)
            Toogle to include atoms coordinates in submolecule.

    Returns
    -------
        submol: rdkit.Chem.rdchem.RWMol
            Submol determined by specified atom list
    """
    if not isinstance(amap, list):
        amap = list(amap)
    submol = Chem.RWMol(Chem.Mol())
    for aix in amap:
        submol.AddAtom(mol.GetAtomWithIdx(aix))
    for i, j in combinations(amap, 2):
        bond = mol.GetBondBetweenAtoms(i, j)
        if bond:
            submol.AddBond(amap.index(i),
                           amap.index(j),
                           bond.GetBondType())
    if includeConformer:
        for conf in mol.GetConformers():
            new_conf = Chem.Conformer(len(amap))
            for i in range(len(amap)):
                new_conf.SetAtomPosition(i, conf.GetAtomPosition(amap[i]))
                new_conf.SetId(conf.GetId())
                new_conf.Set3D(conf.Is3D())
            submol.AddConformer(new_conf)
    return submol



def MolFromPDBBlock(molBlock,
                    sanitize=True,
                    removeHs=True,
                    flavor=0):
    # before 2019.03 pre-sanitization is required
    pre_sanitize = False
    if sanitize and rdkit.__version__ <= '2018.09':
        pre_sanitize = True

    mol = Chem.MolFromPDBBlock(molBlock,
                               sanitize=pre_sanitize,
                               removeHs=removeHs,
                               flavor=flavor)
    if mol is None:
        return None
    # Adjust connectivity
    for atom in mol.GetAtoms():
        res = atom.GetPDBResidueInfo()
        if res is None:
            continue
        res_name = res.GetResidueName()
        atom_name = res.GetName().strip()

        # Fix missing double bonds in RDKit - double bonds
        if atom_name == 'O' and not res.GetIsHeteroAtom() and atom.GetDegree() == 1:
            atom.SetNoImplicit(True)
            atom.GetBonds()[0].SetBondType(Chem.BondType.DOUBLE)

        # Double bonds in sidechains
        if res_name in ['HID', 'HIE', 'HIP']:
            if atom_name == 'CD2':
                for bond in atom.GetBonds():
                    if bond.GetOtherAtom(atom).GetPDBResidueInfo().GetName().strip() == 'CG':
                        bond.SetBondType(Chem.BondType.DOUBLE)
                        break
            if res_name == 'HID':
                if atom_name == 'CE1':
                    for bond in atom.GetBonds():
                        if bond.GetOtherAtom(atom).GetPDBResidueInfo().GetName().strip() == 'ND1':
                            bond.SetBondType(Chem.BondType.DOUBLE)
                            break
            elif res_name in ['HIE', 'HIP']:
                if atom_name == 'CE1':
                    for bond in atom.GetBonds():
                        if bond.GetOtherAtom(atom).GetPDBResidueInfo().GetName().strip() == 'NE2':
                            bond.SetBondType(Chem.BondType.DOUBLE)
                            break
    mol.UpdatePropertyCache(strict=sanitize)

    # Set metal coordination (zero order) bond orders to single to prevent adding Hs
    if rdkit.__version__ >= '2018.03':
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.ZERO:
                a1 = bond.GetBeginAtom()
                a2 = bond.GetEndAtom()
                # single bonds only if there are enough electrons
                if ((a1.GetAtomicNum() in METALS and
                     a2.GetNumImplicitHs() + a2.GetNumExplicitHs() > 0) or
                    (a2.GetAtomicNum() in METALS and
                     a1.GetNumImplicitHs() + a1.GetNumExplicitHs() > 0)):
                    bond.SetBondType(Chem.BondType.SINGLE)

    if sanitize:
        result = Chem.SanitizeMol(mol)
        if result != 0:
            return None

    # Debug
    # for atom in mol.GetAtoms():
    #     res = atom.GetPDBResidueInfo()
    #     if res is None:
    #         continue
    #     res_name = res.GetResidueName()
    #     atom_name = res.GetName().strip()
    #     if atom_name in ['NE2', 'ND1'] and res_name in ['HID', 'HIE', 'HIS']:
    #         print(res_name,
    #               atom_name,
    #               atom.GetDegree(),
    #               atom.GetTotalValence(),
    #               atom.GetNumExplicitHs(),
    #               atom.GetNumImplicitHs(),
    #               sum(n.GetAtomicNum() == 1 for n in atom.GetNeighbors()),
    #               sep='\t')

    return mol



# Mol2 Atom typing
def _sybyl_atom_type(atom):
    """ Asign sybyl atom type
    Reference #1: http://www.tripos.com/mol2/atom_types.html
    Reference #2: http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf
    """
    sybyl = None
    atom_symbol = atom.GetSymbol()
    atomic_num = atom.GetAtomicNum()
    hyb = atom.GetHybridization()-1  # -1 since 1 = sp, 2 = sp1 etc
    hyb = min(hyb, 3)
    degree = atom.GetDegree()
    aromtic = atom.GetIsAromatic()

    # define groups for atom types
    guanidine = '[NX3,NX2]([!O,!S])!@C(!@[NX3,NX2]([!O,!S]))!@[NX3,NX2]([!O,!S])'  # strict
    # guanidine = '[NX3]([!O])([!O])!:C!:[NX3]([!O])([!O])' # corina compatible
    # guanidine = '[NX3]!@C(!@[NX3])!@[NX3,NX2]'
    # guanidine = '[NX3]C([NX3])=[NX2]'
    # guanidine = '[NX3H1,NX2,NX3H2]C(=[NH1])[NH2]' # previous
    #

    if atomic_num == 6:
        if aromtic:
            sybyl = 'C.ar'
        elif degree == 3 and _atom_matches_smarts(atom, guanidine):
            sybyl = 'C.cat'
        else:
            sybyl = '%s.%i' % (atom_symbol, hyb)
    elif atomic_num == 7:
        if aromtic:
            sybyl = 'N.ar'
        elif _atom_matches_smarts(atom, 'C(=[O,S])-N'):
            sybyl = 'N.am'
        elif degree == 3 and _atom_matches_smarts(atom, '[$(N!-*),$([NX3H1]-*!-*)]'):
            sybyl = 'N.pl3'
        elif _atom_matches_smarts(atom, guanidine):  # guanidine has N.pl3
            sybyl = 'N.pl3'
        elif degree == 4 or hyb == 3 and atom.GetFormalCharge():
            sybyl = 'N.4'
        else:
            sybyl = '%s.%i' % (atom_symbol, hyb)
    elif atomic_num == 8:
        # http://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
        if degree == 1 and _atom_matches_smarts(atom, '[CX3](=O)[OX1H0-]'):
            sybyl = 'O.co2'
        elif degree == 2 and not aromtic:  # Aromatic Os are sp2
            sybyl = 'O.3'
        else:
            sybyl = 'O.2'
    elif atomic_num == 16:
        # http://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
        if degree == 3 and _atom_matches_smarts(atom, '[$([#16X3]=[OX1]),$([#16X3+][OX1-])]'):
            sybyl = 'S.O'
        # https://github.com/rdkit/rdkit/blob/master/Data/FragmentDescriptors.csv
        elif _atom_matches_smarts(atom, 'S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])(-[#6])-[#6]'):
            sybyl = 'S.o2'
        else:
            sybyl = '%s.%i' % (atom_symbol, hyb)
    elif atomic_num == 15 and hyb == 3:
        sybyl = '%s.%i' % (atom_symbol, hyb)

    if not sybyl:
        sybyl = atom_symbol
    return sybyl


def _atom_matches_smarts(atom, smarts):
    idx = atom.GetIdx()
    patt = Chem.MolFromSmarts(smarts)
    for m in atom.GetOwningMol().GetSubstructMatches(patt):
        if idx in m:
            return True
    return False


def _amide_bond(bond):
    a1 = bond.GetBeginAtom()
    a2 = bond.GetEndAtom()
    if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7 or
            a2.GetAtomicNum() == 6 and a1.GetAtomicNum() == 7):
        # https://github.com/rdkit/rdkit/blob/master/Data/FragmentDescriptors.csv
        patt = Chem.MolFromSmarts('C(=O)-N')
        for m in bond.GetOwningMol().GetSubstructMatches(patt):
            if a1.GetIdx() in m and a2.GetIdx() in m:
                return True
    return False


def MolFromPDBQTBlock(block, sanitize=True, removeHs=True):
    """Read PDBQT block to a RDKit Molecule

    Parameters
    ----------
        block: string
            Residue name which explicitly pint to a ligand(s).
        sanitize: bool (default=True)
            Should the sanitization be performed
        removeHs: bool (default=True)
            Should hydrogens be removed when reading molecule.

    Returns
    -------
        mol: rdkit.Chem.rdchem.Mol
            Molecule read from PDBQT
    """
    pdb_lines = []
    name = ''
    data = {}
    for line in block.split('\n'):
        # Get all know data from REMARK section
        if line[:12] == 'REMARK  Name':
            name = line[15:].strip()
        elif line[:18] == 'REMARK VINA RESULT':
            tmp = line[19:].split()
            data['vina_affinity'] = tmp[0]
            data['vina_rmsd_lb'] = tmp[1]
            data['vina_rmsd_ub'] = tmp[2]

        # no more data to collect
        if line[:4] != 'ATOM':
            continue

        pdb_line = line[:56]
        pdb_line += '1.00  0.00           '

        # Do proper atom type lookup
        atom_type = line[71:].split()[1]
        if atom_type == 'A':
            atom_type = 'C'
        elif atom_type[:1] == 'O':
            atom_type = 'O'
        elif atom_type[:1] == 'H':
            atom_type = 'H'
            if removeHs:
                continue
        elif atom_type == 'NA':
            atom_type = 'N'
        elif atom_type == 'SA':
            atom_type = 'S'

        pdb_lines.append(pdb_line + atom_type)
    mol = Chem.MolFromPDBBlock('\n'.join(pdb_lines), sanitize=False)
    if sanitize:
        Chem.SanitizeMol(mol)
    else:
        Chem.GetSSSR(mol)
    # reorder atoms using serial
    new_order = sorted(range(mol.GetNumAtoms()),
                       key=lambda i: (mol.GetAtomWithIdx(i)
                                      .GetPDBResidueInfo()
                                      .GetSerialNumber()))
    mol = Chem.RenumberAtoms(mol, new_order)

    # properties must be set on final copy of Mol, RenumberAtoms purges data
    mol.SetProp('_Name', name)
    for k, v in data.items():
        mol.SetProp(str(k), str(v))

    return mol



def PDBQTAtomLines(mol, donors, acceptors):
    """Create a list with PDBQT atom lines for each atom in molecule. Donors
    and acceptors are given as a list of atom indices.
    """

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
            pdbqt_line += 'A'
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
            pdbqt_line += atom.GetSymbol()
        pdbqt_lines.append(pdbqt_line)
    return pdbqt_lines



def MolToPDBQTBlock(mol, flexible=True, addHs=False, computeCharges=False):
    """Write RDKit Molecule to a PDBQT block

    Parameters
    ----------
        mol: rdkit.Chem.rdchem.Mol
            Molecule with a protein ligand complex
        flexible: bool (default=True)
            Should the molecule encode torsions. Ligands should be flexible,
            proteins in turn can be rigid.
        addHs: bool (default=False)
            The PDBQT format requires at least polar Hs on donors. By default Hs
            are added.
        computeCharges: bool (default=False)
            Should the partial charges be automatically computed. If the Hs are
            added the charges must and will be recomputed. If there are no
            partial charge information, they are set to 0.0.

    Returns
    -------
        block: str
            String wit PDBQT encoded molecule
    """
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

    # vina scores
    if (mol.HasProp('vina_affinity') and mol.HasProp('vina_rmsd_lb') and mol.HasProp('vina_rmsd_lb')):
        pdbqt_lines.append('REMARK VINA RESULT:  ' +
                           ('%.1f' % float(mol.GetProp('vina_affinity'))).rjust(8) +
                           ('%.3f' % float(mol.GetProp('vina_rmsd_lb'))).rjust(11) +
                           ('%.3f' % float(mol.GetProp('vina_rmsd_ub'))).rjust(11))

    pdbqt_lines.append('REMARK  Name = ' + (mol.GetProp('_Name') if mol.HasProp('_Name') else ''))
    if flexible:
        # Find rotatable bonds
        '''
        rot_bond = Chem.MolFromSmarts('[!$(*#*)&!D1&!$(C(F)(F)F)&'
                                      '!$(C(Cl)(Cl)Cl)&'
                                      '!$(C(Br)(Br)Br)&'
                                      '!$(C([CH3])([CH3])[CH3])&'
                                      '!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&'
                                      '!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&'
                                      '!$([CD3](=[N+])-!@[#7!D1])&'
                                      '!$([#7!D1]-!@[CD3]=[N+])]-!@[!$(*#*)&'
                                      '!D1&!$(C(F)(F)F)&'
                                      '!$(C(Cl)(Cl)Cl)&'
                                      '!$(C(Br)(Br)Br)&'
                                      '!$(C([CH3])([CH3])[CH3])]')
        '''
        
        #rot_bond = Chem.MolFromSmarts('[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]') # From Chemaxon
        rot_bond  = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]') #single and not ring, really not in ring?
        
        exclude_list = ['[NX3]-[CX3]=[O,N]']
        #, '[!#1]-[C;H3,Cl3,F3,Br3]']
        excluded_atoms = []
        for excluded_smarts in exclude_list:
            excluded_bond = Chem.MolFromSmarts(excluded_smarts)
            excluded_atoms += [(x[0],x[1]) for x in list(mol.GetSubstructMatches(excluded_bond))]

        tertiary_amide_bonds = Chem.MolFromSmarts('[NX3]([!#1])([!#1])-[CX3]=[O,N]')

        bond_atoms = list(mol.GetSubstructMatches(rot_bond))

        tertiary_amide_bond_atoms=[(x[0],x[3]) for x in list(mol.GetSubstructMatches(tertiary_amide_bonds))]
        for excld_atom in excluded_atoms:
            #bond_atoms.remove(amide_bond_atom)
            excld_atom_reverse=(excld_atom[1],excld_atom[0])
            if excld_atom in bond_atoms:
                bond_atoms.remove(excld_atom)
            elif excld_atom_reverse in bond_atoms:
                bond_atoms.remove(excld_atom_reverse)

        if len(tertiary_amide_bond_atoms) > 0:
            for tertiary_amide_bond_atom in tertiary_amide_bond_atoms:
                bond_atoms.append(tertiary_amide_bond_atom)


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
                #else:
                #    tmp_bigger= max(len(tmp_frags_list[0]), len(tmp_frags_list[1]))
                #mol.GetBonds()[b_index].SetProp("large_part", str(tmp_bigger))
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


        #def weigh_frags(frag):
            """sort by the fragment size and the number of bonds (secondary)"""
            #num_bonds = 0
            # bond_weight = 0
            #big_frag_size=0
            #for a1, a2 in bond_atoms:
            #    if a1 in frag or a2 in frag:
            #        num_bonds += 1
                    #big_frag_size = max(big_frag_size, int(mol.GetBondBetweenAtoms(a1, a2).GetProp("large_part")))
            # changed signs are fixing mixed sorting type (ascending/descending)
            #return -num_bonds, -len(frag),  # bond_weight
        def weigh_frags(frag):
            return fg_bigbranch[frag], -fg_num_rotbonds[frag],   # bond_weight
        frags = sorted(frags, key=weigh_frags)

        # Start writting the lines with ROOT
        pdbqt_lines.append('ROOT')
        frag = frags.pop(0)
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

def HelpAndExit():
    print("USAGE: rdkit2pdbqt.py -r protein.pdb")
    print("USAGE: rdkit2pdbqt.py -l ligand.sdf")
    exit()

if __name__ == "__main__":
    #parsing the command
    if len(sys.argv) != 3:
        HelpAndExit()

    if sys.argv[1] == '-r':
        receptor_mol=Chem.MolFromPDBFile(sys.argv[2], removeHs=False)
        #AllChem.AssignBondOrdersFromTemplate(receptor_mol,receptor_mol)
        # detect aromaticity
        #AllChem.AssignAtomChiralTagsFromStructure(receptor_mol)
        #AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmarts('[a]'), receptor_mol)
        #Chem.AssignStereochemistry(receptor_mol)

        receptor_lines=MolToPDBQTBlock(receptor_mol, False, False, True)
        print(receptor_lines)
        exit()

    elif sys.argv[1] == '-l':
        mol=Chem.MolFromMolFile(sys.argv[2],removeHs=False, sanitize=False)
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
        pdbqtlines=MolToPDBQTBlock(mol, True, False, True)
        print(pdbqtlines)
        exit()

    else:
        HelpAndExit()
