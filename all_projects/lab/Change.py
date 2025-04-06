import shutil
import copy
import os
import numpy as np
import math
import urllib.request

other_atoms = ['H', 'F', 'CL', 'BR', 'I']
hetero_atoms = {'N': 14, 'S': 32, 'O': 16, 'P': 31}
pi_heteroatoms = ['N', 'S', 'O']

interest_atom = 'F'

F_C_bond = 1.4
H_C_bond = 1.1

r = 1.2


def download_pdb(ligand_name):
    urllib.request.urlretrieve('http://files.rcsb.org/download/' + ligand_name + '.sdf', ligand_name + '.sdf')


def dict_creator(ligand_name): #создает словарь с атомами, координатами, соседями и кч из прочитанного файла
    file = open(ligand_name, 'r')
    lines = file.readlines()
    inf = lines[3].split()
    num_atoms, num_bonds = int(inf[0]), int(inf[1])
    xyz_lines = lines[4:num_atoms + 4]
    bonds_lines = lines[(num_atoms + 4):(num_atoms + num_bonds + 4)]
    n = 1
    atoms = {}
    int_atoms = []
    bonds = {}
    chir_list = []
    double_bonds = []
    for line in xyz_lines:
        line = line.split()
        x, y, z, atom, chir = line[0], line[1], line[2], line[3], line[6]
        x, y, z = float(x), float(y), float(z)
        chir = int(chir)
        if chir != 0:
            chir_list[n] = chir
        atoms[n] = [atom, [x, y, z], [], chir]
        if atom in hetero_atoms:
            atoms[n][3] += hetero_atoms[atom]
        if atom == interest_atom:
            int_atoms.append(n)
        bonds[n] = []
        n += 1
    for line in bonds_lines:
        line = line.split()
        at_1, at_2, n = int(line[0]), int(line[1]), int(line[2])
        bonds[at_1].append(at_2)
        bonds[at_2].append(at_1)
        atoms[at_1][2].append(at_2)
        atoms[at_2][2].append(at_1)
        atoms[at_1][3] += n * 100
        atoms[at_2][3] += n * 100
        if n == 2:
            double_bonds.append(sorted([at_1, at_2]))
    return atoms, int_atoms, chir_list, bonds, double_bonds


def ident_atoms(atoms): #находит атомы с одинаковыми кч
    cn = [atoms[i][3] for i in atoms]
    rep_cn = []
    for n in rep(cn):
        if n not in rep_cn:
            rep_cn.append(n)
    sim_atoms = []
    for cn in rep_cn:
        numb = []
        for key, val in atoms.items():
            if cn in val:
                numb.append(key)
        sim_atoms.append(numb)
    return sim_atoms


def rep(seq): #находит повторы в списке и удаляет их
    seen = set()
    seen_add = seen.add
    return [x for x in seq if (x in seen or seen_add(x))]


def gen_dir(fold): #создает папку
    if os.path.exists(fold):
        pass
    else:
        os.makedirs(fold)


def isoms(atoms_list, bond, int_atoms, bonds, coords, flag=''): #определяет геометрию двойной связи
    root_list = []
    for parent_atom in bond:
        if flag == '':
            child_atoms = [atom for atom in atoms_list[parent_atom][2] if atom not in bond]
            if child_atoms == []:
                return None
            child_num = [atoms_list[atom][3] for atom in child_atoms]
            if len(child_num) == 2:
                if child_num[0] == child_num[1]:
                    print('112')
                    return None
        else:
            child_atoms = []
            for child in bonds[parent_atom]:
                if child in int_atoms or child in atoms_list:
                    child_atoms.append(child)
            child_num = [0] * len(child_atoms)
            for i in range(len(child_atoms)):
                if child_atoms[i] not in int_atoms:
                    child_num[i] = atoms_list[child_atoms[i]][3]
            for i in range(len(child_atoms)):
                if child_num[i] == 0:
                    child_num[i] = max(child_num) + 1000
        major_atom_index = child_num.index(max(child_num))
        root_list.append([parent_atom, child_atoms[major_atom_index]])
    v = []
    for pair in root_list:
        parent_atom = pair[0]
        child_atom = pair[1]
        vect = np.array(coords[child_atom - 1]) - np.array(coords[parent_atom - 1])
        v.append(vect)
    leng_1 = np.linalg.norm(v[0])
    leng_2 = np.linalg.norm(v[1])
    cos = (np.dot(v[0], v[1]))/(leng_1 * leng_2)
    angle = float(np.arccos(cos))
    angle = angle * 180 / math.pi
    if angle < 100:
        atoms_list[bond[0]][3] += 13
        atoms_list[bond[1]][3] += 13
        return 'E'
    else:
        return 'S'


ligands = {'5dda': '59M'}

for struct in ligands:
    filename = ligands[struct] + '.mol'
    name = ligands[struct]
    print(name)
    atoms, int_atoms, chir_list, bonds, double_bonds = dict_creator(filename)
    gen_dir('sim_atoms/' + name)
    shutil.copy(filename, 'sim_atoms/' + name)
    gen_dir('change/' + name)
    shutil.copy(filename, 'change/' + name)
    for i in atoms:  # убрали хиральность, связанную с наличием фтора
        if i in chir_list:
            H = 0
            F = 0
            for n in bonds[i]:
                if n not in atoms and n not in int_atoms:
                    H += 1
                if n in int_atoms:
                    F += 1
            if H == 1 and F == 1:
                atoms[i][3] -= chir_list[i]
    for i in range(len(atoms)):   #сам алгоритм ч1
        atoms_copy = copy.deepcopy(atoms)
        for p in atoms_copy:
            neig = atoms_copy[p][2]
            for n in neig:
                atoms_copy[p][3] += atoms[n][3]
        atoms = atoms_copy

    sim_atom = ident_atoms(atoms)  # список с парами симметричных атомов
    sim_H = []
    sim_C = []
    for pair in sim_atom:
        if atoms[pair[0]][0] in other_atoms:
            sim_H.append(pair)
        else:
            sim_C.append(pair)

    print(sim_C)
    print(sim_H)