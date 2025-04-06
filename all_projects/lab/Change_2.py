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


def download_sdf(ligand_name): #скачивание sdf файла лиганда
    urllib.request.urlretrieve('https://files.rcsb.org/ligands/download/' + ligand_name + '_model.sdf', ligand_name + '.sdf')


def read_file(file_name): #читает mol-файла, создает словарь с атомами, координатами, соседями и кч
    file_name = r'C:/Users/Дмитрий/Desktop/По лабе/rot/' + file_name
    file = open(file_name, 'r')
    lines = file.readlines()
    inf = lines[3].split()
    num_atoms, num_bonds = int(inf[0]), int(inf[1])
    xyz = lines[4:num_atoms+4]
    bonds = lines[(num_atoms+4):(num_atoms+num_bonds + 4)]
    n = 1
    atoms = {}
    graph = {}
    atoms_coord = []
    int_atoms = []
    all_bonds = {}
    chir_list = {}
    double_bonds = []
    for line in xyz:
        line = line.split()
        x, y, z, atom, chir = line[0], line[1], line[2], line[3], line[6]
        if atom not in other_atoms:
            x, y, z = float(x), float(y), float(z)
            chir = int(chir)
            if chir != 0:
                chir_list[n] = chir
            atoms[n] = [atom, [x, y, z], [], chir]
            graph[n] = []
        if atom in hetero_atoms:
            atoms[n][3] += hetero_atoms[atom]
        if atom == interest_atom:
            int_atoms.append(n)
        atoms_coord.append([float(x), float(y), float(z)])
        all_bonds[n] = []
        n += 1
    for line in bonds:
        line = line.split()
        at_1, at_2, n = int(line[0]), int(line[1]), int(line[2])
        if at_1 in atoms and at_2 in atoms:
            atoms[at_1][2].append(at_2)
            graph[at_1].append(at_2)
            atoms[at_2][2].append(at_1)
            graph[at_2].append(at_1)
            atoms[at_1][3] += n * 100
            atoms[at_2][3] += n * 100
        if n == 2:
            double_bonds.append([at_1, at_2])
        all_bonds[at_1].append(at_2)
        all_bonds[at_2].append(at_1)
    return atoms, double_bonds, graph, int_atoms, atoms_coord, all_bonds, chir_list #


def neighbors(val): #просто выводит соседей
    return val[2]


def element(val): #просто выводит символ элемента
    return val[0]


def rep(seq): #находит повторы в списке и удаляет их
    seen = set()
    seen_add = seen.add
    return [x for x in seq if (x in seen or seen_add(x))]


def gen_dir(fold): #создает папку
    if os.path.exists(fold):
        pass
    else:
        os.makedirs(fold)


def file_new(name, atoms_list, pair, number): #создает папку с исходным файлом и парами симметричных атомов
    file = open("sim_atoms/" + name + "/" + name + "_" + str(number) + ".xyz", "w", encoding="utf-8")
    file.write(str(len(atoms_list)) + "\n\n")
    for atom in pair:
        file.write(atoms_list[atom][0] + "\t")
        coords = [str(i) for i in atoms_list[atom][1]]
        for c in coords:
            file.write(c + "\t")
        file.write("\n")
    file.close()


def isoms(atoms_list, bond, int_atoms, bonds, coords, C_atoms, flag=''): #определяет геометрию двойной связи
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
                if child in int_atoms or child in C_atoms:
                    child_atoms.append(child)
            child_num = [0] * len(child_atoms)
            for i in range(len(child_atoms)):
                if child_atoms[i] not in int_atoms:
                    child_num[i] = C_atoms[child_atoms[i]][3]
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


def dfs_paths(C_atoms, point_1, point_2): #находит путь между парой атомов
    graph = {key: val[2] for key, val in C_atoms.items()}
    stack = [(point_1, [point_1])]
    while stack:
        (vertex, path) = stack.pop()
        for next in set(graph[vertex]) - set(path):
            if next == point_2:
                yield path + [next]
            else:
                stack.append((next, path + [next]))


def joint_double_bond(double_bonds, C_atoms, pair): #определяет геометрию двойной связи, от одного атома которой расходятся цепи, содержащие потенциально симметричные атомы
    paths = list(dfs_paths(C_atoms, pair[0], pair[1]))
    not_suitable = []
    for path in paths:
        l = len(path)
        if l%2 == 0:
            continue
        center_atom = path[l//2]
        int_atom = 0
        for atom in C_atoms[center_atom][2]:
            if (atom not in path) and (sorted([center_atom, atom]) in double_bonds):
                int_atom = atom
        if int_atom == 0:
            continue
        branches = C_atoms[int_atom][2].copy()
        branches.remove(center_atom)
        neig = C_atoms[center_atom][2].copy()
        neig.remove(int_atom)
        if len(branches) == 2 and (C_atoms[branches[0]][3] != C_atoms[branches[1]][3]):
            not_suitable.append(neig)
    return not_suitable


def ident_atoms(C_atoms): #находит атомы с одинаковыми кч
    cn = [C_atoms[i][3] for i in C_atoms]
    rep_cn = rep(cn)
    sim_atom = []
    for cn in rep_cn:
        numb = []
        for key, val in C_atoms.items():
            if cn in val:
                numb.append(key)
        sim_atom.append(numb)
    return sim_atom


def chain_double_bonds(C_atoms, pair, double_bonds, int_atoms, bonds, atoms_coord): #определяет геометрию двойных связей на пути между потенциально симметричными атомами
    paths = list(dfs_paths(C_atoms, pair[0], pair[1]))
    if len(paths) == 1:
        path = paths[0]
        for i in range(len(path)):
            if i == 0:
                continue
            bond = sorted([path[i], path[i-1]])
            if bond in double_bonds:
                isoms(C_atoms, bond, int_atoms, bonds, atoms_coord, C_atoms)


def find_cycle(graph, k, v, rev=1):   # ищет все циклы в молекуле
    global cycle_atoms2
    if k not in v:
        v.append(k)
    for i in graph[k]:
        if i not in v:
            v.append(i)
        else:
            if i == rev:
                continue
            else:
                cycle_atoms2.append(i) if i not in cycle_atoms2 else None
                continue
        v = find_cycle(graph, i, v, k)
    return list(v)


def non_cycle_double_bonds(cycle_atoms, double_bonds, atoms_list, int_atoms, bonds, atoms_coord, C_atoms): #определяет геометрию двойных связей вне циклов
    for bond in double_bonds:
        if bond[0] not in cycle_atoms and bond[1] not in cycle_atoms:
            isoms(atoms_list, bond, int_atoms, bonds, atoms_coord, C_atoms)


def one_root(atoms_list, pair, C_atoms): #определяет геометрию атома, от которого отходят две цепи с потенциально симметричными атомами
    root = 0
    for i in atoms_list[pair[0]][2]:
        if i in atoms_list[pair[1]][2]:
            root = i
    if root == 0:
        return None
    root_root = [i for i in C_atoms[root][2] if i not in pair]
    if len(root_root) == 2:
        if atoms_list[root_root[0]][3] != atoms_list[root_root[1]][3] and atoms_list[root_root[0]][0] == atoms_list[root_root[1]][0]:
            atoms_list[pair[0]][3] += 2


def change_bond(el, atom, root, atoms_coord): #замена атома
    atom_coord = atoms_coord[atom - 1]
    root_coord = atoms_coord[root - 1]
    vec_bond = np.array(atom_coord) - np.array(root_coord)
    l_vec = np.linalg.norm(vec_bond)
    change_vec = 0
    if el == 'F':
        change_vec = vec_bond * (F_C_bond/l_vec) + np.array(root_coord)
    if el == 'H':
        change_vec = vec_bond * (H_C_bond/l_vec) + np.array(root_coord)
    return change_vec


def aromatic_cycle(double_bonds, cycle_atoms, C_atoms, graph): # находит 5- и 6-членные ароматические циклы в молекуле
    help_list = []
    arom_cycles = []
    for i in cycle_atoms:
        f = [j for j in C_atoms[i][2] if j in cycle_atoms]
        for fin in f:
            paths = list(dfs_paths(C_atoms, i, fin))
            help_list.append(paths)
    for i in help_list:
        for path in i:
            if 5 <= len(path) <= 6:
                if sorted(path) not in arom_cycles:
                    arom_cycles.append(sorted(path))

    del_cyc = []
    for cyc in arom_cycles:
        for atom in cyc:
            neig = [i for i in graph[atom] if i in cyc]
            if len(neig) != 2:
                del_cyc.append(cyc)
    arom_cycles = [i for i in arom_cycles if i not in del_cyc]

    all_arom_cyc = []
    pi_bonds = []
    s_pair = []
    for cyc in arom_cycles:
        for atom in cyc:
            for neig in graph[atom]:
                bond = sorted([atom, neig])
                if bond in double_bonds and neig in cyc:
                    pi_bonds.append(atom)
            if atom not in pi_bonds:
                if C_atoms[atom][0] in pi_heteroatoms:
                    s_pair.append(atom)

    for cyc in arom_cycles:
        e = 0
        for atom in cyc:
            if atom in pi_bonds:
                e += 1
            if atom in s_pair:
                e += 2
        n = (e - 2) / 4
        if n % 1 == 0:
            all_arom_cyc.append(cyc)
    return all_arom_cyc  #


def ident_stereo(C_atoms, atom, double_bonds, bonds, int_atoms, atoms_coord): #определяет конфигурацию атома углерода
    n_root_neig = []
    for i in bonds[atom]:
        if i in C_atoms or i in int_atoms:
            n = 0
            if i in C_atoms:
                n = C_atoms[i][3]
            if n not in n_root_neig:
                n_root_neig.append(n)
    double_bonds_atoms = []
    for bond in double_bonds:
        double_bonds_atoms.append(bond[0])
        double_bonds_atoms.append(bond[1])
    n = len(n_root_neig)
    if atom in double_bonds_atoms and 2 <= n <= 3:
        if len(bonds[atom]) <= 3:
            return 'just one'
        bond = []
        for neig in C_atoms[atom][2]:
            if sorted([neig, atom]) in double_bonds:
                bond = sorted([neig, atom])
        atoms_dict = copy.deepcopy(C_atoms)
        conf = isoms(atoms_dict, bond, int_atoms, bonds, atoms_coord, C_atoms, flag='F')
        return conf
    elif 1 <= n <= 2 and atom not in double_bonds_atoms:
        return 'not stereo'
    elif n >= 3 and atom not in double_bonds_atoms:
        neigs = {}
        nums_list = []
        for neig in C_atoms[atom][2]:
            if C_atoms[neig][3] in neigs.keys():
                return 'simple neighbours'
            neigs[C_atoms[neig][3]] = neig
            nums_list.append(C_atoms[neig][3])
        for neig in bonds[atom]:
            if neig in int_atoms:
                neigs[max(nums_list) + 1000] = neig
        neigs = dict(sorted(neigs.items()))
        val = list(neigs.values())
        atom_1 = val[-1]
        atom_2 = val[-2]
        coord_1 = np.array(atoms_coord[atom_1 - 1])
        coord_2 = np.array(atoms_coord[atom_2 - 1])
        root_coord = np.array(atoms_coord[atom - 1])
        vec_1 = coord_1 - root_coord
        vec_2 = coord_2 - root_coord
        cr = np.cross(vec_1, vec_2)
        atom_3 = val[-3]
        coord_3 = np.array(atoms_coord[atom_3 - 1])
        vec_3 = coord_3 - root_coord
        scal = np.dot(cr, vec_3)
        if scal >= 0:
            return 'R'
        else:
            return 'S'


def sim_change(pair, C_atoms, bonds, input_atom, root, int_atoms, double_bonds, atoms_coord): # находит все симметричные атомы(относительно данного)
    neigs_root = bonds[root]
    neigs_root_int = []
    help_list = []
    for i in neigs_root:
        if i in int_atoms:
            int_atoms.remove(i)
            neigs_root_int.append(i)
    int_atoms.append(input_atom)
    root_conf = ident_stereo(C_atoms, root, double_bonds, bonds, int_atoms, atoms_coord)
    int_atoms.remove(input_atom)
    if root_conf == 'not stereo':
        for neig in bonds[root]:
            if neig not in C_atoms:
                help_list.append(neig)
    for atom in pair:
        if atom != root:
            neigs = bonds[atom]
            neigs_int = []
            for i in neigs:
                if i in int_atoms:
                    int_atoms.remove(i)
                    neigs_int.append(i)
            for neig in bonds[atom]:
                neig_conf = ''
                if neig not in C_atoms:
                    if neig not in int_atoms:
                        int_atoms.append(neig)
                        neig_conf = ident_stereo(C_atoms, atom, double_bonds, bonds, int_atoms, atoms_coord)
                        int_atoms.remove(neig)
                    else:
                        neig_conf = ident_stereo(C_atoms, atom, double_bonds, bonds, int_atoms, atoms_coord)
                if neig_conf == root_conf:
                    help_list.append(neig)
            for i in neigs_int:
                int_atoms.append(i)
    for i in neigs_root_int:
        int_atoms.append(i)
    return help_list


def new_file_change(dict_change, lines, bonds, atoms_coord, atom, name): # создает новый файл с заменой
    changes = {}
    lines_1 = copy.deepcopy(lines)
    for F_change in dict_change['F']:
        F_root = bonds[F_change][0]
        F_coord = change_bond('F', F_change, F_root, atoms_coord)
        changes[F_change] = [F_coord, 'F']
    for H_change in dict_change['H']:
        H_root = bonds[H_change][0]
        H_coord = change_bond('H', H_change, H_root, atoms_coord)
        changes[H_change] = [H_coord, 'H']
    for i in changes:
        line = lines_1[i + 3]
        coord, el = changes[i]
        coord = [round(float(i), 4) for i in coord]
        for g in range(len(coord)):
            s_i = str(coord[g])
            n = 0
            dot = len(s_i)
            for j in range(len(s_i)):
                if s_i[j] == '.':
                    dot = j
                elif j > dot:
                    n += 1
            s_i += '0'*(4 - n)
            coord[g] = s_i
        s = line.split(' ')
        s_list = []
        p = 0
        for k in s:
            if k == '':
                p += 1
            else:
                s_list.append(p)
                p = 0
        change_line = ' ' * (s_list[0]) + coord[0] + ' ' * (s_list[1] + 1) + coord[1] + ' ' * (s_list[2] + 1) + coord[2] + ' ' * (s_list[3] + 1) + el + '   0  0  0  0  0  0  0  0  0  0  0  0' + '\n'
        lines_1[i + 3] = change_line
    new_file = open("change/" + name + "/" + name + "_" + str(atom) + ".mol", "w", encoding="utf-8")
    for line in lines_1:
        new_file.write(line)
    new_file.close()


def change(file_name, sim_atoms, C_atoms, bonds, int_atoms, double_bonds, atoms_coord): # ввод атомов, которые нужно заменить(для ввода в консоли)
    file = open(file_name, 'r')
    lines = file.readlines()
    file.close()
    change_atoms = {}
    help_list = []
    for pair in sim_atoms:
        for i in pair:
            help_list.append(i)
    while True:
        i = input('Номер атома:')
        if i == '':
            break
        i = int(i)
        el = input('Поменять на:')
        if el == 'F' and i in int_atoms:
            print('*** Это уже F ***')
            continue
        if el == 'H' and i not in int_atoms:
            print('*** Это уже H ***')
            continue
        if i not in change_atoms:
            change_atoms[i] = el
        print('___')
    for atom in change_atoms:
        root = bonds[atom][0]
        all_change = [atom]
        if root in help_list:
            for pair in sim_atoms:
                if root in pair:
                    sim_pair = pair
                    sim_atom_list = sim_change(sim_pair, C_atoms, bonds, atom, root, int_atoms, double_bonds, atoms_coord)
                    for i in sim_atom_list:
                        all_change.append(i)
        if change_atoms[atom] == 'F':
            for i in all_change:
                change_dict = {'F': [i], 'H': []}
                if i in int_atoms:
                    continue
                for j in all_change:
                    if i != j:
                        if j in int_atoms:
                            change_dict['F'].append(j)
                        else:
                            change_dict['H'].append(j)
                print(change_dict)
                #new_file_change(change_dict, lines, bonds, atoms_coord, i)
        if change_atoms[atom] == 'H':
            for i in all_change:
                change_dict = {'H': [i], 'F': []}
                if i not in int_atoms:
                    continue
                for j in all_change:
                    if i != j:
                        if j in int_atoms:
                            change_dict['H'].append(j)
                        else:
                            change_dict['F'].append(j)
                print(change_dict)
                #new_file_change(change_dict, lines, bonds, atoms_coord, i)
        if root not in help_list:
            H_list = []
            F_list = []
            for neig in bonds[root]:
                if neig not in C_atoms and neig not in int_atoms:
                    H_list.append(neig)
                if neig in int_atoms:
                    F_list.append(neig)
            if change_atoms[atom] == 'F':
                if len(H_list) == 3:
                    for i in H_list:
                        change_dict = {'F': [i], 'H': []}
                        for j in all_change:
                            if i != j:
                                if j in int_atoms:
                                    change_dict['F'].append(j)
                                else:
                                    change_dict['H'].append(j)
                        print(change_dict)
                        #new_file_change(change_dict, lines, bonds, atoms_coord, i)
            if change_atoms[atom] == 'H':
                if len(F_list) == 3:
                    for i in F_list:
                        change_dict = {'H': [i], 'F': []}
                        for j in all_change:
                            if i != j:
                                if j not in int_atoms:
                                    change_dict['H'].append(j)
                                else:
                                    change_dict['F'].append(j)
                        print(change_dict)
                        #new_file_change(change_dict, lines, bonds, atoms_coord, i)


def change_2(sim_atoms, C_atoms, bonds, int_atoms, number_atom, double_bonds, atoms_coord): #поиск симметричных атомов
    help_list = []
    for pair in sim_atoms:
        for i in pair:
            help_list.append(i)
    root = bonds[number_atom][0]
    all_change = [number_atom]
    if root in help_list:
        for pair in sim_atoms:
            if root in pair:
                sim_pair = pair
                sim_atom_list = sim_change(sim_pair, C_atoms, bonds, number_atom, root, int_atoms, double_bonds, atoms_coord)
                for i in sim_atom_list:
                    all_change.append(i)
    if root not in help_list:
        all_change = [neig for neig in bonds[root] if neig not in C_atoms]
        if len(all_change) < 3:
            all_change = [number_atom]
    return all_change


def all_ident(ligand_name, flag='only sim_atoms'):
    """ Outputs a list with lists of all identical atoms in ligand's sceleton"""
    filename = ligand_name + '.sdf'
    download_sdf(ligand_name)
    C_atoms, double_bonds, graph, int_atoms, atoms_coord, bonds, chir_list = read_file(
        filename)  # прочитали файл, создали словарь с атомами и информацией о них, список двойных связей
    for i in C_atoms:  # убрали хиральность, связанную с наличием фтора
        if i in chir_list:
            H = 0
            F = 0
            for n in bonds[i]:
                if n not in C_atoms and n not in int_atoms:
                    H += 1
                if n in int_atoms:
                    F += 1
            if H == 1 and F == 1:
                C_atoms[i][3] -= chir_list[i]

    for i in range(len(C_atoms)):  # сам алгоритм ч1
        C_copy = copy.deepcopy(C_atoms)
        for p in C_copy:
            neig = neighbors(C_copy[p])
            for n in neig:
                C_copy[p][3] += C_atoms[n][3]
        C_atoms = C_copy

    sim_atom = ident_atoms(C_atoms)  # список с парами симметричных атомов

    for atom in C_atoms:
        C_atoms[atom][3] += int((ident_stereo(C_atoms, atom, double_bonds, bonds, int_atoms, atoms_coord) == 'R'))

    for k in graph.keys():
        vis = find_cycle(graph, k, [])

    non_cycle_double_bonds(cycle_atoms2, double_bonds, C_atoms, int_atoms, bonds, atoms_coord,
                           C_atoms)  # разделили по геометрии двойные связи вне циклов

    for pair in sim_atom:  # проверили потенциально симметричные атомы
        chain_double_bonds(C_atoms, pair, double_bonds, int_atoms, bonds, atoms_coord)
        one_root(C_atoms, pair, C_atoms)

    aromacycles = aromatic_cycle(double_bonds, cycle_atoms2, C_atoms,
                                 graph)  # обнулили все атомы в ароматических циклах
    for cycle in aromacycles:
        for atom in cycle:
            C_atoms[atom][3] = 0

    for i in range(len(C_atoms) + 3):  # сам алгоритм ч2
        if i % 3 == 0:
            sim_atom = ident_atoms(C_atoms)
            add_atoms = []
            for pair in sim_atom:
                if len(pair) == 2:
                    not_suit = joint_double_bond(double_bonds, C_atoms, pair)
                    for p in not_suit:
                        add_atoms.append(p)
            for pair in add_atoms:
                C_atoms[pair[0]][3] += 1
                C_atoms[pair[1]][3] += 2
        C_copy = copy.deepcopy(C_atoms)
        for p in C_copy:
            neig = neighbors(C_copy[p])
            for n in neig:
                C_copy[p][3] += C_atoms[n][3]
        C_atoms = C_copy

    sim_atoms = []  # окончательно определили симметричные атомы
    for i in ident_atoms(C_atoms):
        if i not in sim_atoms:
            sim_atoms.append(i)
    os.remove(filename)
    if flag =='only sim_atoms':
        return sim_atoms
    if flag =='all lists':
        return sim_atoms, C_atoms, bonds, int_atoms, double_bonds, atoms_coord


def ident_place(ligand_name, number_atom): #выводит номера атомов, симметричных данному(только водород/фтор)
    """ Outputs the numbers of atoms identical to the one entered. """
    filename = ligand_name + '.sdf'
    download_sdf(ligand_name)
    sim_atoms, C_atoms, bonds, int_atoms, double_bonds, atoms_coord = all_ident(ligand_name, flag='all lists')

    atoms_list_0 = change_2(sim_atoms, C_atoms, bonds, int_atoms, number_atom, double_bonds, atoms_coord)
    atoms_list = []
    for i in atoms_list_0:
        atoms_list.append(i) if i not in atoms_list else None
    os.remove(filename)
    return atoms_list


cycle_atoms2 = []
print(ident_place('59Q', 35))
print(all_ident('59Q'))
