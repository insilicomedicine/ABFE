about = """Removes misplaced water molecules"""

import argparse
from shutil import copyfile

import networkx as nx
import numpy as np
# Use MDAnalysis instead of pmx
from MDAnalysis import Universe
from MDAnalysis.analysis.distances import distance_array


def find_indices(l, m):
    return np.where(np.in1d(l, m))[0]


def valid_bb_hb(atom1, atom2):
    n1, n2 = atom1.type, atom2.type
    i1, i2 = atom1.resnum, atom2.resnum
    if abs(i1 - i2) >= 3:
        return n1 == 'N' and n2 == 'O' or n1 == 'O' and n2 == 'N'
    else:
        return False


def valid_bbsc_hb(atom1, atom2):
    n1, n2 = atom1.type, atom2.type
    return n1 in ['N', 'O'] and n2 not in ['N', 'O']


def valid_hb(atom1, atom2):
    term1 = valid_bb_hb(atom1, atom2)
    term2 = valid_bbsc_hb(atom1, atom2)
    term3 = valid_bbsc_hb(atom2, atom1)

    return term1 or term2 or term3


def score_waters(gro_path, dw):
    """Read shell PDB and calculate number of h-bonds for each water"""
    gro = Universe(gro_path)
    "Select all water molecules 7.5A around protein"
    shell = gro.select_atoms('(around 7.5 protein) and name OW')
    polar = gro.select_atoms('(type N or type O) and not resname SOL')
    shell = shell
    # Added to take into account water on the inner outer shell edge
    inner_shell = gro.select_atoms('(around 5 not resname SOL) and name OW')
    I_shell = inner_shell.atoms.ix
    i_to_ri = dict(zip(shell.atoms.ix, list(shell.resindices)))

    D_polar = distance_array(polar.positions, polar.positions, box=gro.dimensions)
    global I_polar
    I_polar = polar.atoms.ix_array
    G_polar = nx.Graph()

    for i, x in enumerate(polar):
        G_polar.add_node(i)
        for j in range(i + 1, len(polar)):
            y = polar[j]
            if D_polar[i, j] <= 3.5 and valid_hb(x, y):
                G_polar.add_edge(i, j)

    not_hbonding = {
        I_polar[i]
        for i, d in G_polar.degree
        if polar[i].type in ['N', 'O'] and d >= 1
    }
    global nhb
    nhb = sorted(not_hbonding)
    # Advanced indexing
    assel = find_indices(np.array(I_polar), np.array(list(nhb)))
    mask = np.ones(len(polar), dtype=bool)
    mask[assel[0]] = False
    hbonding = polar[mask]
    system = hbonding + shell
    D_sys = distance_array(system.positions, system.positions, box=gro.dimensions)
    I_sys = system.atoms.ix
    G_sys = nx.Graph()

    for i, x in enumerate(system):
        G_sys.add_node(i)
        for j in range(i + 1, len(system)):
            if 2.5 <= D_sys[i, j] <= dw:
                G_sys.add_edge(i, j)

    scores = []
    for i, n_hb in G_sys.degree:
        idx = I_sys[i]
        if idx in I_shell:
            scores.append((i_to_ri[idx], n_hb))

    return scores


def renumber_waters(waters_ri, renumber_dict):
    return [renumber_dict[x] for x in waters_ri]


def remove_waters_gro(wet_gro, new_gro, to_remove):
    removed = []
    lines = open(wet_gro, 'r').readlines()
    n_atoms = 3 * len(to_remove)

    with open(new_gro, 'w') as fh:
        fh.write(lines[0])
        N = int(lines[1].strip())
        N -= n_atoms
        fh.write(' %d\n' % N)

        for line in lines[2:-1]:
            resnr = int(line[:5])
            if resnr in to_remove:
                removed.append(resnr)
                continue
            else:
                fh.write(line)

        fh.write(lines[-1])
    print('Not removed')
    print(set(to_remove) - set(removed))

    return len(to_remove)


def remove_waters_top(top, to_remove):
    copyfile(top, f'{top}.bkp')
    N_sol_remove = len(to_remove)

    lines = open(f'{top}.bkp', 'r').readlines()
    with open(top, 'w') as out:
        N_sol = _remove(lines, out, N_sol_remove)
    return N_sol - N_sol_remove


def _remove(lines, out, N_sol_remove):
    in_molecules = False
    sol_block = False
    after_sol = False
    result = 0
    for line in lines:
        if line.startswith('[ molecules ]'):
            in_molecules = True
        if not in_molecules:
            out.write(line)
        elif line.startswith('SOL'):
            sol_block = True
            tk = line.strip().split()
            result += int(tk[1])
        else:
            if sol_block:
                after_sol = True
                out.write('SOL %d\n' % (result - N_sol_remove))
                sol_block = False
            out.write(line)
    if not after_sol:
        out.write('SOL %d\n' % (result - N_sol_remove))
    return result


def main():
    parser = argparse.ArgumentParser(description=about)
    parser.add_argument('wet',
                        metavar='wet_gro',
                        type=str,
                        help='''Path to fully solvated gro. \
                                End-to-end consecutive numbering is required.''')
    parser.add_argument('out',
                        metavar='new_gro',
                        type=str,
                        help='Path to output gro')
    parser.add_argument('top',
                        metavar='topology',
                        type=str,
                        help='Path to topology')
    parser.add_argument('-d',
                        type=float,
                        help='''Max length of hbond to water O. \
                                Default is 3.2.''',
                        default=3.2)
    parser.add_argument('-n',
                        type=int,
                        help='''Min number of hbonds to water to keep it. \
                                Default is 2''',
                        default=2)

    args = parser.parse_args()

    water_scores = score_waters(args.wet, args.d)
    remove = [x[0] for x in water_scores if x[1] < args.n]
    N_removed = remove_waters_gro(args.wet, args.out, remove)
    N_present = remove_waters_top(args.top, remove)
    print("Removed %d waters, now %d waters." % (N_removed, N_present))


if __name__ == '__main__':
    main()