import argparse
import os
import re
from tqdm import tqdm
import logging

import pmx
from pmx.AbsoluteDG import *
from pmx.gmx import get_gmx

import MDAnalysis as mda

from restraints import search
from restraints.restraints import FindBoreschRestraint
from AbsoluteDG_child import AbsoluteDG_child

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
os.environ['GMXLIB']='/home/jovyan/pmx/src/pmx/data/mutff'


# To add doc string and typing
def find_restraint_atoms(fe: AbsoluteDG_child, eq_tpr: str, eq_xtc: str, ligand_name: str) -> None:
    """_summary_

    Args:
        fe (AbsoluteDG): _description_
        eq_tpr (str): _description_
        eq_xtc (str): _description_
        ligand_name (str): The name of the ligand in the input PDB file
    """
    u = mda.Universe(eq_tpr, eq_xtc,  refresh_offsets=True)
    l_sel = f'resname {ligand_name} and not name H*'
    ligand_atoms = search.find_ligand_atoms(u, l_selection=l_sel, p_align="protein and name CA")
    atom_set = []
    for l_atoms in ligand_atoms:
        psearch = search.FindHostAtoms(u, l_atoms[0], p_selection="protein and name CA")
        psearch.run()
        atom_set.extend([(l_atoms, p) for p in psearch.host_atoms])
    boresch = FindBoreschRestraint(u, atom_set)
    boresch.run()
    l1, p1 = boresch.restraint.bond.atomgroup.ids
    l2 = boresch.restraint.angles[0].atomgroup.ids[0]
    p2 = boresch.restraint.angles[1].atomgroup.ids[-1]
    dih1 = boresch.restraint.dihedrals[0].atomgroup.ids
    dih3 = boresch.restraint.dihedrals[2].atomgroup.ids
    set1 = {l1, l2, p1, p2}
    for i in dih1:
        if i not in set1:
            l3 = i
    set2 = {l1, l2, p1, p2, l3}
    for i in dih3:
        if i not in set2:
            p3 = i
    p1 -= len(u.select_atoms(f'resname {ligand_name}'))
    p2 -= len(u.select_atoms(f'resname {ligand_name}'))
    p3 -= len(u.select_atoms(f'resname {ligand_name}'))
    fe.restraint_ligatoms = [l1, l2, l3]
    fe.restraint_recatoms = [p1, p2, p3]


def main():
    parser = argparse.ArgumentParser(description='Script for preparing and running ABFE calculations.')
    parser.add_argument('--file_prefix', type=str, required=True, help='Prefix of all the files in the working directory (e.g., 3HTB_clean).')
    parser.add_argument('--ligand_name', type=str, default='LIG', help='The name of the ligand in the input PDB file. Defaults to "LIG".')
    parser.add_argument('--processes', type=int, default=32, help='The number of processes to use for the calculations. Defaults to 16.')
    parser.add_argument('--lig', type=str, default='ligand', help='Ligand for the AbsoluteDG function. Defaults to "ligand".')
    parser.add_argument('--apoCase', type=str, default='apo', help='ApoCase for the AbsoluteDG function. Defaults to "apo".')
    parser.add_argument('--replicas', type=int, default=1, help='The number of replicas to use for the calculations. Defaults to 1.')
    parser.add_argument('--workdir', type=str, required=True, help='The working directory path.')
    args = parser.parse_args()
    
    fe = AbsoluteDG_child(ligList=[args.lig], apoCase=args.apoCase, bDSSB=True)
    fe.kBond = 2092.0
    fe.kAngle = 20.92
    fe.kDihedral = 125.52
    fe.frameNum = 80
    fe.bLigMidSpread = False
    fe.bLigMaxSpread = False
    fe.rdssb = 1.0
    fe.boxd = 0.6
    fe.workPath = args.workdir
    fe.mdpPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'alchemistry_mdps')
    fe.replicas = args.replicas
    fe.structTopPath = os.path.join(args.workdir, 'struct_top') # No need to change it

    find_restraint_atoms(fe, 
                         os.path.join(os.path.join(args.workdir, 'eq_files'), f'{args.file_prefix}_eq.tpr'), 
                         os.path.join(os.path.join(args.workdir, 'eq_files'), f'{args.file_prefix}_eq.xtc'),
                         args.ligand_name)

    fe.simTypes = ['em','eq_posre','eq','transitions']
    fe.prepareFreeEnergyDir()
    # assemble the systems 
    fe.assemble_systems()

    # box/water/ions
    fe.boxWaterIons()
    # em_posre
    fe.prepare_simulation(simType='em')
    fe.run_simulation(simType='em', nt=args.processes)

    # Do the same for all other simTypes consecutively. When running transitions, make sxure to set fe.equilTime, fe.bGenTiTpr
    # eq_posre
    fe.prepare_simulation( simType='eq_posre', prevSim='em' )
    fe.run_simulation(simType='eq_posre', nt=args.processes)

    # eq
    fe.prepare_simulation( simType='eq', prevSim='eq_posre' )
    fe.run_simulation(simType='eq', nt=args.processes)

    fe.equilTime = 1080.0 # ps to discard as equilibration
    fe.bGenTiTpr = True
    
    # !!!!!!!!!!!!!! We need to adapt the paths through the arguments !!!!!!!!!!!!
    tprpaths = fe.prepare_simulation( simType='transitions')
    gmxexec = get_gmx()

    for path in tprpaths:
        parentdir = os.path.dirname(os.path.realpath(path))
        basename = os.path.basename(os.path.dirname(os.path.realpath(path)))
        pattern = r'\d+'
        subprocess.run(f'{gmxexec} mdrun -nb gpu -s {parentdir}/tpr.tpr -dhdl {parentdir}/dhdl{re.findall(pattern, basename)[0]}.xvg -e {parentdir}/ener.edr -g {parentdir}/md.log -cpo {parentdir}/state.cpt -nt {args.processes}', shell=True, check=True, capture_output=False)
    

    correction_dir = f"{args.workdir}/../.."
    os.system(f"python {correction_dir}/StandardState.py > {correction_dir}/energy_correction.log")
    fe.run_analysis( ligs=[args.lig] )

    
if __name__ == "__main__":
    # chown current directory to user
    # os.system("sudo chown -R 1000:100 .")
    main()
