import pmx
from pmx import *
#from pmx.utils import create_folder
from pmx import gmx, ligand_alchemy, jobscript
from pmx.parser import *
from pmx.forcefield import *
from pmx.utils import * #create_folder, clean_gromacs_backup_files
from pmx.DoubleBox import DoubleBox
#from pmx.forcefield import TopolBase
from pmx.geometry import *
import sys
import os,shutil
import re
import subprocess
import glob
import random
import numpy as np
import pandas as pd
import copy as cp
from pmx.Restraints import *
import argparse
from tqdm import tqdm
import logging

from pmx.AbsoluteDG import *
from pmx.gmx import get_gmx

import MDAnalysis as mda

from restraints import search
from restraints.restraints import FindBoreschRestraint


import MDAnalysis as mda
from MDAnalysis.analysis.distances import dist
from MDAnalysis.lib.distances import calc_dihedrals
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt

class AbsRestraints_pairwise:

    def __init__(self, topology_path, traj_path, outII_path):

        self.sysTraj = traj_path
        self.sysTop = topology_path

        self._findAnomPairs()
        self.outii = outII_path

        self._write_ii()


    def _findAnomPairs(self):

        u = mda.Universe(self.sysTop, self.sysTraj, topology_format='ITP', include_dir="/usr/local/share/gromacs/top/")
        lig_heavy = u.select_atoms("resname LIG and not name H*")

        anchors_dict = {}
        for lig_atom in lig_heavy:
            for prot_atom in u.select_atoms(f"(protein or resname PRT) and (around 8 index {lig_atom.index}) and (not name H*) and backbone"): # protein does not recognise PRT
                anchors_dict[(lig_atom.index,prot_atom.index)]={}
                anchors_dict[(lig_atom.index, prot_atom.index)]["dists"]=[]


        for frame in u.trajectory:
            for lig_atom_index, prot_atom_index in anchors_dict.keys():
                distance = dist(mda.AtomGroup([u.atoms[lig_atom_index]]), mda.AtomGroup([u.atoms[prot_atom_index]]), box=frame.dimensions)[2][0]
                anchors_dict[(lig_atom_index,prot_atom_index)]["dists"].append(distance)


        # change lists to numpy arrays
        for pair in anchors_dict.keys():
            anchors_dict[pair]["dists"] = np.array(anchors_dict[pair]["dists"])

        # calculate average and SD
        for pair in anchors_dict.keys():
            anchors_dict[pair]["avg_dist"] = anchors_dict[pair]["dists"].mean()
            anchors_dict[pair]["sd_dist"] = anchors_dict[pair]["dists"].std()
        
# get n pairs with lowest SD
        pairs_ordered_sd=[]
        for item in sorted(anchors_dict.items(), key=lambda item: item[1]["sd_dist"]):
            pairs_ordered_sd.append(item[0])
            #print(f'Pair: {item[0]}, av distance: {item[1]["avg_dist"]:.2f}, SD: {item[1]["sd_dist"]:.2f}')
            

        dists_precor = []
        atom_pairs_dist = []
        #Prepare pairs for correlations search

        for i in range(80):
            atom1 = u.atoms[pairs_ordered_sd[i][0]]
            atom2 = u.atoms[pairs_ordered_sd[i][1]]
            atom_pairs_dist.append((atom1.index, atom2.index))
            dists_precor.append(anchors_dict[(atom1.index, atom2.index)]["dists"])

        dist_array = np.array(dists_precor)

        #Correlations search

        correlations = []

        for i in range(len(dist_array)):
            for j in range(i+1, len(dist_array)):
                corr_value = np.corrcoef(dist_array[i], dist_array[j])[0, 1]
                correlations.append(((atom_pairs_dist[i], atom_pairs_dist[j]), corr_value))

        correlations.sort(key=lambda x: abs(x[1]))

        pairs_ordered_corr = []

        #Choose part of correlated pairs

        for i in correlations:
            if abs(i[1]) < 0.5:
                pairs_ordered_corr.append(i[0][0])
                pairs_ordered_corr.append(i[0][1])

        unique_pairs_ordered_sd = []
        lig_ats = []
        recept_ats = []

        #Choose only pairs with unique atoms

        for pair in pairs_ordered_corr:
            if pair[0] not in lig_ats:
                if pair[1] not in recept_ats:
                    unique_pairs_ordered_sd.append(pair)
                    lig_ats.append(pair[0])
                    recept_ats.append(pair[1])
                    if len(unique_pairs_ordered_sd) == 22:
                        break
                    
        self.restraints = {}

        #Generate restraints

        for pair in unique_pairs_ordered_sd:

            r0 = round(anchors_dict[pair]["avg_dist"], 2)

            kl = 2.478 / ((anchors_dict[pair]["sd_dist"]*0.1)**2)
            kl = round(kl, 2)

            if kl > 8000: kl = 8000

            dists = np.array(anchors_dict[pair]['dists'])
            lower_bound = np.percentile(dists, 2.5)
            upper_bound = np.percentile(dists, 97.5)


            dl = 0
            
            new_pair = (u.atoms[pair[0]].id, u.atoms[pair[1]].id)
            self.restraints[new_pair] = (r0, kl, dl, lower_bound, upper_bound)



    def _write_ii(self):

        with open(self.outii, "w") as file:
            file.write("[ intermolecular_interactions ]\n")
            file.write("[ bonds ]\n")
            for key, value in self.restraints.items():
                atom1, atom2 = key
                low = value[3] / 10
                up1 = value[4] / 10
                k_dr = value[1]
                line = f"{atom1:<5d}{atom2:<5d}10 {low:<11.6f}{up1:<11.6f}{1000:<11.0f}{k_dr}\n"
                file.write(line)


class AbsoluteDG_child(AbsoluteDG):
    def _assemble_system_topologies( self, lig='', case='holo', molList=[] ):
        
        strPath = self._get_specific_path( lig=lig, case=case, bStrTop=True )
        
        if self.bDSSB==True:
            outPath = self._get_specific_path( lig=lig, state=case, wp='dssb' )
            outPathHolo = self._get_specific_path( lig=lig, state='holo', wp='dssb' )
        else:
            outPath = self._get_specific_path( lig=lig, state=case, wp='protein' )                    
            outPathHolo = self._get_specific_path( lig=lig, state='holo', wp='protein' )                    
            

        ##### HYBRID TOPOLOGIES (in both cases DSSB/noDSSB)
        # create hybrid topologies
        ffitpfile = self.ffitp[lig][case]
#        if self.bDSSB==True:
        itpfile = '{0}/{1}'.format(strPath,self.ligands[lig]['ligand'].topItp)
        itpfileA = '{0}A.itp'.format(self._name_without_extension(itpfile))
        itpfileB = '{0}B.itp'.format(self._name_without_extension(itpfile))
        ffitpfile = self.ffitp[lig][case]
        topDecoupleObj = TOPdecouple( ff=self.ff, itpfile=itpfile, itpfileA=itpfileA, itpfileB=itpfileB, ffitpfile=ffitpfile )
        newMolNames = topDecoupleObj.newMolNames # for DSSB the molecules are renamed in the ITP
        newMolItps = [itpfileA,itpfileB]    
        
        
        ####### CREATE TOP ######
        topFname = '{0}/topol.top'.format(outPath)
        itps = [ os.path.relpath( self.ffitp[lig][case], start=outPath ) ]
        mols = []
        systemName = 'system {0}'.format(case)
        water = None # water is special, i.e. append to the end
        
        for mol in molList:
            # water comes last, here only save it
            if mol.molItpName=='SOL' or mol.molItpName=='HOH' or mol.molItpName=='WAT' or \
               mol.molItpName=='Water' or mol.molItpName=='TIP3':
                water = mol
            else:
                molName = mol.molItpName
                # in case of DSSB need to replace molecule names, as there are now two different molecules
                if molName==self.ligands[lig]['ligand'].molItpName: # and self.bDSSB==True:
                    molName = newMolNames.pop(0)
                    molItp = os.path.relpath( newMolItps.pop(0),start=outPath )
                    itps.append( molItp )
                else:
                    if mol.topItp!=None:
                        molItp = os.path.relpath( '{0}/{1}'.format(strPath,mol.topItp),start=outPath )                   
                        itps.append( molItp )                    
                mols.append([molName,mol.molNumber])   

        # water comes last
        if water!=None:
            mols.append([water.molItpName,water.molNumber])
        
        self._create_top(fname=topFname,itp=itps,mols=mols,systemName=systemName)         
        
        ####### INTERMOLECULAR RESTRAINTS ######
        # for holo state generate the restraints
        self.iiRestrFile[case] = None
        if case=='holo':
            strFile = '{0}/system.pdb'.format(outPath)
            ligItpFile = self.ligands[lig]['ligand'].topItpPath          
            protItpFile = self.ligands[lig]['protein'].topItpPath
            indLig = np.arange(0,self.ligands[lig]['ligand'].natoms)
            indProt = np.arange(self.ligands[lig]['ligand'].natoms,self.ligands[lig]['protein'].natoms)
            self.iiRestrFile[case] = '{0}/ii.itp'.format(strPath)
            restrDgFile = '{0}/restr_dG.dat'.format(strPath)
            restrAtomsFile = '{0}/restr_atoms.pdb'.format(strPath)
            # may require apo file for restraint generation
            strFileApo = None
            if self.apoCase!=None: 
                outPathApo = self._get_specific_path( lig=lig, state='apo', wp='protein' ) 
                strFileApo = '{0}/system.pdb'.format(outPathApo)
            

            top_path_restr = self.workPath + "/complex.top"
            traj_path_restr = self.workPath + "/traj.xtc"
            AbsRestraints_pairwise(topology_path = top_path_restr, traj_path = traj_path_restr, outII_path = self.iiRestrFile[case])
        
        
        else: # for apo state, copy or transfer the restraints
            strPathHolo = self._get_specific_path( lig=lig, case='holo', bStrTop=True )
            iiRestrFileHolo = '{0}/ii.itp'.format(strPathHolo)
            self.iiRestrFile[case] = '{0}/ii.itp'.format(strPath)

            cmd = 'cp {0} {1}'.format(iiRestrFileHolo,self.iiRestrFile[case])
            os.system(cmd)            
            
        # save the ii file to later append it to the final top file
        self.iiRestrFile[case] = os.path.relpath( self.iiRestrFile[case], start=outPath )        
        
        ####### IF DSSB, RESTRAIN TWO COMs: PROTEIN AND LIGAND ######
        if self.bDSSB==True:        
            # restrain ligand (topology B)
            ligandStrFile = self.ligands[lig]['ligand'].structPdbPath
            ligandItpFile = itpfileB
            outLigandItpFile = itpfileB                    
            PositionRestraints( itpname=ligandItpFile, outname=outLigandItpFile, 
                                pdbfile=ligandStrFile, bPosreIfdef=False, stateBonded='A',
                                bIgnoreIfdefs=True)
            
            # restrain protein
            protStrFile = self.ligands[lig]['protein'].structPdbPath
            protItpFile = self.ligands[lig]['protein'].topItpPath
            if case=='apo' and self.apoCase!=None:
                protStrFile = self.ligands[self.apoCase]['protein'].structPdbPath
                protItpFile = self.ligands[self.apoCase]['protein'].topItpPath                
            outProtItpFile = '{0}/{1}'.format(strPath,self.ligands[lig]['protein'].topItp) 
            PositionRestraints( itpname=protItpFile, outname=outProtItpFile, 
                                pdbfile=protStrFile, name='CA', bPosreIfdef=False)


