import sys, os
import copy as cp
import numpy as np
from pmx import *
from pmx import geometry
from pmx import ndx
from pmx import library
from pmx.forcefield import *
import scipy
from scipy.special import erf
from string import digits
from pmx import library as pmxlib

class AbsRestraints:

    def __init__( self, **kwargs ):

        # constants
        R = 8.31445985*0.001  # Gas constant in kJ/mol/K
        self.T = 298.15
        self.RT = R*self.T

        # structure: one of those need to be present
        self.system = None # pmx model
        self.strFile = None 
        self.strFileApo = None

        # booleans
        self.bLigMaxSpread = False
        self.bLigMidSpread = True
        self.bMartini = False
        self.bManLig = False
        self.bManProt = False
        self.bVerbose = False

        # force constants
        self.kBond = 4184.0
        self.kAngle = 41.84
        self.kDihedral = 41.84

        # optional max distance allowed between apo and holo atoms (in nm)
        self.apoDistanceLimit = 0.2 # nm

        # files
        self.topLigFile = None
        self.topProtFile = None
        self.ndxFile = None
        self.dgFile = None
        self.outiiFile = None
        self.pdbRestrAtoms = None # an optional output with the restrained atoms only (convenient for viewing)

        # ligandID: may provide manually
        self.indLig1 = False
        self.indLig2 = False
        self.indLig3 = False

        # proteinID: may provide manually
        self.indProt1 = False
        self.indProt2 = False
        self.indProt3 = False

        # center ligand name
        self.ligCenterByName = False

        # indeces
        self.indLig = []
        self.indProt = []

        for key, val in kwargs.items():
            setattr(self,key,val)

        # pdb file is necessary
        if self.strFile==None:
            sys.stdout.write('Need to provide structure file')
            sys.exit(0)
        self.system = Model(self.strFile) 

        # if apo pdb is given
        self.systemApo = None
        if self.strFileApo!=None:
            self.systemApo = Model(self.strFileApo) 

        # topologies
        self.protTop = None
        self.ligTop = None
        if len(self.topLigFile)!=0:
            self.ligTop = TopolBase( self.topLigFile )
        if len(self.topProtFile)!=0:
            self.protTop = TopolBase( self.topProtFile )

        # merge topologies
        # if self.ligTop != None and self.protTop != None:
#            mergedTop = self._merge_top( )
#            renumLigTop = self._renum_lig_top( )

        if len(self.indProt)==0:
            self.indProt = self._select_ndx(self.ndxFile,message='Select index group for protein:\n')
        if len(self.indLig)==0:
            self.indLig = self._select_ndx(self.ndxFile,message='Select index group for ligand:\n')

        # ligand indeces
        if self.bManLig==True:
            if self.indLig1==False:
                self.indLig1 = self._select_ndx(self.ndxFile,message='Select first ligand atom for distance restraint:\n')
                self._check_ndx_size( len(self.indLig1), desiredSize=1 )
            if self.indLig2==False:
                self.indLig2 = self._select_ndx(self.ndxFile,message='Select second ligand atom for angle restraint:\n')
                self._check_ndx_size( len(self.indLig2), desiredSize=1 )
            if self.indLig3==False:
                self.indLig3 = self._select_ndx(self.ndxFile,message='Select third ligand atom for dihedral restraint:\n')
                self._check_ndx_size( len(self.indLig3), desiredSize=1 )
        elif self.ligCenterByName!=False: # central ligand atom provided manually, the rest of atoms automatically
            self.indLig1 = self._ligAtomByName( self.ligCenterByName, self.system )
            if self.indLig1==False:
                sys.stdout.write('Cannot find the ligand atom ',self.ligCenterByName)
                sys.stdout.write('Exiting...')
                sys.exit(1)
#            indLig1 = self._select_ndx(ndxFile,message='Select one ligand atom for distance restraint:\n')
#            self._check_ndx_size( len(indLig1), desiredSize=1 )
            (self.indLig1,self.indLig2,self.indLig3) = self._identify_ligand_atoms( self.system, self.indLig1 )
        else: # find all ligand atoms automatically
            (self.indLig1,self.indLig2,self.indLig3) = self._identify_ligand_atoms( self.system, self.indLig1 )

        # atoms
        self.aLig1 = self.system.atoms[self.indLig1[0]]
        self.aLig2 = self.system.atoms[self.indLig2[0]]
        self.aLig3 = self.system.atoms[self.indLig3[0]]
        self.aLig1.a2nm()
        self.aLig2.a2nm()
        self.aLig3.a2nm()

        # check if all went fine
        if self.indLig1[0]==-42 or self.indLig2[0]==-42 or self.indLig3[0]==-42:
            print('Something went wrong in ligand atom selection. Exiting...')
            sys.exit(0)

        # protein indeces
        if self.bManProt==True:
            if self.indProt1==False:
                self.indProt1 = self._select_ndx(self.ndxFile,message='Select first protein atom for distance restraint:\n')
                self._check_ndx_size( len(self.indProt1), desiredSize=1 )
            if self.indProt2==False:
                self.indProt2 = self._select_ndx(self.ndxFile,message='Select second protein atom for angle restraint:\n')
                self._check_ndx_size( len(self.indProt2), desiredSize=1 )
            if self.indProt3==False:
                self.indProt3 = self._select_ndx(self.ndxFile,message='Select third protein atom for dihedral restraint:\n')
                self._check_ndx_size( len(self.indProt3), desiredSize=1 )
        else: # find protein indeces
            bRandProtAtoms = False
            if bRandProtAtoms==True:
                (self.indProt1,self.indProt2,self.indProt3) = self._random_protein_atoms( self.system )
            else:
                (self.indProt1,self.indProt2,self.indProt3) = self._identify_protein_atoms( self.system, self.indLig1, self.indLig2, self.indLig3 )
        # atoms
        self.aProt1 = self.system.atoms[self.indProt1[0]]
        self.aProt2 = self.system.atoms[self.indProt2[0]]
        self.aProt3 = self.system.atoms[self.indProt3[0]]
        self.aProt1.a2nm()
        self.aProt2.a2nm()
        self.aProt3.a2nm()

        # check if all went fine
        if self.indProt1[0]==-42 or self.indProt2[0]==-42 or self.indProt3[0]==-42:
            print('Something went wrong in protein atom selection.')
            print('Atoms identified: {0} {1} {2}'.format(self.indProt1[0],self.indProt2[0],self.indProt3[0]))
            print('Exiting...')
            sys.exit(0)

        # distance
        self.dist = self.aProt1 - self.aLig1

        # angles
        # aLig2 - aLig1 - aProt1
        self.ang1 = self.aLig1.angle(self.aLig2,self.aProt1,degree=True)
        self._check_angle( self.ang1 )
        # aLig1 - aProt1 - aProt2
        self.ang2 = self.aProt1.angle(self.aLig1,self.aProt2,degree=True)
        self._check_angle( self.ang2 )

        # dihedrals
        # aLig3 - aLig2 - aLig1 - aProt1
        self.dih1 = self.aLig3.dihedral(self.aLig2,self.aLig1,self.aProt1,degree=True)
        # aLig2 - aLig1 - aProt1 - aProt2
        self.dih2 = self.aLig2.dihedral(self.aLig1,self.aProt1,self.aProt2,degree=True)
        # aLig1 - aProt1 - aProt2 - aProt3
        self.dih3 = self.aLig1.dihedral(self.aProt1,self.aProt2,self.aProt3,degree=True)


        # generate [ intermolecular_interactions ]
#        mergedTop.has_ii = True
        self.ii = {}
        self._gen_ii_bonds()
        self._gen_ii_angles()
        self._gen_ii_dihedrals()

        # optional output of restrained atoms as .pdb
        if self.pdbRestrAtoms!=None:
            mout = Model()
            mout.atoms.extend( (self.aLig1,self.aLig2,self.aLig3,self.aProt1,self.aProt2,self.aProt3) )
            mout.write(self.pdbRestrAtoms)

        # write intermolecular interaction
        self._write_ii( self.outiiFile )

        # calculate restraint contribution to the free energy
        if self.dgFile!=None: 
            fp = open(self.dgFile,'w')
        self._calc_restraint_dg( )
        if self.bVerbose==True:
            sys.stdout.write('Restraint contribution to free energy (according to B-K): %3.4f kJ/mol\n' % self.dg )
            sys.stdout.write('Restraint contribution to free energy (according to B-K): %3.4f kcal/mol\n' % self.dgkcal )
            sys.stdout.write('\n')
        if self.dgFile!=None: 
            fp.write('Restraint contribution to free energy (according to B-K): %3.4f kJ/mol\n' % self.dg )
            fp.write('Restraint contribution to free energy (according to B-K): %3.4f kcal/mol\n' % self.dgkcal )
            fp.write('\n')

        self._calc_restraint_dg_w_limits( )
        if self.bVerbose==True:
            sys.stdout.write('Restraint contribution to free energy (w ideal limits): %3.4f kJ/mol\n' % self.dg )
            sys.stdout.write('Restraint contribution to free energy (w ideal limits): %3.4f kcal/mol\n' % self.dgkcal )
            sys.stdout.write('\n')
        if self.dgFile!=None: 
            fp.write('Restraint contribution to free energy (w ideal limits): %3.4f kJ/mol\n' % self.dg )
            fp.write('Restraint contribution to free energy (w ideal limits): %3.4f kcal/mol\n' % self.dgkcal )
            fp.write('\n')
        
        self._calc_restraint_dg_w_gromacs_limits( )
        if self.bVerbose==True:
            sys.stdout.write('Restraint contribution to free energy (w gmx limits): %3.4f kJ/mol\n' % self.dg )
            sys.stdout.write('Restraint contribution to free energy (w gmx limits): %3.4f kcal/mol\n' % self.dgkcal )
        if self.dgFile!=None: 
            fp.write('Restraint contribution to free energy (w gmx limits): %3.4f kJ/mol\n' % self.dg )
            fp.write('Restraint contribution to free energy (w gmx limits): %3.4f kcal/mol\n' % self.dgkcal )
            fp.close()

        # molecule name
#        mergedTop.name = 'merged'

        # write topology
 #       if outItpFile!=False:
  #          renumLigTop.write(outItpFile,stateBonded='A')
        #    mergedTop.write(outItpFile,stateBonded='A')

    def _check_angle(self, angle):
        extr1 = 0.5*self.kAngle*np.power((angle-0.0)/180.0*np.pi,2)
        extr2 = 0.5*self.kAngle*np.power((angle-180.0)/180.0*np.pi,2)
        if extr1/self.RT<5.0:
            if self.bVerbose==True:
                sys.stdout.write('Warning: one of the angle restraints is less than 5 kT away from 0 degree angle\n')
            return(False)
        elif extr2/self.RT<5.0:
            if self.bVerbose==True:
                sys.stdout.write('Warning: one of the angle restraints is less than 5 kT away from 180 degree angle\n')
            return(False)
        return(True)

    def _calc_restraint_dg(self):
        """Calculates effect of restraints in gas phase following https://doi.org/10.1021/jp0217839.
        
        This is the approximate result with no limits placed on the restraint coordinates.
        """
        
        V0 = 1.66            # standard volume in nm^3
        dgPrefactor = ( 8.0*np.power(np.pi,2.0)*V0/(np.power(self.dist,2.0)*np.sin(self.ang1*np.pi/180.0)*np.sin(self.ang2*np.pi/180.0)) )
        dgForceConstants = np.sqrt(self.kBond*self.kAngle*self.kAngle*self.kDihedral*self.kDihedral*self.kDihedral)/np.power(2.0*np.pi*self.RT,3.0)
        self.dg = -self.RT * np.log(dgPrefactor*dgForceConstants)
        self.dgkcal = self.dg/4.184
        
        
    def _calc_restraint_dg_w_limits(self):
        """Calculates effect of restraints in gas phase taking limits into account.
        
        This is the "ideal" case with limits symmetric around the equilibrium angles.
        dist       [0, inf]
        ang-ang_0  [-Pi, Pi]
        dih-dih_0  [-Pi, Pi]
        """
        
        V0 = 1.66            # standard volume in nm^3
        
        #-1^6 = 1;
        #Boresch, ... and Karplus put their Is in here
        # and they work out to I=-2 given infinite limits of constrant variables.
        # So in their case there is a sqrt(2*Pi*RT)^6 instead of sqrt(0.5*Pi*RT)^6
        integral_multiplier = np.power(np.pi*self.RT/2.0, 3.0)
        
        I_dist = -1 - erf( np.sqrt(self.kBond/(2.0*self.RT)) * self.dist )
        I_ang = erf( np.sqrt(self.kAngle/(2.0*self.RT)) * (-np.pi) ) - \
                erf( np.sqrt(self.kAngle/(2.0*self.RT)) * (np.pi) )
        I_dih = erf( np.sqrt(self.kDihedral/(2.0*self.RT)) * (-np.pi) ) - \
                erf( np.sqrt(self.kDihedral/(2.0*self.RT)) * (np.pi) )
        
        I_ang = [I_ang, I_ang]
        I_dih = [I_dih, I_dih, I_dih]
        
        forceConstants = np.sqrt(self.kBond*self.kAngle*self.kAngle*self.kDihedral*self.kDihedral*self.kDihedral)
        partition_func = np.power(self.dist,2.0)*np.sin(self.ang1*np.pi/180.0)*np.sin(self.ang2*np.pi/180.0) * \
                         integral_multiplier * I_dist * I_ang[0] * I_ang[1] * \
                         I_dih[0] * I_dih[1] * I_dih[2] / forceConstants
        
        self.dg = -self.RT * np.log(8.0*np.power(np.pi,2.0)*V0 / partition_func)
        self.dgkcal = self.dg/4.184
        
        
    def _calc_restraint_dg_w_gromacs_limits(self):
        """Calculates effect of restraints in gas phase using limits of gromacs harmonic restraints.
        
        This is the case for harmonic restraint potentials in gromacs.
        dist       [0, inf]
        ang        [0, Pi]
        dih-dih_0  [-Pi, Pi]
        """
        
        V0 = 1.66            # standard volume in nm^3
        
        #-1^6 = 1;
        #Boresch, ... and Karplus put their Is in here
        # and they work out to I=-2 given infinite limits of constrant variables.
        # So in their case there is a sqrt(2*Pi*RT)^6 instead of sqrt(0.5*Pi*RT)^6
        integral_multiplier = np.power(np.pi*self.RT/2.0, 3.0)
        
        I_dist = -1 - erf( np.sqrt(self.kBond/(2.0*self.RT)) * self.dist )
        I_angA = erf( np.sqrt(self.kAngle/(2.0*self.RT)) * ((self.ang1*np.pi/180.0)-np.pi) ) - \
                 erf( np.sqrt(self.kAngle/(2.0*self.RT)) * (self.ang1*np.pi/180.0) )
        I_angB = erf( np.sqrt(self.kAngle/(2.0*self.RT)) * ((self.ang2*np.pi/180.0)-np.pi) ) - \
                 erf( np.sqrt(self.kAngle/(2.0*self.RT)) * ((self.ang2*np.pi/180.0)) )
        I_dih = erf( np.sqrt(self.kDihedral/(2.0*self.RT)) * (-np.pi) ) - \
                erf( np.sqrt(self.kDihedral/(2.0*self.RT)) * (np.pi) )
                
        #print I_angA, I_angB
#        print erf( sqrt(self.kAngle/(2.0*self.RT)) * (self.ang1-np.pi) ), sqrt(self.kAngle/(2.0*self.RT)) * (self.ang1-np.pi), \
#              erf( sqrt(self.kAngle/(2.0*self.RT)) * (self.ang1) ), sqrt(self.kAngle/(2.0*self.RT)) * (self.ang1)
                
        I_ang = [I_angA, I_angB]
        I_dih = [I_dih, I_dih, I_dih]
        
        forceConstants = np.sqrt(self.kBond*self.kAngle*self.kAngle*self.kDihedral*self.kDihedral*self.kDihedral)
        partition_func = np.power(self.dist,2.0)*np.sin(self.ang1*np.pi/180.0)*np.sin(self.ang2*np.pi/180.0) * \
                         integral_multiplier * I_dist * I_ang[0] * I_ang[1] * \
                         I_dih[0] * I_dih[1] * I_dih[2] / forceConstants
        
        self.dg = -self.RT * np.log(8.0*np.power(np.pi,2.0)*V0 / partition_func)
        self.dgkcal = self.dg/4.184

    def _write_ii( self, outiiFile ):
        fp = open(outiiFile, 'w')
        fp.write('\n [ intermolecular_interactions ]\n')
        # bonds
        if 'bonds' in self.ii.keys():
           fp.write(' [ bonds ]\n')
           for b in self.ii['bonds']:
               fp.write('%6d %6d %6d' % ( b[0].id, b[1].id, b[2] ))
               if len(b)>3:
                   for x in b[3]:
                       fp.write(' %14.6f' % x)
               fp.write('\n')
           fp.write('\n')
        # angles
        if 'angles' in self.ii.keys():
           fp.write(' [ angles ]\n')
           for ang in self.ii['angles']:
               fp.write('%6d %6d %6d %6d' % ( ang[0].id, ang[1].id, ang[2].id, ang[3] ))
               if len(ang)>4:
                   for x in ang[4]:
                       fp.write(' %14.6f' % x)
               fp.write('\n')
           fp.write('\n')
        # dihedrals
        if 'dihedrals' in self.ii.keys():
           fp.write(' [ dihedral_restraints ]\n')
           for dih in self.ii['dihedrals']:
               fp.write('%6d %6d %6d %6d %6d' % ( dih[0].id, dih[1].id, dih[2].id, dih[3].id, dih[4] ))
               if len(dih)>5:
                   for x in dih[5]:
                       fp.write(' %14.6f' % x)
               fp.write('\n')
           fp.write('\n')
        fp.close()

    def _gen_ii_bonds( self ):
        #bond1 = [ self.aLig1, self.aProt1, 6, [self.dist, 0.0, self.dist, self.kBond] ]
        bond1 = [ self.aLig1, self.aProt1, 6, [self.dist, self.kBond] ]
        self.ii['bonds'] = [bond1]

    def _gen_ii_angles( self ):
        # angle1 = [ self.aLig2, self.aLig1, self.aProt1, 1, [self.ang1, 0.0, self.ang1, self.kAngle] ]
        # angle2 = [ self.aLig1, self.aProt1, self.aProt2, 1, [self.ang2, 0.0, self.ang2, self.kAngle] ]
        angle1 = [ self.aLig2, self.aLig1, self.aProt1, 1, [self.ang1, self.kAngle] ]
        angle2 = [ self.aLig1, self.aProt1, self.aProt2, 1, [self.ang2, self.kAngle] ]
        self.ii['angles'] = [ angle1, angle2 ]

    def _gen_ii_dihedrals( self ):
        # dihedral1 = [ self.aLig3, self.aLig2, self.aLig1, self.aProt1, 2, [self.dih1, 0.0, self.dih1, self.kDihedral] ]
        # dihedral2 = [ self.aLig2, self.aLig1, self.aProt1, self.aProt2, 2, [self.dih2, 0.0, self.dih2, self.kDihedral] ]
        # dihedral3 = [ self.aLig1, self.aProt1, self.aProt2, self.aProt3, 2, [self.dih3, 0.0, self.dih3, self.kDihedral] ]
        dihedral1 = [ self.aLig3, self.aLig2, self.aLig1, self.aProt1, 1, [self.dih1, 0.0, self.kDihedral] ]
        dihedral2 = [ self.aLig2, self.aLig1, self.aProt1, self.aProt2, 1, [self.dih2, 0.0, self.kDihedral] ]
        dihedral3 = [ self.aLig1, self.aProt1, self.aProt2, self.aProt3, 1, [self.dih3, 0.0, self.kDihedral] ]
        self.ii['dihedrals'] = [ dihedral1, dihedral2, dihedral3 ]

    def _renum_lig_top( self ):
        renumTop = cp.deepcopy(self.ligTop)
        protAtomNum = len(self.protTop.atoms)
        # atoms
        for a in renumTop.atoms:
            a.id += protAtomNum
            a.cgnr += protAtomNum
        return(renumTop)

    def _merge_top( self ):
        mergedTop = cp.deepcopy(self.protTop)
        protAtomNum = len(mergedTop.atoms)
        # atomtypes
        if self.ligTop.has_atomtypes==True:
            mergedTop.has_atomtypes=True
            for atype in self.ligTop.atomtypes:
                mergedTop.atomtypes.append(atype)
        if self.protTop.has_atomtypes==True:
            mergedTop.has_atomtypes=True
            for atype in self.protTop.atomtypes:
                mergedTop.atomtypes.append(atype)
        # atoms
        for a in self.ligTop.atoms:
            a.id += protAtomNum
            a.cgnr += protAtomNum
            mergedTop.atoms.append(a)
        # bonds
        for b in self.ligTop.bonds:
            mergedTop.bonds.append(b)
        # pairs
        for p in self.ligTop.pairs:
            mergedTop.pairs.append(p)
        # angles
        for ang in self.ligTop.angles:
            mergedTop.angles.append(ang)
        # dihedrals
        for dih in self.ligTop.dihedrals:
            mergedTop.dihedrals.append(dih)

        return(mergedTop)

    def _change_outfile_format(filename, ext):
        head, tail = os.path.split(filename)
        name, ex = os.path.splitext(tail)
        new_name = os.path.join(head, name+'.'+ext)
        return new_name

    def _random_protein_atoms( self, system ):
        n = len(system.atoms)
        rand = np.random.choice(n, size=3, replace=False)
        at = system.atoms
        a1 = at[rand[0]]
        a2 = at[rand[1]]
        a3 = at[rand[2]]
 
        indProt1 = a1.id-1
        indProt2 = a2.id-1
        indProt3 = a3.id-1

        return([indProt1],[indProt2],[indProt3])

    def _get_mass( self, aname ):
        aname = aname.translate(str.maketrans('','',digits))
        if aname.startswith('BR') or aname.startswith('Br') or aname.startswith('br'):
            return(pmxlib._atommass['BR'])
        elif aname.startswith('CL') or aname.startswith('Cl') or aname.startswith('cl'):
            return(pmxlib._atommass['CL'])
#        elif aname.startswith('D'):
#            return(pmxlib._atommass[aname[1]])
        elif aname.startswith('H'):
            return(0.0)
        elif aname.startswith('E'):
            return(0.0)
        else:
            try:
                return(pmxlib._atommass[aname[0]])
            except:
                return(1.0)
        return(1.0)

    def _get_atom_name( self, aname ):
        if self.bMartini==True:
            return('C')
        aname = aname.translate(str.maketrans('','',digits))
        if ('CL' in aname) or ('Cl' in aname) or ('cl' in aname):
            return('CL')
        if ('BR' in aname) or ('Br' in aname) or ('br' in aname):
            return('BR')
        return(aname[0])

    def _get_lig_COM( self, resname, system ):
        com = [0.0, 0.0, 0.0]
        mass = 0.0
        for a in system.atoms:
            if a.resname==resname and (a.id-1 in self.indLig):
                a.a2nm()
                amass = self._get_mass( a.name )
                com[0] += a.x[0]*amass
                com[1] += a.x[1]*amass
                com[2] += a.x[2]*amass
                mass += amass
        com[0] /= mass
        com[1] /= mass
        com[2] /= mass
        return(com)#,mass)

    def _closest_lig_atom_to_COM( self, system, resname, ligCOM ):
        indLig = -42
        mindist = 99999.999
        for a in system.atoms:
            if a.resname==resname and (a.id-1 in self.indLig):
                aname = self._get_atom_name( a.name )
                if aname.startswith('C') and not aname.startswith('CL'):
                    d = np.power(a.x[0]-ligCOM[0],2)+np.power(a.x[1]-ligCOM[1],2)+np.power(a.x[2]-ligCOM[2],2)
                    if d<mindist:
                        mindist = d
                        indLig = a.id-1
        return(indLig)

    def _furthest_lig_atom( self, system, resname, indLig1, indLig2=False, usedID=[], bReturnAtom=False ):
        indLig = -42
        maxdist = -99999.999
        maxatom = ''
        for a in system.atoms:
            if a.resname==resname and (a.id-1 in self.indLig):
                aname = self._get_atom_name( a.name )
                if aname.startswith('C') and not aname.startswith('CL') and (a.id-1 not in usedID):
                    d = a - system.atoms[indLig1]
                    if indLig2!=False:
                        d += a - system.atoms[indLig2]
                    if d>maxdist:
                        maxdist = d
                        indLig = a.id-1
                        maxatom = a
        if bReturnAtom==True:
            return(maxatom)
        return(indLig)

    def _closest_lig_atom( self, system, resname, indLig1, usedID=[], bReturnAtom=False ):
        indLig = -42
        mindist = 99999.999
        minatom = ''
        for a in system.atoms:
            if a.resname==resname and (a.id-1 in self.indLig):
                aname = self._get_atom_name( a.name )
                if aname.startswith(('C','N')) and not aname.startswith('CL') and (a.id-1 not in usedID):
                    d = a - system.atoms[indLig1]
                    if d<mindist:
                        mindist = d
                        indLig = a.id-1
                        minatom = a
        if bReturnAtom==True:
            return(minatom)
        return(indLig)

    def _mid_lig_atom( self, system, resname, ligCOM, indLig1, indLig2=False, usedID=[] ):
        # closest atom
#        minatom = self._closest_lig_atom( system, resname, indLig1, usedID, bReturnAtom=True)
        # furthest atom
        maxatom = self._furthest_lig_atom( system, resname, indLig1, indLig2=indLig2, usedID=usedID, bReturnAtom=True)
        # mid point
        ax = (ligCOM[0]+maxatom.x[0])*0.5
        ay = (ligCOM[1]+maxatom.x[1])*0.5
        az = (ligCOM[2]+maxatom.x[2])*0.5

        indLig = -42
        mindist = 99999.999
        for a in system.atoms:
            if a.resname==resname and (a.id-1 in self.indLig):
                aname = self._get_atom_name( a.name )
                if aname.startswith(('C','N')) and not aname.startswith('CL') and (a.id-1 not in usedID):
                    d = (a.x[0]-ax)**2 + (a.x[1]-ay)**2 + (a.x[2]-az)**2
                    if d<mindist:
                        mindist = d
                        indLig = a.id-1
        return(indLig)

    def _get_resname( self, system ):
        # don't use topology resname, as it may not match the one in the pdb
#        if self.ligTop != None:
#            return(self.ligTop.residues[0].resname)
#        else:
        for a in system.atoms:
            if (a.id-1 in self.indLig):
                return(a.resname)

    def _identify_ligand_atoms( self, system, indLig1 ):
        if indLig1==False:
            indLig1 = -42
        indLig2 = -42
        indLig3 = -42
        usedID = []

        # first atom closest to COM
#        resname = self.ligTop.residues[0].resname
        resname = self._get_resname( system )
        if indLig1==-42:
            ligCOM = self._get_lig_COM( resname, system ) # conversion a2nm is performed
            indLig1 = self._closest_lig_atom_to_COM( system, resname, ligCOM )
        usedID.append(indLig1)            

        if self.bLigMaxSpread==True:
            # second atom furthest from first atom
            indLig2 = self._furthest_lig_atom( system, resname, indLig1, usedID=usedID )
        elif self.bLigMidSpread==True:
            # second atom in the middle from the furthest and closest atom
            indLig2 = self._mid_lig_atom( system, resname, ligCOM, indLig1, usedID=usedID )
        else:
            # second atom closest to the first atom
            indLig2 = self._closest_lig_atom( system, resname, indLig1, usedID=usedID )
        usedID.append(indLig2)

        if self.bLigMaxSpread==True:
            # third atom furthest from first and second atoms
            indLig3 = self._furthest_lig_atom( system, resname, indLig1, indLig2=indLig2, usedID=[] )
        elif self.bLigMidSpread==True:
            indLig3 = self._mid_lig_atom( system, resname, ligCOM, indLig1, indLig2=indLig2, usedID=[] )
        else:
            indLig3 = self._closest_lig_atom( system, resname, indLig1, usedID=usedID )
       
        return([indLig1],[indLig2],[indLig3])

    def _ligAtomByName( self, atomName, system ):
        resname = self._get_resname( system )
        #resname = self.ligTop.residues[0].resname
        for a in system.atoms:
            if a.resname==resname and a.name==atomName:
                return( a.id-1 )
        return(False)

    def _dist_between_atoms( self, a1, a2 ):
        d = np.power(a1.x[0]-a2.x[0],2)+np.power(a1.x[1]-a2.x[1],2)+np.power(a1.x[2]-a2.x[2],2)
        return(np.sqrt(d))

    def _closest_prot_atom( self, system, indLig, usedID, bAngleCheck=0, indProt1=-42 ):
        mindist = 999999.999
        indProt = -42
        ligAtom = system.atoms[indLig]
        for a in system.atoms:
            a.a2nm()
            if (a.resname in library._protein_residues) or self.bMartini==True: # for martini it can be anything
                if (a.name=='CA' or a.name=='BB' or self.bMartini==True) and (a.id-1 not in usedID) and (a.id-1 in self.indProt):
                    d = a - ligAtom 
                    bAngle = True
                    if self.strFileApo != None: # check the distance against the closest apo structure atom
                        atomApoInd,atomApo,apoDist = self._closest_apo_atom_dist(a)
                        if apoDist > self.apoDistanceLimit:
                            continue
                    if bAngleCheck==1: # check for angle 1
                        ang = ligAtom.angle(self.aLig2,a,degree=True)
                        bAngle = self._check_angle( ang )
                    elif bAngleCheck==2: # check for angle 2
                        aProt1 = system.atoms[indProt1]
                        ang = aProt1.angle(self.aLig1,a,degree=True)
                        bAngle = self._check_angle( ang )
                    if d<mindist and bAngle:
                        mindist = d
                        indProt = a.id-1
        return(indProt)

    def _closest_apo_atom_dist( self, aholo ):
        # find the distance to the closest apo atom
        d = _find_atom_inY( aholo.id, self.system, self.systemApo, bReturnDist=True ) 
        return(d)

    def _identify_protein_atoms( self, system, indLig1, indLig2, indLig3 ):
        indProt1 = -42
        indProt2 = -42
        indProt3 = -42
        indLig1 = indLig1[0]
        indLig2 = indLig2[0]
        indLig3 = indLig3[0]
        usedID = []

        # indProt1 prot atom closest to indLig1
        indProt1 = self._closest_prot_atom( system, indLig1, usedID, bAngleCheck=1 ) 
        usedID.append(indProt1)

        # indProt2 prot atom closest to indLig2
#        indProt2 = self._closest_prot_atom( system, indLig2, usedID ) 
        indProt2 = self._closest_prot_atom( system, indProt1, usedID ) 
        usedID.append(indProt2)

        # indProt3 prot atom closest to indLig3
#        indProt3 = self._closest_prot_atom( system, indLig3, usedID, bAngleCheck=2, indProt1=indProt1 ) 
        indProt3 = self._closest_prot_atom( system, indProt2, usedID, bAngleCheck=2, indProt1=indProt1 ) 
        usedID.append(indProt3)
#        print indLig1, indProt1, indProt2, indProt3

        return([indProt1],[indProt2],[indProt3])

    def _identify_protein_atomsAlt( self, system, indLig ):
        indProt1 = -42
        indProt2 = -42
        indProt3 = -42
        ligAtom = system.atoms[indLig[0]]
        mindist = 99999.999
        minres = -42
        for a in system.atoms:
            if a.resname in library._protein_residues:
                if a.name=='N' or a.name=='CA' or a.name=='C' or a.name=='BB':
                    d = a-ligAtom
                    if d<mindist:
                        mindist = d
                        indProt1 = a.id-1
                        minres = a.resnr-1

        for a in system.residues[minres].atoms:
            if a.name=='N' or a.name=='CA' or a.name=='C' or a.name=='BB':
                if indProt2==-42:
                    if indProt1!=a.id-1:
                        indProt2 = a.id-1
                elif indProt3==-42:
                    if indProt1!=a.id-1:
                        indProt3 = a.id-1

        return([indProt1],[indProt2],[indProt3])


    def _check_ndx_size( self,ndxGroupSize, desiredSize ):
        if ndxGroupSize!=desiredSize:
            sys.stdout.write('Need to select a group with %d atoms, but the selected group contained %d atoms. Exiting...\n' % (desiredSize,ndxGroupSize) )
            sys.exit(0)


    def _select_ndx( self, fname="index.ndx", message=False ):
        if not os.path.isfile( fname ):
            return(False)
        ndx_file = ndx.IndexFile( fname )
        names = ndx_file.names
        i = 0
        ndxDict = {}
        for name in names:
            atomNum = len(ndx_file[name].ids)
            sys.stdout.write('%d %s: %d atoms\n' % (i,name,atomNum) )
            ndxDict[i] = ndx_file[name].ids
            i+=1
        sys.stdout.write('\n')

        if message==False:
            sys.stdout.write('Select a group for analysis:\n')
        else:
            sys.stdout.write(message+'\n')

        ndxNum = -1
        while ndxNum==-1:
            ndxNum = input()
            if ndxNum.isdigit()==False:
                sys.stdout.write('Wrong index group number selected (use an integer number)\n')
                ndxNum = -1
                continue

            ndxNum = int(ndxNum)
            if (ndxNum >= i) or (ndxNum < 0):
                sys.stdout.write('Wrong index group number selected\n')
                ndxNum = -1
        sys.stdout.write('Selected group %d\n\n' % ndxNum)

        res = []
        res = np.asarray(ndxDict[ndxNum])-1 # starting from 0
        return(res)

########################################################################
####### functions to transfer the restraints between two systems #######
########################################################################
def _find_atom_inY( ind, strX, strY, bName=True, bResName=False, indMap={}, bReturnDist=False ):
    """Given an atom index in strX, identify the closest atom in strY.
    Optionally, can also require atom names to match.

    Parameters
    ----------
    ind: int
        index in strX    
    strX: Model
        reference structure
    stryY: Model
        structure to search in
    indMap: dict
        dictionary of index maps indMap[indX]=indY

    Attributes
    ----------
    bName: bool
        should atom names also match. Default True
    bResName: bool
        should residue names also match. Default False
    bReturnDist: bool
        return distance between atoms. Default False

    Returns
    -------
    index: int
        index of the identified atom
    atom: atoms
        identified atom in the Model
    """

    aX = strX.atoms[ind-1]
    aX.a2nm()
    aY = 0
    minDist = 999.99
    for a in strY.atoms:
        a.a2nm()
        if bName==True and a.name!=aX.name:
            continue
        if bResName==True and not _compare_resnames(a.resname,aX.resname):
            continue
        if a.id in indMap.values():
            continue
        d = a-aX
        if d<minDist:
            minDist = d
            aY = a
    if bReturnDist==True:
        return(aY.id,aY,minDist)
    else:
        return(aY.id,aY)

def _compare_resnames(rname1,rname2):
    if rname1==rname2:
        return(True)
    elif rname1.startswith('H') and rname2.startswith('H'):
        return(True)
    else:
        return(False)

def transfer_ii( pdbX, pdbY, iiInpFname, iiOutFname, bName=True, bResName=True, pdbRestrAtoms=None ):
    """Transfers intermolecular restraints from one structure to another,
    i.e. from pdbX to pdbY. This is done by searching for closest atoms
    in the structures.

    Parameters
    ----------
    pdbX: pdb file
        structure for which restraints are known
    pdbY: pdb file
        structure for which to create restraints
    iiInpFname: restraint file
        restraints for pdbX
    iiOutFname: restraint file
        output restraint file
    pdbRestrAtoms: pdb output file, optional
        output atoms that are restrained

    Attributes
    ----------
    bName: bool
        should atom names also match. Default True
    bResName: bool
        should residue names also match. Default True
    """

    strX = Model(pdbX)
    strY = Model(pdbY)

    fp = open(iiInpFname,'r')
    lines = fp.readlines()
    fp.close()

    fp = open(iiOutFname,'w')
    atomsY = [] 
    indMap = {} # map of identified indeces: indMap{indX} = indY
    for l in lines:
        if '[' in l:
            fp.write(l)
        elif l=='\n':
            fp.write(l)
        else:
            l = l.rstrip()
            l = l.lstrip()
            foo = l.split()
#             print(foo,len(foo))

            # go over the array in reverse order to get atomID
            out = []
            for i in range(1,len(foo)+1):
                if i>5:
                    atomIndX = strX.atoms[int(foo[-i])-1].id
                    if atomIndX in indMap.keys():
                        atomIndY = indMap[atomIndX]
                    else:
                        atomIndY,atomY = _find_atom_inY( int(foo[-i]), strX, strY, bName=bName, bResName=bResName, indMap=indMap )
                        atomsY.append(atomY)
                        indMap[atomIndX] = atomIndY
                    out.insert(0,atomIndY)
                else:
                    out.insert(0,foo[-i])

            # go over the elements again to output them
            fp.write('   ')
            for i in range(0,len(foo)):
                fp.write('%s   ' % str(out[i]))
            fp.write('\n')
             
    fp.close()

    # optional output of restrained atoms as .pdb
    if pdbRestrAtoms!=None:
        mout = Model()
        mout.atoms.extend( atomsY )
        mout.write(pdbRestrAtoms)


class PositionRestraints:
    def __init__(self, **kwargs):
        self.num = -1
        self.resnum = -1
        self.name = ''
        self.fcx = 1000.0
        self.fcy = 1000.0
        self.fcz = 1000.0
        self.itpname = None
        self.outname = None
        self.residfile = None
        self.atomidfile = None
        self.pdbfile = None
        self.bHeavy = True # if to consider only heavy atoms
        self.bAllHeavy = False # restrain all heavy atoms
        self.bPosreIfdef = True # should ifdef around posres be added
        self.bVerbose = False
        self.stateBonded = 'A' # may need to define this to be stateBonded='A' when dealing with cgenff
        self.bIgnoreIfdefs = False # add posre ignoring that it may already be defined inside ifdef

        for key, val in kwargs.items():
            setattr(self,key,val)

        self.position_restraints(self.num, self.resnum, self.name, self.fcx, self.fcy, self.fcz, self.itpname, self.outname, self.residfile, self.atomidfile, self.pdbfile)

    def _posre_exists(self,topol,anum,resnum=-1):
        if topol.has_posre:
            for posre in topol.posre:
                at = posre[0]
                if anum==at.id:
                    if resnum==-1:
                        return(True)
                    elif resnum==at.resnr:
                        return(True)
        if self.bIgnoreIfdefs==True:
            return(False)

        # also, posre can be in the footer, so check there
        bPosre = False
        for foot in topol.footer:
            if foot.startswith('#') or foot.startswith(';') or foot=='':
                continue
            if foot.startswith('['):
                if 'position_restraints' in foot:
                    bPosre = True
                    continue
                else:
                    bPosre = False
            if bPosre==True:
               foo = foot.lstrip().rstrip().split()
               if anum==int(foo[0]):
                   return(True)

        return(False)

    def _append_posre(self,topol,at,fvec):
        rest = ' '.join([str(x) for x in fvec])
        posre =  [at, 1, rest] 
        topol.posre.append(posre)
        topol.has_posre = True

    def _add_posre(self, topol,fvec,anum=-1,aname='',resnum=-1,residlist=[],atomidlist=[]):
        ### for posre based on id list files
        # resid
        if len(residlist)>0:
            for a in topol.atoms:
                if len(aname)>0 and a.name not in aname:
                    continue
                for rnum in residlist:
                    if rnum==a.resnr:
                        if self._posre_exists(topol,a.id,rnum)==False:
                            self._append_posre(topol,a,fvec)
        # atomid 
        if len(atomidlist)>0:
            for aid in atomidlist:
                a = topol.atoms[aid-1]
                if self._posre_exists(topol,a.id)==False:
                    self._append_posre(topol,a,fvec)

        # all heavy atoms
        if self.bAllHeavy==True:
            for a in topol.atoms:
                if a.name.startswith('H') or a.name.startswith('D'):
                    continue
                if not self._posre_exists(topol,a.id):
                    self._append_posre(topol,a,fvec)

        ### for single posre: by atom name and residue id
        if anum!=-1: # by atom number
            a = topol.atoms[anum-1]
            if resnum!=-1:
                if resnum==a.resnr:
                    if not self._posre_exists(topol,anum,resnum):
                        self._append_posre(topol,a,fvec)
            else:
                if not self._posre_exists(topol,anum,resnum):
                    self._append_posre(topol,a,fvec)
        elif len(residlist)==0 and len(atomidlist)==0: # by atom name
            for a in topol.atoms:
                if a.name in aname:
                    if resnum!=-1:
                        if resnum==a.resnr:
                            if self._posre_exists(topol,a.id,resnum):
                                break
                            else:
                                self._append_posre(topol,a,fvec)
                    else:
                        if not self._posre_exists(topol,a.id):
                            self._append_posre(topol,a,fvec)

    def _is_perturbed_res( self, residue ):
        for atom in residue.atoms:
            if atom.atomtypeB is not None and (atom.q!=atom.qB or atom.m != atom.mB or atom.atomtype != atom.atomtypeB): return True
        return False
 
    def _read_id_file( self, fname ):
        fp = open(fname,'r')
        lines = fp.readlines()
        fp.close()
        out = []
        for l in lines:
            l.rstrip()
            foo = l.split()
            for i in foo:
                out.append(int(i))
        return(out)

    def _is_heavy( self, a ):
        mass = self._get_mass( a.name )
        if mass > 1.0:
            return(True)
        else:
            return(False)

    def _get_mass( self, aname ):
        aname = aname.translate(str.maketrans('','',digits)) 
        if aname.startswith('Br') or aname.startswith('BR'):
            return(pmxlib._atommass['BR'])
        elif aname.startswith('Cl') or aname.startswith('CL'):
            return(pmxlib._atommass['CL'])
        elif aname.startswith('Mg') or aname.startswith('MG'):
            return(pmxlib._atommass['MG'])
        elif aname.startswith('D'):
            return(pmxlib._atommass[aname[1]])
        elif aname.startswith('EP'):
            return(0.0)
        else:
            return(pmxlib._atommass[aname[0]])
        return(1.0)

    def _get_com(self,m):
        com = [0.0, 0.0, 0.0]
        mass = 0.0

        # com
        for a in m.atoms:
            a.a2nm()
            if ('HOH' in a.resname) or ('SOL' in a.resname) or ('WAT' in a.resname):
                continue
            aname = a.name[0]
            if aname.isdigit():
                aname = a.name[1]
            amass = self._get_mass( a.name )

            com[0] += a.x[0]*amass
            com[1] += a.x[1]*amass
            com[2] += a.x[2]*amass
            mass += amass

        com[0] /= mass
        com[1] /= mass
        com[2] /= mass

        return(com)

    def _find_closest_to_com( self, com, m, bName=False, name='' ):
        mind = 99999.999
        anum = -1
        resnum = -1
        for a in m.atoms:
            if bName==True:
                if a.name != name[0]:
                    continue
            if self.bHeavy==True:
                if a.name.startswith('H') or a.name.startswith('D'):
                    continue
            d = np.sqrt( (a.x[0]-com[0])**2 + (a.x[1]-com[1])**2 + (a.x[2]-com[2])**2)
            if d < mind:
                if self._is_heavy( a ):
                    mind = d
                    anum = a.id
                    resnum = a.resnr
        return(anum)#,resnum)

    def position_restraints(self, num, resnum, name, fcx, fcy, fcz, itpname, outname, residfile, atomidfile, pdbfile):
        bNum = False
        bName = False
        residlist = []
        atomidlist = []
        if num > -1:
            bNum = True
        elif name != '':
            bName = True
            name = name.split()
        if residfile is not None:
            residlist = self._read_id_file( residfile )
        if atomidfile is not None:
            atomidlist = self._read_id_file( atomidfile )
        if pdbfile is not None:
            pdb = Model( pdbfile,bPDBTER=True,renumber_residues=False )
            if self.bVerbose==True:
                sys.stdout.write("pdb file was provided (-ipdb option), therefore other options are ignored.")
                sys.stdout.write("Will calculate COM and identify the heavy atom closest to it.")
                sys.stdout.write("If atom name is provided, will identify atom by name closest to the COM.")

        if pdbfile is None:
            if bNum==False and bName==False and len(residlist)==0 and len(atomidlist)==0:
                sys.stdout.write("You need to provide an atom number or an atom name or id list files")
                sys.exit(0)

        if self.bVerbose==True:
            sys.stdout.write('Reading input .itp file "%s""' % (itpname) )
#        topol = Topology( itpname,  topfile = None, version = 'new', assign_types=False )
        topol = TopolBase( itpname )

        fvec = [fcx, fcy, fcz]
        if pdbfile is not None:
            com = self._get_com( pdb )
            anum = self._find_closest_to_com( com, pdb, bName, name )
            self._add_posre( topol,fvec,anum=anum,resnum=resnum )
        elif len(residlist)>0 or len(atomidlist)>0:
            if len(residlist)>0:
                self._add_posre(topol,fvec,aname=name,residlist=residlist)
            if len(atomidlist)>0:
#                print atomidlist
                self._add_posre(topol,fvec,aname=name,atomidlist=atomidlist)
        elif bNum:
            self._add_posre(topol,fvec,anum=num,resnum=resnum)
        else:
            self._add_posre(topol,fvec,aname=name,resnum=resnum)
          

#        chargeB = sum([a.qB if a.qB is not None else a.q for a in topol.atoms])
        chargeB = 0.0 
        for r in topol.residues:
            if self._is_perturbed_res( r ):
                for a in r.atoms:
                    if a.qB is not None:
                        chargeB += a.qB
                    else:
                        chargeB +=  a.q
        topol.write(outname, target_qB=[chargeB], posre_ifdef=self.bPosreIfdef, stateBonded=self.stateBonded)

