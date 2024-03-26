import sys, os, time, math
import numpy as np
from copy import deepcopy
from pmx import *
from pmx import library as pmxlib
from string import digits
from pmx.utils import doLog

class DoubleBox:
    def __init__(self, **kwargs):
        self.r = 2.5
        self.d = 1.5
        self.pdb1 = None
        self.pdb2 = None
        self.outfile = None
        self.bLongestAxis = None
        self.logfile = None

        self.bPython3 = False
        if (sys.version_info > (3, 0)):
            self.bPython3 = True

        for key, val in kwargs.items():
            setattr(self,key,val)

#        self.double_box(self.r, self.d, self.pdb1, self.pdb2, self.outfile )
        self.double_box( ) 

    def _translate(self,m,v,fact=1.0):
        for a in m.atoms:
            a.x[0] = a.x[0] + fact*v[0]
            a.x[1] = a.x[1] + fact*v[1]
            a.x[2] = a.x[2] + fact*v[2]

    def _get_mass( self, a ):
        if self.bPython3:
            aname = a.name.translate(str.maketrans('','',digits))
        else:
            aname = a.name.translate(None, digits)

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
            try:
                return(pmxlib._atommass[aname[0]])
            except:
                return(1.0)
        return(1.0)

    def _principal_axes( self, m ):
        tensor = np.matrix([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        for a in m.atoms:
            if ('HOH' in a.resname) or ('SOL' in a.resname) or ('WAT' in a.resname):
                continue
            mass = self._get_mass( a )
            tensor[0,0] += mass*( np.power(a.x[1],2) + np.power(a.x[2],2) )
            tensor[1,1] += mass*( np.power(a.x[0],2) + np.power(a.x[2],2) )
            tensor[2,2] += mass*( np.power(a.x[0],2) + np.power(a.x[1],2) )
            tensor[0,1] -= mass*a.x[0]*a.x[1]
            tensor[1,0] = tensor[0,1]
            tensor[0,2] -= mass*a.x[0]*a.x[2]
            tensor[2,0] = tensor[0,2]
            tensor[1,2] -= mass*a.x[1]*a.x[2]
            tensor[2,1] = tensor[1,2]
        evals,evecs = np.linalg.eig(tensor)
        idx = evals.argsort()[::1] # sort descending, because later the rotation matrix will be calculated as transpose
        evals = evals[idx]
        evecs = evecs[:,idx]
        # check if one axis needs to be flipped
        crossprod = np.cross(np.transpose(evecs[:,0]),np.transpose(evecs[:,1]))
        dotprod = np.dot(np.transpose(evecs[:,2]),np.transpose(crossprod))
        if dotprod < 0.0:
            evecs[:,2] *= -1.0
        rotmat = np.transpose(evecs)
        return(rotmat)

    def _rotate( self, atoms, R):
        for atom in atoms:
            x_old = map(lambda x: x, atom.x)
            if self.bPython3==True:
                x_old = list(x_old)
            for i in range(3):
                atom.x[i] = 0.0
                for j in range(3):
                    atom.x[i] += x_old[j]*R[i,j]
#            for r in range(3):
#                atom.x[r] = 0
#                for c in range(3):
#                    atom.x[r]+=R[r,c]*x_old[c]

    def _get_com_radius(self,m):
        com = [0.0, 0.0, 0.0]
        radius = 0.0
        mass = 0.0

        # com
        for a in m.atoms:
            a.a2nm()
            if ('HOH' in a.resname) or ('SOL' in a.resname) or ('WAT' in a.resname):
                continue
            aname = a.name[0]
            if aname.isdigit():
                aname = a.name[1]
            amass = self._get_mass( a )
	
            com[0] += a.x[0]*amass
            com[1] += a.x[1]*amass
            com[2] += a.x[2]*amass
            mass += amass

        com[0] /= mass
        com[1] /= mass
        com[2] /= mass

        # radius
        for a in m.atoms:
            if ('HOH' in a.resname) or ('SOL' in a.resname) or ('WAT' in a.resname):
                continue
            foo = (a.x[0]-com[0])**2 + (a.x[1]-com[1])**2 + (a.x[2]-com[2])**2
            if foo>radius:
                radius = foo
        radius = np.sqrt(radius)

        return(com,radius)

    def _get_rect_dist( self, m ):
        minx = 9999.999
        maxx = -9999.999
        miny = 9999.999
        maxy = -9999.999
        minz = 9999.999
        maxz = -9999.999
        for a in m.atoms:
            if ('HOH' in a.resname) or ('SOL' in a.resname) or ('WAT' in a.resname):
                continue
            if a.x[0] < minx:
                minx = a.x[0]
            if a.x[0] > maxx:
                maxx = a.x[0]
            if a.x[1] < miny:
                miny = a.x[1]
            if a.x[1] > maxy:
                maxy = a.x[1]
            if a.x[2] < minz:
                minz = a.x[2]
            if a.x[2] > maxz:
                maxz = a.x[2]
        a = maxx-minx
        b = maxy-miny
        c = maxz-minz
        return(a,b,c)

    def _get_extr( self, m, i=0 ):
        minx = 9999.999
        maxx = -9999.999
        for a in m.atoms:
            if ('HOH' in a.resname) or ('SOL' in a.resname) or ('WAT' in a.resname):
                continue
            if a.x[i] < minx:
                minx = a.x[i]
            if a.x[i] > maxx:
                maxx = a.x[i]
        return(minx,maxx)

    def _translate_sticking_out( self, m, l, i=0 ):
        minx = 0.0
        maxx = 0.0
        for a in m.atoms:
            if ('HOH' in a.resname) or ('SOL' in a.resname) or ('WAT' in a.resname):
                continue
            if a.x[i] < minx:
                minx = a.x[i]
            if (a.x[i]-l) > maxx:
                maxx = a.x[i]-l
        return( -1.0*(minx+maxx) )
        

    def double_box(self): #, c, d, pdb1, pdb2, outfile,bLongestAxis):
        m1 = Model(self.pdb1,bPDBTER=False,renumber_residues=False)
        m2 = Model(self.pdb2,bPDBTER=False,renumber_residues=False)

        # com and radii
        com1,rad1 = self._get_com_radius(m1)
        com2,rad2 = self._get_com_radius(m2)

        # remove COM from the systems
        self._translate(m1,com1,fact=-1.0)
        self._translate(m2,com2,fact=-1.0)
 
        # if needed, calculate principal axes
        if self.bLongestAxis:
            princRot1 = self._principal_axes( m1 )
            princRot2 = self._principal_axes( m2 )
            self._rotate( m1.atoms, princRot1 )
            self._rotate( m2.atoms, princRot2 )
            a_rect1,b_rect1,c_rect1 = self._get_rect_dist( m1 )
            a_rect2,b_rect2,c_rect2 = self._get_rect_dist( m2 )

        # estimate cube's edge for the larger structure
        a_cube = 0.0
        if rad1>rad2:
            a_cube = 2*rad1+2*self.d
        else:
            a_cube = 2*rad2+2*self.d
        if self.bLongestAxis:
            a_rect = a_rect1 + a_rect2 + self.r + 2.0*self.d
            b_rect = b_rect1 + b_rect2 + 2.0*self.d
            c_rect = c_rect1 + c_rect2 + 2.0*self.d
            if self.logfile!=None:
                doLog(self.logfile,"Cuboid dimensions (nm): {0} {1} {2}".format(np.round(a_rect),np.round(b_rect),np.round(c_rect,2)))
        else:
            if self.logfile!=None:
                doLog(self.logfile,"Cube's edge (nm): {0}".format(a_cube))

        # cube's diagonal
        if self.bLongestAxis:
            d_cube = np.sqrt( np.power(a_rect,2) + np.power(b_rect,2) + np.power(c_rect,2) )
            if self.logfile!=None:
                doLog(self.logfile,"Cuboid's diagonal (nm): ".format(np.round(d_cube,2)))
        else:
            d_cube = a_cube*np.sqrt(3.0)#/2.0
            if self.logfile!=None:
                doLog(self.logfile,"Cube's diagonal (nm): {0}".format(np.round(d_cube,2)))

        # check if cube is enough
        if self.bLongestAxis:
            a_box = a_rect
            b_box = b_rect
            c_box = c_rect
        else:
            dist_cube = d_cube - 2.0*rad1 - 2.0*rad2 - self.r - 2.0*self.d
            a_box = a_cube
            b_box = a_cube
            c_box = a_cube
            if dist_cube > 0.0:
                if self.logfile!=None:
                    doLog(self.logfile,"Having a cube with the diagonal {0} nm is good enough".format(np.round(d_cube,2)))
            else:
                if self.logfile!=None:
                    doLog(self.logfile,"Need to extend the cube")
                #d_rect = 2.0*rad1 + 2.0*rad2 + self.r + 2.0*self.d
                d_rect = 2.0*rad1 + 2.0*rad2 + self.r + 1.0*self.d
                #delta_a_cube = -a_cube + np.sqrt(d_rect**2 - 2.0*a_cube**2)
                delta_a_cube = (d_rect/np.sqrt(3.0)) - a_cube
                a_box = a_cube + delta_a_cube
                if self.logfile!=None:
                    doLog(self.logfile,"Rectangle's edge (nm): ".format(np.round(a_box,2)))

        # translate the larger structure to the middle of the box
        if rad1>rad2:
#            print "Translating structure 1 to the middle of the box"
            # translate smaller to (0,b/2,c/2) and larger to (maxx2+dist+abs(minx2),b/2,c/2)
            if self.bLongestAxis:
                self._translate(m2,[0.0,b_box/2.0,c_box/2.0])
                minx1,maxx1 = self._get_extr( m1, 0 )
                minx2,maxx2 = self._get_extr( m2, 0 )
                self._translate(m1,[maxx2+c+abs(minx1),b_box/2.0,c_box/2.0])
            else:
                #self._translate(m1,[a_box/2.0,b_box/2.0,c_box/2.0])
                self._translate(m1,[(rad1+rad2+self.r)/np.sqrt(3.0),(rad1+rad2+self.r)/np.sqrt(3.0),(rad1+rad2+self.r)/np.sqrt(3.0)])

        else:
#            print "Translating structure 2 to the middle of the box"
            if self.bLongestAxis:
                self._translate(m1,[0.0,b_box/2.0,c_box/2.0])
                minx1,maxx1 = self._get_extr( m1, 0 )
                minx2,maxx2 = self._get_extr( m2, 0 )
                self._translate(m2,[maxx1+c+abs(minx2),b_box/2.0,c_box/2.0])
            else:
                #self._translate(m2,[a_box/2.0,b_box/2.0,c_box/2.0])
                self._translate(m2,[(rad1+rad2+self.r)/np.sqrt(3.0),(rad1+rad2+self.r)/np.sqrt(3.0),(rad1+rad2+self.r)/np.sqrt(3.0)])

        # create output
        mout = m1
        chainList = []
        waterIonsString = ['HOH','SOL','WAT','K','Na','KJ','NaJ','NA','SOD','Cl','CL','CLJ','MG','Mg','Ca','CA']
        for ch in m1.chains:
            if ch.id not in chainList:
                bWaterIons = True
                for r in ch.residues:
                    if r.resname not in waterIonsString:
                        bWaterIons = False
                        break
                if bWaterIons==False:
                    chainList.append(ch.id)
        chainIDstring = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
        chainIDchanges = {}
        for ch in m2.chains:
            if ch.id in chainList:
	        # pick a new chain ID
                bFound = False
                while bFound==False:
                    foo = chainIDstring[0]
                    chainIDstring = chainIDstring.lstrip(chainIDstring[0])
                    if foo not in chainList:
                        bFound=True
                        chainList.append(foo)
                        chainIDchanges[ch.id] = foo
        for a in m2.atoms:
            chID = a.chain_id
            if chID in chainIDchanges:
                a.chain_id = chainIDchanges[chID]
        mout.atoms.extend(m2.atoms)

        # create the box
        mout.box = [[a_box,0.0,0.0],[0.0,b_box,0.0],[0.0,0.0,c_box]]

        # determine translation exactly by checking how much the atoms are sticking out
        xtransl = self._translate_sticking_out( mout, a_box, 0 )
        ytransl = self._translate_sticking_out( mout, b_box, 1 )
        ztransl = self._translate_sticking_out( mout, c_box, 2 )
        self._translate(mout,[3.00*xtransl,3.00*ytransl,3.00*ztransl]) # increase translation by a factor 1.01 to move a bit further from walls

        mout.write(self.outfile,bPDBTER=True)

        # volume
        V = a_box*b_box*c_box
        if self.logfile!=None:
            doLog("Volume: {0} nm^3".format(np.round(V,2)))


