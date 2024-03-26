#!/usr/bin/env python

"""This module contains wrapper functions for some of the most used Gromacs
tools.
"""

import os
from .utils import which
import subprocess
from .library import mdps
import subprocess


def get_gmx():
    """Gets path to gmx executable, and throws error if not found.
    """
    gmx = which('gmx_mpi') if which('gmx') is None else which('gmx')
    if gmx is not None:
        return gmx
    else:
        raise EnvironmentError('gmx executable not found')


def set_gmxlib():
    """Sets the environment variable GMXLIB to the default/expected pmx
    location.
    """
    path = os.path.abspath(__file__)
    dir_path = os.path.dirname(path)
    gmxlib = os.path.join(dir_path, 'data/mutff')
    os.environ['GMXLIB'] = gmxlib


def editconf(f, o='editconf.gro', bt='cubic', d=1.2, other_flags='', verbose=True, gmxexec=None):
    """Simple ``gmx editconf`` wrapper.

    Parameters
    ----------
    f : str
        input structure file
    o : str, optional
        name of output structure file. Default is "solvate.gro"
    bt : str
        box type:triclinic, cubic, dodecahedron, or octahedron
    d : float
        distance between the solute and the box (nm)
    verbose : bool
        full gromacs output (True)
    gmxexec : str
        gmx with the full path
    other_flags : str, optional
        additional flags to pass as you would type them in the shell

    Returns
    -------
    None

    """

    if gmxexec is None:
        gmxexec = get_gmx()

    if verbose:
        capture_output = False
    else:
        capture_output = True

    subprocess.run('{gmxexec} editconf -f {f} -o {o} -bt {bt} -d {d} '
         '{other_flags}'.format(gmxexec=gmxexec, f=f, o=o, bt=bt, d=d, other_flags=other_flags),
         shell=True, check=True, capture_output=capture_output)


def pdb2gmx(f, o='pdb2gmx.gro', p='topol.top', ff='amber99sb-star-ildn-mut',
            water='tip3p', other_flags='', verbose=True, gmxexec=None):
    """Simple ``gmx pdb2gmx`` wrapper.

    Parameters
    ----------
    f : str
        input structure file
    o : str, optional
        output structure file. Default is "pdb2gmx.gro".
    p : str, optional
        name of topology file. Default is "topol.top".
    ff : str, optional
        forcefield. Default is "amber99sb-star-ildn-mut".
    water : str, optional
        water model. Default is "tip3p".
    verbose : bool
        full gromacs output (True)
    gmxexec : str
        gmx with the full path
    other_flags : str, optional
        additional flags to pass as you would type them in the shell

    Returns
    -------
    None

    """

    if gmxexec is None:
        gmxexec = get_gmx()

    if verbose:
        capture_output = False
    else:
        capture_output = True

    subprocess.run('{gmxexec} pdb2gmx -f {f} -o {o} -p {p} -ff {ff} -water {water} '
         '{other_flags}'.format(gmxexec=gmxexec, f=f, o=o, p=p, ff=ff, water=water, other_flags=other_flags), 
         shell=True, check=True, capture_output=capture_output)


def solvate(cp, cs='spc216.gro', p='topol.top', o='solvate.gro',
            other_flags='', verbose=True, gmxexec=None):
    """Simple ``gmx solvate`` wrapper.

    Parameters
    ----------
    cp : str
        input structure file
    cs : str, optional
        structure file of the solvent. Default is "spc216.gro".
    p : str, optional
        name of topology file. Default is "topol.top"
    o : str, optional
        name of output structure file. Default is "solvate.gro"
    verbose : bool
        full gromacs output (True)
    gmxexec : str
        gmx with the full path
    other_flags : str, optional
        additional flags to pass as you would type them in the shell

    Returns
    -------
    None

    """

    if gmxexec is None:
        gmxexec = get_gmx()

    if verbose:
        capture_output = False
    else:
        capture_output = True

    subprocess.run('{gmxexec} solvate -cp {cp} -cs {cs} -p {p} -o {o} '
         '{other_flags}'.format(gmxexec=gmxexec, cp=cp, cs=cs, p=p, o=o, other_flags=other_flags),
         shell=True, check=True, capture_output=capture_output)


def grompp(f, c, p, o='grompp.tpr', maxwarn=0, other_flags='', verbose=True, gmxexec=None):
    """Simple ``gmx grompp`` wrapper.

    Parameters
    ----------
    f : str
        input mdp file
    c : str
        input structure file
    p : str
        input topology file. Default is "topol.top"
    o : str, optional
        output tpr file. Default is "grompp.tpr"
    verbose : bool
        full gromacs output (True)
    gmxexec : str
        gmx with the full path
    maxwarn : int, optional
        number of allowed warnings. Default is 0.
    other_flags : str, optional
        additional flags to pass as you would type them in the shell

    Returns
    -------
    None

    """

    if gmxexec is None:
        gmxexec = get_gmx()

    if verbose:
        capture_output = False
    else:
        capture_output = True

    subprocess.run('{gmxexec} grompp -f {f} -c {c} -r {c} -p {p} -o {o} -maxwarn {maxwarn}'
         '{other_flags}'.format(gmxexec=gmxexec, f=f, c=c, p=p, o=o, maxwarn=maxwarn, other_flags=other_flags), 
         shell=True, check=True, capture_output=capture_output)


def genion(s, p, o='genion.gro', np=0, nn=0, conc=0.15, neutral=True,
           other_flags='', verbose=True, gmxexec=None):

    """Simple ``gmx genion`` wrapper. By default, group 'SOL' will be replaced
    by ions.

    Parameters
    ----------
    s : str
        input tpr file
    p : str
        input topology file.
    o : str, optional
        name of output structure file. Default is "genion.gro"
    np : int, optional
        number of positive ions. Default is 0.
    nn : int, optional.
        number of negative ions. Default is 0.
    conc : float, optional
        specify salt concentration (mol/liter). Default is 0.15 M.
    neutral : bool
        whether to add enough ions to neutralise the system. These
        ions are added on top of those specified with -np/-nn or -conc.
        Default is True.
    verbose : bool
        full gromacs output (True)
    gmxexec : str
        gmx with the full path
    other_flags : str, optional
        additional flags to pass as you would type them in the shell

    Returns
    -------
    None

    """

    if gmxexec is None:
        gmxexec = get_gmx()

    if neutral is True:
        other_flags += ' -neutral'

    if verbose:
        capture_output = False
    else:
        capture_output = True

    
    subprocess.run('echo "SOL" | {gmxexec} genion -s {s} -p {p} -o {o} -np {np} -nn {nn} -conc {conc} '
         '{other_flags}'.format(gmxexec=gmxexec, s=s, p=p, o=o, np=np, nn=nn, conc=conc, other_flags=other_flags), 
         shell=True, check=True, capture_output=capture_output)


def trjconv(f, s, o='trjconv.xtc', ur='compact', pbc='none', fit='none',
            out_grp='System', fit_grp='C-alpha', sep=False, other_flags='', verbose=True, gmxexec=None):
    """Simple ``gmx trjconv`` wrapper.

    Parameters
    ----------
    f : str
        input structure of trajectory file
    s : str
        input tpr file
    o : str, optional
        output trajectory/structure file. Default is "trjconv.xtc"
    ur : str, optional
        unit-cell representation: rect, tric, compact. Default is 'compact'.
    pbc : str, optional
        PBC treatment: none, mol, res, atom, nojump, cluster, whole.
        Default is 'none'.
    fit : str, optional
        fit molecule to ref structure in the structure file: none, rot+trans,
        rotxy+transxy, translation, transxy, progressive.
        Default is 'none'.
    out_grp : str, optional
        output group. Defauls is 'System'.
    fit_grp : str, optional
        group to use for the fitting if 'fit' is not none.
        Default is 'C-alpha'.
    sep : bool, optional
        write each frame to a separate .gro, .g96 or .pdb file.
        Default is False.
    verbose : bool
        full gromacs output (True)
    gmxexec : str
        gmx with the full path
    other_flags : str, optional
        additional flags to pass as you would type them in the shell

    Returns
    -------
    None

    """

    if gmxexec is None:
        gmxexec = get_gmx()

    if sep is True:
        other_flags += ' -sep'

    if verbose:
        capture_output = False
    else:
        capture_output = True

    if fit == 'none':
        subprocess.run('echo "{out_grp}" | {gmxexec} trjconv -f {f} -s {s} -o {o} -ur {ur} -pbc {pbc}'
             '{other_flags}'.format(gmxexec=gmxexec, f=f, s=s, o=o, ur=ur, pbc=pbc, out_grp=out_grp, other_flags=other_flags),
             shell=True, check=True, capture_output=capture_output)
    else:
        
        subprocess.run('echo "{fit_grp}" "{out_grp}" | {gmxexec} trjconv -f {f} -s {s} -o {o} -ur {ur} -pbc {pbc} -fit {fit}'
             '{other_flags}'.format(gmxexec=gmxexec, f=f, s=s, o=o, ur=ur, pbc=pbc, fit=fit,
                                   out_grp=out_grp, fit_grp=fit_grp, other_flags=other_flags), 
             shell=True, check=True, capture_output=capture_output)


def mdrun(s, deffnm='md', other_flags='', verbose=True, gmxexec=None):
    """Simple ``gmx mdrun`` wrapper.

    Parameters
    ----------
    s : str
        input tpr file
    deffnm : str, optional
        set the default filename for all file options. Default is 'md'.
    verbose : bool, optional
        whether to activate verbose flag in Gromacs mdrun. Default is False.
    gmxexec : str
        gmx with the full path
    other_flags : str, optional
        additional flags to pass as you would type them in the shell.

    Returns
    -------
    None

    """

    if gmxexec is None:
        gmxexec = get_gmx()

    if verbose is True:
        other_flags += ' -v'

    if verbose:
        capture_output = False
    else:
        capture_output = True

    subprocess.run('{gmxexec} mdrun -s {s} -deffnm {deffnm} {other_flags}'.format(gmxexec=gmxexec, s=s, deffnm=deffnm, other_flags=other_flags),
                   shell=True, check=True, capture_output=capture_output)


def write_mdp(mdp, fout='mdpfile.mdp', nsteps=10000, cutoff=1.0, T=300):
    """Writes a few standard mdp files.

    With the argument ``mdp`` you can choose from a few standard predefined
    mdp file: "enmin" (energy minimisation); "npt" (simulation in NPT ensemble);
    "npt-restr" (simulation in NPT ensemble with -DPOSRES defined).

    Parameters
    ----------
    mdp : str
        what type of mdp file to load and write.
        Options available are: 'enmin', 'npt-restr', 'npt'.
    fout : str
        filename of the mdp file to be written.
    nsteps : int
        number of steps.
    cutoff : float
        short-range cutoff (nm) for vdw and coulomb interactions.
    T : float
        temperature in Kelvin.

    Returns
    -------
    None
    """

    lines = mdps[mdp]
    with open(fout, 'w') as f:
        f.write(lines.format(nsteps=nsteps, cutoff=cutoff, T=T))
