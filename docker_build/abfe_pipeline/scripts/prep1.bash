#!/bin/bash
set -e -x

# Usage:
# ./prep1.bash PDB_FILE_PATH MDP_PATH WORKDIR_PATH CPU_NUM

GPU=0

GRO=$( which gmx|| which gmx_mpi )

PDB_FILE=$( realpath $1 )
MDPPATH=$( realpath $2 )
WORKDIR=$( realpath $3 )
CPU_NUM=$4

BASENAME=$( basename $PDB_FILE )
BASENAME=${BASENAME%.pdb}

PROTEIN_PDB=${WORKDIR}/${BASENAME}_protein.pdb
LIGAND_PDB=${WORKDIR}/${BASENAME}_ligand.pdb
COMPLEX_PDB=${WORKDIR}/${BASENAME}_complex.pdb

FILE_PREFIX=$( basename $LIGAND_PDB )
FILE_PREFIX=${FILE_PREFIX%_ligand.pdb}

# Separate protein and ligand from the complex PDB
grep -v "^HETATM" $PDB_FILE > $PROTEIN_PDB
grep "^HETATM" $PDB_FILE > $LIGAND_PDB

# Add hydrogens to the ligand
#obabel -ipdb $LIGAND_PDB -O $LIGAND_PDB -h

# Generate ligand topology and coordinates using ACPYPE
export OMP_NUM_THREADS=12
acpype -i $LIGAND_PDB
unset OMP_NUM_THREADS

# Rename output and clean useless files
cp $( pwd )/${FILE_PREFIX}_ligand.acpype/${FILE_PREFIX}_ligand_GMX.gro ${WORKDIR}/${FILE_PREFIX}_ligand.gro
cp $( pwd )/${FILE_PREFIX}_ligand.acpype/${FILE_PREFIX}_ligand_GMX.itp ${WORKDIR}/${FILE_PREFIX}_MOL.itp
cp $( pwd )/${FILE_PREFIX}_ligand.acpype/posre_${FILE_PREFIX}_ligand.itp ${WORKDIR}/${FILE_PREFIX}_posre_ligand.itp
#rm -r $( pwd )/${FILE_PREFIX}_ligand.acpype
# Remove atomtypes, since in amber force fields they are shared
sed -n '/^\[ atomtypes \]/,/^$/p' ${WORKDIR}/${FILE_PREFIX}_MOL.itp > ${WORKDIR}/${FILE_PREFIX}_ffMOL.itp
cp ${WORKDIR}/${FILE_PREFIX}_MOL.itp ${WORKDIR}/${FILE_PREFIX}_MOL_pmx.itp
sed -i '/^\[ atomtypes \]/,/^$/d' ${WORKDIR}/${FILE_PREFIX}_MOL.itp


# $GRO pdb2gmx -f "$PROTEIN_PDB" -ff amber99sb-ildn -o protein.gro -p protein.top -water tip3p -merge all -renum -ignh
# $GRO editconf -f protein.gro -o protein.pdb
# $GRO editconf -f ligand.gro -o ligand.pdb

$GRO pdb2gmx -f $PROTEIN_PDB -ff amber99sb-ildn -o ${WORKDIR}/${FILE_PREFIX}_protein.gro -p ${WORKDIR}/${FILE_PREFIX}_protein.top -water tip3p -merge all -renum -ignh -i ${WORKDIR}/${FILE_PREFIX}_posre.itp
$GRO editconf -f ${WORKDIR}/${FILE_PREFIX}_protein.gro -o ${WORKDIR}/${FILE_PREFIX}_protein.pdb
$GRO editconf -f ${WORKDIR}/${FILE_PREFIX}_ligand.gro -o ${WORKDIR}/${FILE_PREFIX}_ligand.pdb

APP_DIR=$( realpath $( dirname -- "${BASH_SOURCE[0]}" ) )
if test -f ${APP_DIR}/mdout.mdp; then
    rm ${APP_DIR}/mdout.mdp
fi

# Combine protein and ligand topologies
awk -v ffMOL_itp_file=${WORKDIR}/${FILE_PREFIX}_ffMOL.itp '
  {
    if ($0 ~ /#include ".*forcefield.itp"/) {
      print;
      print "; Include ligand parameters";
      print "#include \""ffMOL_itp_file"\"";
    } else {
      print;
    }
  }
  ' ${WORKDIR}/${FILE_PREFIX}_protein.top > ${WORKDIR}/${FILE_PREFIX}_protein_tmp.top


python $APP_DIR/correct_topol.py ${WORKDIR}/${FILE_PREFIX}_protein_tmp.top ${WORKDIR}/${FILE_PREFIX}_posre.itp ${WORKDIR}/${FILE_PREFIX}_MOL.itp

mv ${WORKDIR}/${FILE_PREFIX}_protein_tmp.top ${WORKDIR}/${FILE_PREFIX}_complex.top

sed -i "/; Compound[[:space:]]*#mols/a ${FILE_PREFIX}_ligand  1" ${WORKDIR}/${FILE_PREFIX}_complex.top

# Combine protein and ligand .gro files
head -n 1 ${WORKDIR}/${FILE_PREFIX}_protein.gro > ${WORKDIR}/${FILE_PREFIX}_complex.gro
x1=$( grep -c ^ ${WORKDIR}/${FILE_PREFIX}_protein.gro )
x2=$( grep -c ^ ${WORKDIR}/${FILE_PREFIX}_ligand.gro )
echo $x1 $x2 | awk '{print $1+$2-6}' >> ${WORKDIR}/${FILE_PREFIX}_complex.gro
tail -n +3 ${WORKDIR}/${FILE_PREFIX}_ligand.gro | head -n -1 >> ${WORKDIR}/${FILE_PREFIX}_complex.gro # Swap the order
tail -n +3 ${WORKDIR}/${FILE_PREFIX}_protein.gro | head -n -1 >> ${WORKDIR}/${FILE_PREFIX}_complex.gro # Swap the order
tail -n 1 ${WORKDIR}/${FILE_PREFIX}_protein.gro >> ${WORKDIR}/${FILE_PREFIX}_complex.gro

echo 0 | $GRO editconf -f ${WORKDIR}/${FILE_PREFIX}_complex.gro  -o ${WORKDIR}/${FILE_PREFIX}_complex_c.gro -d 1 -princ -bt dodecahedron

$GRO grompp -f ${MDPPATH}/em.mdp -p ${WORKDIR}/${FILE_PREFIX}_complex.top  -c ${WORKDIR}/${FILE_PREFIX}_complex_c.gro -o ${WORKDIR}/${FILE_PREFIX}_em.tpr -maxwarn 1

$GRO mdrun -s ${WORKDIR}/${FILE_PREFIX}_em.tpr -nt $CPU_NUM -deffnm ${WORKDIR}/${FILE_PREFIX}_em -v 

$GRO solvate -cs spc216 -cp ${WORKDIR}/${FILE_PREFIX}_em.gro -p ${WORKDIR}/${FILE_PREFIX}_complex.top -o ${WORKDIR}/${FILE_PREFIX}_complex_s.gro
$GRO editconf -f ${WORKDIR}/${FILE_PREFIX}_complex_s.gro -o ${WORKDIR}/${FILE_PREFIX}_cs.gro -resnr 1
