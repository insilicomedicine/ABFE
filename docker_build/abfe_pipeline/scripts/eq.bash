#!/bin/bash

set -e -x

# Usage:
# ./eq.bash MDP_PATH WORKDIR_POSTFIX CPU_NUM

MDPPATH=$( realpath $1 )
WORKDIR_POSTFIX=$( realpath $2 )
CPU_NUM=$3

GRO=$( which gmx || which gmx_mpi )

sed -i 's/1000.*1000.*1000/ FC  FC  FC /g' ${WORKDIR_POSTFIX}_posre.itp
sed -i 's/^   FC  FC  FC /  1000     1   FC  FC  FC /g' ${WORKDIR_POSTFIX}_posre.itp
$GRO grompp -f ${MDPPATH}/nvt.mdp -p ${WORKDIR_POSTFIX}_complex.top -c ${WORKDIR_POSTFIX}_em_w.gro  -o ${WORKDIR_POSTFIX}_nvt.tpr -r ${WORKDIR_POSTFIX}_em_w.gro
$GRO mdrun -s ${WORKDIR_POSTFIX}_nvt.tpr -nt $CPU_NUM -cpt 30 -v -deffnm ${WORKDIR_POSTFIX}_nvt

$GRO grompp -f ${MDPPATH}/npt0.mdp -p ${WORKDIR_POSTFIX}_complex.top -c ${WORKDIR_POSTFIX}_nvt.gro  -o ${WORKDIR_POSTFIX}_npt0.tpr -r ${WORKDIR_POSTFIX}_em_w.gro -maxwarn 1
$GRO mdrun -s ${WORKDIR_POSTFIX}_npt0.tpr -nt $CPU_NUM -cpt 30 -v -deffnm ${WORKDIR_POSTFIX}_npt0

for i in {1..7}; do
    $GRO grompp -f ${MDPPATH}/npt${i}.mdp -p ${WORKDIR_POSTFIX}_complex.top -c ${WORKDIR_POSTFIX}_npt$(( i-1 )).gro -t ${WORKDIR_POSTFIX}_npt$(( i-1 )).cpt -o ${WORKDIR_POSTFIX}_npt$i.tpr -r ${WORKDIR_POSTFIX}_em_w.gro -maxwarn 1
    $GRO mdrun -s ${WORKDIR_POSTFIX}_npt$i.tpr -nt $CPU_NUM -cpt 30 -v -deffnm ${WORKDIR_POSTFIX}_npt$i
done
$GRO grompp -f ${MDPPATH}/eq.mdp -p ${WORKDIR_POSTFIX}_complex.top -c ${WORKDIR_POSTFIX}_npt7.gro -o ${WORKDIR_POSTFIX}_eq.tpr
$GRO mdrun -s ${WORKDIR_POSTFIX}_eq.tpr -nt $CPU_NUM -cpt 30 -v -deffnm ${WORKDIR_POSTFIX}_eq
