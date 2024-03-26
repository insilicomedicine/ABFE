#!/bin/bash
set -e -x

# Usage:
# ./prep2.bash MDP_PATH WORKDIR_POSTFIX CPU_NUM

MDPPATH=$( realpath $1 )
WORKDIR_POSTFIX=$( realpath $2 )
CPU_NUM=$3

GRO=$( which gmx || which gmx_mpi )

$GRO grompp -f ${MDPPATH}/em.mdp -p ${WORKDIR_POSTFIX}_complex.top  -c ${WORKDIR_POSTFIX}_complex_w.gro -o ${WORKDIR_POSTFIX}_em_s.tpr -maxwarn 1
echo 'SOL' | $GRO genion -s ${WORKDIR_POSTFIX}_em_s.tpr -p ${WORKDIR_POSTFIX}_complex.top  -o ${WORKDIR_POSTFIX}_em_si.gro -neutral -conc 0.15

$GRO grompp -f ${MDPPATH}/em.mdp -p ${WORKDIR_POSTFIX}_complex.top -c ${WORKDIR_POSTFIX}_em_si.gro  -o ${WORKDIR_POSTFIX}_em_w.tpr
$GRO mdrun -s ${WORKDIR_POSTFIX}_em_w.tpr -nt $CPU_NUM -deffnm ${WORKDIR_POSTFIX}_em_w -c 
