#!/bin/bash

# Usage:
# ./ABFE_cli.sh --input_file [INPUT_FILE] --mdp_path [MDP_PATH] --cpu_num [CPU_NUM] --workdir [WORKDIR] 
#               --ligand_name [LIGAND_NAME] --lig [LIG_CASE] --apo [APO_CASE] --replicas [REPLICAS]

# Read the argument values
while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -i|--input_file) INPUT_FILE=$( realpath $2 ); shift;; # The input complex file in PDB format
      -m|--mdp_path) MDP_PATH=$( realpath $2 ); shift;;
      -n|--cpu_num) CPU_NUM=$2; shift;; # Number of available CPUs
      -w|--workdir) WORKDIR=$( realpath $2 ); shift;; # Working directory. It must be empty
      -l|--ligand_name) LIGAND_NAME=$2; shift;; # The name of the ligand in the input PDB file. Defaults to "LIG"
      -g|--lig) LIG_CASE=$2; shift;; # Ligand for the AbsoluteDG function. Defaults to "ligand"
      -a|--apo) APO_CASE=$2; shift;; # ApoCase for the AbsoluteDG function. Defaults to "apo"
      -r|--replicas) REPLICAS=$2; shift;; # The number of replicas to use for the calculations. Defaults to 1
    esac
    shift
done

# Check the argument values
if [ -z $INPUT_FILE ]
then
    echo "ERROR: --input_file argument is missed!"
    exit
else
    if [ ! -f "$INPUT_FILE" ]; then
        echo "ERROR: The input file does not exist!"
        exit
    fi
fi
if [ -z $MDP_PATH ]
then
    MDP_PATH=$( pwd )/mdps
fi
if [ -z $WORKDIR ]
then
    echo "ERROR: --workdir argument is missed!"
    exit
fi
if [ -z $CPU_NUM ]
then
    CPU_NUM=16
fi
if [ -z $LIGAND_NAME ]
then
    LIGAND_NAME="LIG"
fi
if [ -z $LIG_CASE ]
then
    LIG_CASE="ligand"
fi
if [ -z $APO_CASE ]
then
    APO_CASE="apo"
fi
if [ -z $REPLICAS ]
then
    REPLICAS=1
fi

# Run preparation
APP_DIR=$( realpath $( dirname -- "${BASH_SOURCE[0]}" ) )

chmod +x ${APP_DIR}/scripts/prep1.bash
${APP_DIR}/scripts/prep1.bash $INPUT_FILE $MDP_PATH $WORKDIR $CPU_NUM

FILE_PREFIX=$( basename $INPUT_FILE )
FILE_PREFIX=${FILE_PREFIX%.pdb}
python ${APP_DIR}/scripts/watertool.py ${WORKDIR}/${FILE_PREFIX}_cs.gro ${WORKDIR}/${FILE_PREFIX}_complex_w.gro ${WORKDIR}/${FILE_PREFIX}_complex.top

chmod +x ${APP_DIR}/scripts/prep2.bash
${APP_DIR}/scripts/prep2.bash $MDP_PATH ${WORKDIR}/${FILE_PREFIX} $CPU_NUM

chmod +x ${APP_DIR}/scripts/eq.bash
${APP_DIR}/scripts/eq.bash $MDP_PATH ${WORKDIR}/${FILE_PREFIX} $CPU_NUM

sed -n -e '/^\s*\[\s*system\s*\]\s*$/q' -e '/^\s*\[\s*moleculetype\s*\]\s*$/,$p' ${WORKDIR}/${FILE_PREFIX}_protein.top > ${WORKDIR}/${FILE_PREFIX}_protein.itp

# Create directory structure
for (( i = 1 ; i <= $REPLICAS ; i++))
do
    mkdir ${WORKDIR}/replicas_${i}
    mkdir ${WORKDIR}/replicas_${i}/struct_top
    mkdir ${WORKDIR}/replicas_${i}/struct_top/${LIG_CASE}
    mkdir ${WORKDIR}/replicas_${i}/struct_top/${APO_CASE}
    mkdir ${WORKDIR}/replicas_${i}/struct_top/${LIG_CASE}/protein
    mkdir ${WORKDIR}/replicas_${i}/struct_top/${LIG_CASE}/ligand
    mkdir ${WORKDIR}/replicas_${i}/struct_top/${APO_CASE}/protein
    mkdir ${WORKDIR}/replicas_${i}/eq_files

    # Copying files to the created directories
    cp ${WORKDIR}/${FILE_PREFIX}_protein.pdb ${WORKDIR}/replicas_${i}/struct_top/${APO_CASE}/protein/
    cp ${WORKDIR}/${FILE_PREFIX}_protein.itp ${WORKDIR}/replicas_${i}/struct_top/${APO_CASE}/protein/
    cp ${WORKDIR}/${FILE_PREFIX}_posre.itp ${WORKDIR}/replicas_${i}/struct_top/${APO_CASE}/protein/

    cp ${WORKDIR}/${FILE_PREFIX}_protein.pdb ${WORKDIR}/replicas_${i}/struct_top/${LIG_CASE}/protein/
    cp ${WORKDIR}/${FILE_PREFIX}_protein.itp ${WORKDIR}/replicas_${i}/struct_top/${LIG_CASE}/protein/
    cp ${WORKDIR}/${FILE_PREFIX}_posre.itp ${WORKDIR}/replicas_${i}/struct_top/${LIG_CASE}/protein/

    cp ${WORKDIR}/${FILE_PREFIX}_MOL.itp ${WORKDIR}/replicas_${i}/struct_top/${LIG_CASE}/ligand/
    cp ${WORKDIR}/${FILE_PREFIX}_ffMOL.itp ${WORKDIR}/replicas_${i}/struct_top/${LIG_CASE}/ligand/ffMOL.itp
    cp ${WORKDIR}/${FILE_PREFIX}_ligand.pdb ${WORKDIR}/replicas_${i}/struct_top/${LIG_CASE}/ligand/
    cp ${WORKDIR}/${FILE_PREFIX}_posre_ligand.itp ${WORKDIR}/replicas_${i}/struct_top/${LIG_CASE}/ligand/

    cp ${WORKDIR}/${FILE_PREFIX}_eq.tpr ${WORKDIR}/replicas_${i}/eq_files/
    cp ${WORKDIR}/${FILE_PREFIX}_eq.xtc ${WORKDIR}/replicas_${i}/eq_files/
    cp ${WORKDIR}/${FILE_PREFIX}_eq.xtc ${WORKDIR}/replicas_${i}/traj.xtc
    cp ${WORKDIR}/${FILE_PREFIX}_complex.top ${WORKDIR}/replicas_${i}/complex.top
    
done

#rm ${WORKDIR}/*
cp ${WORKDIR}/replicas_1/struct_top/${LIG_CASE}/protein/${FILE_PREFIX}_posre.itp ${WORKDIR}/${FILE_PREFIX}_posre.itp

# Run ABFE_cli.py script
export GMXLIB="/home/jovyan/pmx/src/pmx/data/mutff"
for (( i = 1 ; i <= $REPLICAS ; i++))
do
    python ${APP_DIR}/ABFE_cli.py --file_prefix $FILE_PREFIX --workdir ${WORKDIR}/replicas_${i} --ligand_name $LIGAND_NAME --lig $LIG_CASE --apoCase $APO_CASE --processes $CPU_NUM --replicas 1 &
done
wait
