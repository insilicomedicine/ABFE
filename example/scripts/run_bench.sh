n=0
for i in 3UNJ.pdb 3S00.pdb 3R8V.pdb 3PIY.pdb; do
        k=${i%.*}
        cp -r /home/jovyan/work /home/jovyan/data/work_${k}
        mkdir -p /home/jovyan/data/work_${k}/output
        cd /home/jovyan/data/work_${k}/ 
        export GMX_MAXBACKUP=500
        export CUDA_VISIBLE_DEVICES=$n
        /home/jovyan/data/work_${k}/ABFE_cli.sh -i /home/jovyan/data/pdbs/${i} -n 16 -w /home/jovyan/data/work_${k}/output/ -l LIG -r 1 > /home/jovyan/data/work_${k}/log 2>&1 &
        n=$((n+1))
done
wait < <(jobs -p)
