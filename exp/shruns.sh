#!/bin/bash
dry="-"
if [ $# -eq 1 ]; then
    dry=$1
fi


nprocs=" 16 32 64 128 256"
nprocs="16"
nprocs="32"
nprocs="64 128 256" ##TODO
nprocs="128 256 512" ##TODO
nprocs="128"
#nprocs="512 1024"

#que="workq"; timelimit="23:50:00"
que="workq"; timelimit="04:00:00"
#que="workq"; timelimit="01:00:00"
#que="workq"; timelimit="00:20:00"
#que="debug"; timelimit="00:30:00"; echo "DEBUG DEBUG DEBUG DEBUG DEBUG "
#prio eager random ws lws dm dmda dmdar dmdas dmdasd heft pheft peager
for nodes in $nprocs; do
    if [ "$dry" == "dry" ]; then
        cmdbatch="echo "
        cmdbatch=""
    else
        cmdbatch=""
    fi
    maxsubs[1]=500; minsubs[1]=250
    maxsubs[2]=1000; minsubs[2]=500
    maxsubs[3]=2000; minsubs[3]=1000
    maxsubs[4]=4000; minsubs[4]=2000
    maxsubs[5]=6000; minsubs[5]=3000
#    for isub in 1 2 3 4 5;do
    for isub in 1;do
        maxsub=${maxsubs[isub]}
        minsub=${minsubs[isub]}
        echo "#Job for max/min sub:$maxsub $minsub on $nodes nodes"
        #$cmdbatch exp/trial.sh $nodes $nodes 31 - hic - "`seq 1 1`" $dry exp/cases/biggersize-allnb-2.txt "13:00:00" workq "prio" $maxsub $minsub
        endt=434
        #endt=1
        sched="prio eager random ws lws dm dmda dmdar dmdas dmdasd heft pheft peager" 
        sched="prio"
        #startt=1; endt=182
        #cases="exp/cases/mpi32.txt"
        #startt=1; endt=182
        #cases="exp/cases/1M.txt"
        
        cases="exp/cases/1M-few-case.txt"
        cases="exp/cases/shmem-1.txt"
        cases="exp/cases/shmem-2.txt"
        startt=1; endt=84
        startt=1; endt=1 #TODO
        step=$endt
        ct=$startt
        while [ $ct -le $endt ]; do
            et=$((ct+step-1))
            if [ $et -gt $endt ]; then
                et=$endt
            fi
            ids="`seq $ct $et`" 
            ct=$((ct+step))
            #echo $ids
            #$cmdbatch exp/shmem.sh $nodes $nodes 56 - hic - "$ids" $dry $cases "13:00:00" $que $sched $maxsub $minsub
            #$cmdbatch exp/shmem.sh $nodes $nodes 112 - hic - "$ids" $dry $cases "13:00:00" $que $sched $maxsub $minsub
        nts="9 18 36 72" #shihab haswell
        nts="36" #shihab haswell
        nts="28" #jasmine broadwell

        nts="16" #oqab sandy bridge
        nts="56" #skylake

	ids="84"
        #nts="11" #uwork
	    for nt in $nts; do
            	$cmdbatch $HOME/hicma/exp/shmem.sh $nodes $nodes $nt - hic - "$ids" $dry $HOME/hicma/$cases "13:00:00" $que $sched $maxsub $minsub
	    done 
        done
    done
done




    #sbatch --nodes=$nodes exp/trial.sh $nodes $nodes 31 - hic - "`seq 1 125`" -
    #sbatch --nodes=$nodes exp/trial.sh $nodes $nodes 31 - hic minmaxsub "`seq 1 125`" -
    #sbatch --nodes=$nodes exp/trial.sh $nodes $nodes 31 - hic minmaxsub "`seq 1 5`" - exp/cases/t2.txt
    #sbatch --nodes=$nodes exp/trial.sh $nodes $nodes 31 - hic minmaxsub "`seq 1 21`" - exp/cases/t2.txt
    #sbatch --nodes=$nodes exp/trial.sh $nodes $nodes 31 - hic minmaxsub "`seq 70 125`" - exp/cases/allsize-allnb-1.txt "13:00:00" workq
    #sbatch --nodes=$nodes exp/trial.sh $nodes $nodes 31 - hic -         "`seq 70 125`" - exp/cases/allsize-allnb-1.txt "13:00:00" workq
    #sbatch --nodes=$nodes exp/trial.sh $nodes $nodes 31 - hic minmaxsub "`seq 1 125`" - exp/cases/allsize-allnb-1.txt "13:00:00" workq
    #sbatch --nodes=$nodes exp/trial.sh $nodes $nodes 31 - hic minmaxsub "`seq 1 61`" - exp/cases/biggersize-allnb-1.txt "13:00:00" workq

    #bash exp/trial.sh $nodes $nodes 31 - hic minmaxsub "1" dry exp/cases/biggersize-allnb-1.txt "13:00:00" workq "prio"
    
    #sbatch --nodes=$nodes exp/trial.sh $nodes $nodes 31 "1  11 12  21  31  41 51  61  71 81"
    #sbatch --nodes=$nodes exp/trial.sh $nodes $nodes 31 "0 1 10 11 12 20 21 30 31 40 41 42 50 51 60 61 70 71 72 80 81"
    #bash exp/trial.sh $nodes $nodes 31 "0 1 10 11 12 20 21 30 31 40 41 42 50 51 60 61 70 71 72 80 81"
    #sbatch --nodes=$nodes exp/trial.sh $nodes
