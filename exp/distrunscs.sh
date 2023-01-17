#!/bin/bash
dry="-"
nprocs="16 32 64 128 256"
nprocs="64"
cases=
trace="-"
op="potrf"
#que="workq"; 
timelimit="06:00:00"
#que="debug"
if [ $# -eq 1 -o $# -eq 2 ]; then
    cases=$1
    if [ ! -f "$cases" ];then
        echo "Cases file does not exist: $cases" 
        exit
    else
        . $cases
    fi
    if [ $# -eq 2 ]; then
        dry=$2
    fi
else 
    echo "Usage: $0 casesFile [dry]" 1>&2
    exit
fi
cdir=`basename $PWD`
if [[ $cdir != "hicma-dev" ]]; then
    echo "run this script in hicma-dev folder. The name of the current folder is: $cdir"
    exit -1
fi
# numthreads=39 #skylake
numthreads=31 #shaheen
numthreads="39" #cascadlake
# numthreads="1 2 4 8 16 31" #shaheen
hn="$HOSTNAME"
if [[ "$hn" = xci* ]]; then
    numthreads="16 32 60 61 62 63 64" #isambard
    numthreads="63" #isambard
    #numthreads="32 63" #isambard
fi
if [[ "$hn" = flamingo ]]; then
    numthreads="55" #flamingo
fi
if [[ "$hn" = vulture ]]; then
    numthreads="40" #vulture
fi
if [[ "$hn" = tuwaiq ]]; then
    numthreads="63" #tuwaiq
fi
if [[ "$hn" = shihab ]]; then
    numthreads="35" #shihab
fi
if [[ "$hn" = jasmine ]]; then
    numthreads="27" #jasmine
fi
if [[ "$hn" = kw60319 ]]; then
    numthreads="27" #kw60319
fi

#que="debug"; timelimit="00:30:00"; echo "DEBUG DEBUG DEBUG DEBUG DEBUG "
for nodes in $nprocs; do
    echo "#Number of nodes: $nodes ============================="
    hn="$HOSTNAME"
    if [[ "$hn" = xci* ]]; then
        que="arm"; 
        #cmdbatch="aprun -n $nodes -d 64 -j 1 " 
        cmdbatch="qsub -q $que -l select=$nodes -l walltime=$timelimit -j oe -N hicma-$nodes "
    elif [[ "$hn" = flamingo ]]; then
        cmdbatch=""
    elif [[ "$hn" = vulture ]]; then
        cmdbatch=""
    elif [[ "$hn" = tuwaiq ]]; then
        cmdbatch=""
    elif [[ "$hn" = shihab ]]; then
        cmdbatch=""
    elif [[ "$hn" = jasmine ]]; then
        cmdbatch=""
    elif [[ "$hn" = kw60319 ]]; then
        cmdbatch=""
    else
        cmdbatch="sbatch --parsable --nodes=$nodes --time=$timelimit  --job-name=hicma-$nodes "  #shaheen
    fi
    echo "# $cmdbatch"
    if [ "$dry" == "dry" ]; then
        cmdbatch=""
    fi
    maxsubs[1]=500; minsubs[1]=250
    maxsubs[2]=1000; minsubs[2]=500
    maxsubs[3]=2000; minsubs[3]=1000
    maxsubs[4]=4000; minsubs[4]=2000
    maxsubs[5]=6000; minsubs[5]=3000
    #for isub in 1 2 3 4 5;do
    for isub in 0;do
        #maxsub=${maxsubs[isub]}
        #minsub=${minsubs[isub]}
        sched="prio eager random ws lws dm dmda dmdar dmdas dmdasd heft pheft peager" 
        sched="prio"

        #cases="exp/cases/cham32.txt";       startt=107; endt=108; step=200
        #ids="`seq 1 5`"
        #$cmdbatch exp/trial.sh $nodes $nodes 31 - cham - "$ids" $dry $cases "13:00:00" $que $sched $maxsub $minsub;continue


        caseids=${allcaseids[$nodes]}      
        ncases=`echo "$caseids" | wc -w`    
        startt=0;   endt=$((ncases-1)); 
    #    echo -n "#`date` on $nodes nodes. $note - $prog - $cases "
    #    echo -n \"$caseids\" 
    #    echo " $ncases ($startt-$endt-$step)"

        ct=$startt
        while [ $ct -le $endt ]; do
            et=$((ct+step-1))
            if [ $et -gt $endt ]; then
                et=$endt
            fi
            #ids="`seq $ct $et`"
            arrcaseids=($caseids)
            ids=
            for icase in `seq $ct $et`; do
                ids="$ids ${arrcaseids[icase]}"
            done
            ct=$((ct+step))
    #        echo "#case ids: $ids"
            for nt in $numthreads; do
                #if [[ "$hn" = xci* ]]; then
                #    cmdbatch2="$cmdbatch -o $PWD/exp/out/$prog-$nodes-$nt-$PBS_JOBID"
                #fi
                corespernode=$nt
                maxsub=$((nodes*(corespernode+10)))
                minsub=$((nodes*corespernode))
                cmd="$PWD/exp/distmemcs.sh $nodes $nodes $nt $trace $prog - \"$ids\" $dry $PWD/$cases $timelimit  $sched $maxsub $minsub $op" 
                if [ "$dry" == "dry" ]; then
                    eval $cmd 
                    continue
                fi
                if [[ "$hn" = xci* ]]; then
                    cmd="$PWD/exp/distmemcs.sh  $nodes $nodes $nt $trace $prog - \"$ids\" $dry $PWD/$cases $timelimit $sched $maxsub $minsub $op" 
                    if [ "$dry" == "dry" ]; then
                        eval $cmd 
                        continue
                    fi
                    jobid=$(echo $cmd | $cmdbatch) 
                    echo "$jobid"
                    date=$(date '+%Y-%m-%d--%H-%M')
                    qalter -o $PWD/exp/out/$prog-$nodes-$nt--$date--$jobid $jobid
                elif [[ "$hn" = flamingo ]]; then
                    eval $cmd
                elif [[ "$hn" = shihab ]]; then
                    eval $cmd
                elif [[ "$hn" = jasmine ]]; then
                    eval $cmd
                elif [[ "$hn" = kw60319 ]]; then
                    eval $cmd
                else
                    $cmdbatch $PWD/exp/distmemcs.sh $nodes $nodes $nt $trace $prog - "$ids" $dry $cases $timelimit  $sched $maxsub $minsub $op
                fi
            done
        done
    done
done

