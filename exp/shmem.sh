#!/bin/bash
#SBATCH --job-name=hicma
#SBATCH --account=k1205
#SBATCH --output=/project/k1205/akbudak/hicma/exp/out/%j
#SBATCH --error=/project/k1205/akbudak/hicma/exp/err/%j
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kadir.akbudak@kaust.edu.sa
#SBATCH --threads-per-core=1
#SBATCH --hint=nomultithread
##SBATCH --time=13:00:00
##SBATCH --time=23:00:00
#SBATCH --time=00:01:00


#3teraflops
#export STARPU_MIN_PRIO=-1
#export STARPU_MAX_PRIO=250

#These do not affect perf for dpotrf M=36000 mb=1000
export STARPU_CALIBRATE=0
#export STARPU_TARGET_AVAILABLE_MEM=8
#export STARPU_MINIMUM_AVAILABLE_MEM=0
#export STARPU_DISK_SWAP_SIZE=0
#export STARPU_LIMIT_CPU_MEM=135000
#export STARPU_LIMIT_MAX_SUBMITTED_TASKS=20000
#export STARPU_LIMIT_MIN_SUBMITTED_TASKS=18000


#FOR DEBUGGING
#export STARPU_WATCHDOG_CRASH=1
#export STARPU_WATCHDOG_TIMEOUT=100000000 #100s
#export STARPU_GENERATE_TRACE=1
#export STARPU_MEMORY_STATS=1
#export STARPU_MAX_MEMORY_USE=1
#export STARPU_BUS_STATS=1
#export STARPU_WORKER_STATS=1
#export STARPU_STATS=1
#export STARPU_IDLE_TIME=/project/k1205/akbudak/hicma/exp/starpulog/$SLURM_JOBID.idle
#export STARPU_PROFILING=1
#STARPU_SILENT=1  

#echo STARPU_SCHED                     $STARPU_SCHED
#echo STARPU_CALIBRATE                 $STARPU_CALIBRATE
#echo STARPU_TARGET_AVAILABLE_MEM      $STARPU_TARGET_AVAILABLE_MEM 
#echo STARPU_MINIMUM_AVAILABLE_MEM     $STARPU_MINIMUM_AVAILABLE_MEM
#echo STARPU_LIMIT_CPU_MEM             $STARPU_LIMIT_CPU_MEM
#echo STARPU_LIMIT_MAX_SUBMITTED_TASKS $STARPU_LIMIT_MAX_SUBMITTED_TASKS
#echo STARPU_LIMIT_MIN_SUBMITTED_TASKS $STARPU_LIMIT_MIN_SUBMITTED_TASKS
#echo STARPU_FXT_TRACE                 $STARPU_FXT_TRACE
#echo STARPU_FXT_PREFIX                $STARPU_FXT_PREFIX
#echo STARPU_LOGFILENAME               $STARPU_LOGFILENAME

numnodes=1
numtasks=1
numthreads=1
tracearg=""
prog=""
minmaxsub=""
dry=""
if [ $# -eq 14 ]; then
    numnodes=$1
    nummpi=$2
    numthreads=$3
    trace=$4
    prog=$5
    minmaxsub=$6
    exps="$7"
    dry=$8
    sizefile=$9
    time_limit=${10}
    queue=${11}
    scheds=${12}
    maxsub=${13}
    minsub=${14}
else
    echo "usage: $0 numnodes nummpi numthreads [trace,-] [hic,cham] [enable minmaxsub depending on #tasks,-] exps [dry,-] sizefile time_limit queue scheds [maxsub,-] [minsub,]"
    echo
    echo "Your input:"
    echo $*
    exit
fi
#Calculate number of nodes and number of mpi tasks per node
sr=$numnodes #$(echo "scale=0;sqrt ( $numnodes ) / 1" | bc -l) ;
p=$(echo " ( l($sr) / l(2) )/2" | bc -l) ;
p=$(echo " scale=0;( $p / 1 ) " | bc -l) ;
sqrt_numnodes=$(echo " scale=0;( 2 ^ $p )/1 " | bc -l) ;
sqrt_numnodesQ=$((numnodes/sqrt_numnodes))
#echo $sqrt_numnodes $sqrt_numnodesQ; exit
ntasks_per_node=$((nummpi/numnodes))


#Parameters for HiCMA
_factor=2 #TODO #ratio of mb/maxrank
_acc=8          #fixed accuracy


#Sizes of matrices for experimental cases
if [ ! -f "$sizefile" ]; then
    echo "Size file does not exist: $sizefile"
    exit
fi
. $sizefile


#cham gives this error for prio:
#[starpu][_starpu_priority_push_task][assert failure] task priority 1464 is not between minimum -5 and maximum 5

sruncmd="" 

#Loop over scheduling algorithms of StarPU
for sched in $scheds;do
    export STARPU_SCHED=$sched
    #Loop over experimental cases
    for iexp in $exps;do 
        echo Experiment case:$iexp nrows:${nrows[iexp]} mb:${nb[iexp]}
        _m=${nrows[iexp]} 
        _b=${nb[iexp]}
        if [ -z "$_m" ]; then 
            continue
        fi
        if [ -z "$_b" ]; then 
            continue
        fi
        _mb=$_b;

        _storage_maxrk=200
        _n=$((_m/_mb*_storage_maxrk));
        _nb=$_storage_maxrk
	#ss
        _maxrank=$_nb
	#edsin
        _maxrank=400


        if [ "$prog" == "hic" ]; then
	    rankfile=$HOME/hicma/exp/ranks/$prog-$sched-$_m-$_mb-$_nb-$numnodes-$nummpi-$numthreads-$SLURM_JOBID-SS-$_acc-$_maxrank
            cmd="$HOME/hicma/build/timing/time_zpotrf_tile \
            --m=$_m \
            --n_range=$_n:$_n \
            --k=$_m \
            --mb=$_mb \
            --nb=$_nb \
            --nowarmup \
            --threads=$numthreads \
            --p=1 \
            $tracearg \
            --rk=0 \
            --acc=$_acc \
            --starshdecay=2 \
            --starshmaxrank=$_maxrank \
            --starshwavek=100 \
	    --ss \
	    --rankfile=$rankfile \
	    "
            
#--ss \
        elif [ "$prog" == "cham" ]; then
            if [ "$trace" != "-" ]; then
                cmd="/project/k1205/akbudak/hicma/chameleon/build-fxt-s121/timing/time_dpotrf_tile --nowarmup --p=$sqrt_numnodes --m=$_m --n_range=$_m:$_m:$_m --nb=$_b --threads=$numthreads $tracearg"
            else
                cmd="/project/k1205/akbudak/hicma/chameleon/build/timing/time_dpotrf_tile --nowarmup --p=$sqrt_numnodes --m=$_m --n_range=$_m:$_m:$_m --nb=$_b --threads=$numthreads"
            fi
        fi
        minmaxsubinfo=
        if [ "$minmaxsub" != "-" ]; then
            _mt=$((_m/_mb))
            _ntasks=$((_mt*_mt*_mt/3))
            _maxsub=$((_ntasks/8))
            _minsub=$((_ntasks/10))
            #_maxsub=10000; _minsub=8000
            #_maxsub=1000;  _minsub=500
            minmaxsubinfo="MT:$_mt NTASKS:$_ntasks MAXSUB:$_maxsub MINSUB:$_minsub"
            export STARPU_LIMIT_MAX_SUBMITTED_TASKS=$_maxsub
            export STARPU_LIMIT_MIN_SUBMITTED_TASKS=$_minsub
        fi
        if [ "$maxsub" != "-" -a "$minsub" != "-" ]; then
            minmaxsubinfo="MAXSUB:$maxsub MINSUB:$minsub"
            export STARPU_LIMIT_MAX_SUBMITTED_TASKS=$maxsub
            export STARPU_LIMIT_MIN_SUBMITTED_TASKS=$minsub
        fi

        msg="M:$_m N:$_n MB:$_mb NB:$_nb MAXRANK:$_maxrank DATE:`date` SCHED:$STARPU_SCHED CMD:$cmd $minmaxsubinfo" 
        echo "!BEGIN:" $msg 
        if [ "$trace" != "-" ]; then
            traceprefix=`pwd`/exp/trace/$prog-$sched-$_m-$_mb-$_nb-$numnodes-$nummpi-$numthreads/$SLURM_JOBID
            mkdir -p $traceprefix
            export STARPU_FXT_PREFIX=$traceprefix
        fi
        if [ "$dry" == "dry" ]; then
            echo $cmd
        else
            echo "!BEGIN:" $msg 1>&2 
	    export MKL_NUM_THREADS=1 
	    export KMP_AFFINITY=granularity=fine,scatter
            $sruncmd $cmd
            echo "!END:" $msg 
            echo "!END:" $msg 1>&2 
        fi
        if [ "$trace" != "-" ]; then
            combinedtrace=${traceprefix}trace
            /project/k1205/akbudak/codes/starpu-1.2.1-fxt/tools/starpu_fxt_tool -i ${traceprefix}prof_file_akbudak_* -o ${combinedtrace}
            mv activity.data  dag.dot  data.rec  distrib.data  tasks.rec  trace.html $traceprefix/ #paje.trace    
            echo $SLURM_JOBID > $traceprefix/0_jobid 
            #echo "Dot is starting:"
            #dot -Tpng $traceprefix/dag.dot -o $traceprefix/$prog.png
            #echo "Dot ended"
        fi
        echo
        date
    done
done

exit 0
    --printindex \
    --check \
    --printmat \
    --printindex \
    --trace \
--tag-output --timestamp-output 

$SLURM_JOBID    Job ID  5741192 $PBS_JOBID
$SLURM_JOB_NAME Job Name    myjob   $PBS_JOBNAME
$SLURM_SUBMIT_DIR   Submit Directory    /lustre/payerle/work    $PBS_O_WORKDIR
$SLURM_JOB_NODELIST Nodes assigned to job   compute-b24-[1-3,5-9],compute-b25-[1,4,8]   cat $PBS_NODEFILE
$SLURM_SUBMIT_HOST  Host submitted from login-1.deepthought2.umd.edu    $PBS_O_HOST
$SLURM_JOB_NUM_NODES    Number of nodes allocated to job    2   $PBS_NUM_NODES
$SLURM_CPUS_ON_NODE Number of cores/node    8,3 $PBS_NUM_PPN
$SLURM_NTASKS   Total number of cores for job???    11  $PBS_NP
$SLURM_NODEID   Index to node running on
relative to nodes assigned to job   0   $PBS_O_NODENUM
$PBS_O_VNODENUM Index to core running on
within node 4   $SLURM_LOCALID
$SLURM_PROCID   Index to task relative to job   0   $PBS_O_TASKNUM - 1

