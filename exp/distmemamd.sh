#!/bin/bash
###PBS -o /home/ri-kakbudak/hicma-dev/exp/out/
#SBATCH --job-name=hicma
##SBATCH --account=k1205
#SBATCH --output=/ibex/scratch/omairyrm/hicma-amd/hicma-dev/amd-exp/out/%j
#SBATCH --error=/ibex/scratch/omairyrm/hicma-amd/hicma-dev/amd-exp/err/%j
#SBATCH --cpus-per-task=128
#SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=rabab.alomairy@gmail.com
#SBATCH --partition=users
#SBATCH --threads-per-core=1
#SBATCH --hint=nomultithread


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
#export STARPU_MAX_MEMORY_USE=1
#export STARPU_IDLE_TIME=/project/k1205/akbudak/hicma/exp/starpulog/$SLURM_JOBID.idle
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
numthreads=31
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
#    queue=${11}
    scheds=${11}
    maxsub=${12}
    minsub=${13}
    op=${14}
else
    echo "usage: $0 numnodes nummpi numthreads [trace,-] [hic,cham] [enable minmaxsub depending on #tasks,-] exps [dry,-] sizefile time_limit queue scheds [maxsub,-] [minsub,] [potrf posv]"
    echo
    echo "Your input has $# arguments as follows:"
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


hn=$HOSTNAME
echo $HOSTNAME
HICMADIR=$HOME/hicma-amd/hicma-dev
BINDIR=$HICMADIR/build/timing
CHAMDIR=$HICMADIR/chameleon/build
if [[ "$hn" = xci* ]]; then #isambard
    sruncmd="aprun -n $numnodes -d 64 -j 1 " 
    cham_block_size_arg="--nb=300"
elif [[ "$hn" = flamingo ]]; then
    sruncmd="numactl --interleave=all  "
elif [[ "$hn" = vulture ]]; then
    sruncmd="numactl --interleave=all  "
elif [[ "$hn" = tuwaiq ]]; then
    sruncmd="numactl --interleave=all  "
elif [[ "$hn" = shihab ]]; then
    sruncmd="numactl --interleave=all  "
elif [[ "$hn" = jasmine ]]; then
    sruncmd="numactl --interleave=all  "
elif [[ "$hn" = kw60319 ]]; then
    sruncmd="numactl --interleave=all  "
else
    HICMADIR=/ibex/scratch/omairyrm/hicma-amd/hicma-dev
    sruncmd="\
        srun \
        --job-name=hicma-$_m-$_mb-$SLURM_JOB_NUM_NODES --hint=nomultithread \
        --nodes=$numnodes \
        --ntasks=$nummpi \
        --ntasks-per-node=$ntasks_per_node \
        numactl --interleave=all  \
        " 
    BINDIR=$HICMADIR/build/timing
    CHAMDIR=$HICMADIR/chameleon/build
fi
echo $sruncmd
echo


if [ "$op" == "potrf" ]; then
    BINNAME=time_z${op}_tile
elif [ "$op" == "posv" ]; then
    #BINDIR=/project/k1205/akbudak/hicma/build-only-prob/testing  #TODO
    BINNAME=testing_z${op}
elif [ "$op" == "genc" ]; then
    #BINDIR=/project/k1205/akbudak/hicma/build-only-prob/testing  #TODO
    BINNAME=time_z${op}_tile
else
    echo "Wrong op: $op"
    exit
fi
BIN=$BINDIR/$BINNAME

#Set LD_LIBRARY_PATH
if [ "$trace" != "-" ]; then
    tracearg=" --trace --profile --dag "
    export STARPU_FXT_TRACE=1
    export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64:/project/k1205/akbudak/codes/starpu13-install-cdt-fxt/lib/:/project/k1205/akbudak/hicma/chameleon/build-cdt-fxt-s13/lib/:/project/k1205/akbudak/codes/hwloc-1.11.7/install-cdt/lib/:/project/k1205/akbudak/codes/fxt-0.3.7/install/lib:$LD_LIBRARY_PATH
    BINDIR=/project/k1205/akbudak/hicma/build-cdt-fxt-s13/timing
    BIN=$BINDIR/$BINNAME
    ##extra flags to obtain more info
    export STARPU_GENERATE_TRACE=1
    export STARPU_MEMORY_STATS=1
    export STARPU_BUS_STATS=1
    export STARPU_WORKER_STATS=1
    export STARPU_STATS=1
    export STARPU_PROFILING=1
else
    #PERFORMANCE RUNS, WITHOUT FXT
    #export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64:/project/k1205/akbudak/hicma/chameleon/build/lib/:/project/k1205/akbudak/codes/hwloc-1.11.7/install/lib/:/project/k1205/akbudak/codes/starpu-1.2.2/install/lib/:$LD_LIBRARY_PATH
    if [ "$dry" != "dry" ]; then
        hn=$HOSTNAME
        if [[ "$hn" = xci* ]]; then #isambard
            . $HOME/hicma-dev/scripts/modules-isambard-allinea.sh
        else
            module load cdt
            #export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64:/project/k1205/toae/2018-01-26/hicma/chameleon/build/lib/:/project/k1205/akbudak/codes/hwloc-1.11.7/install-cdt/lib/:/project/k1205/akbudak/codes/starpu13-install-cdt/lib/:$LD_LIBRARY_PATH
            export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64:/project/k1205/akbudak/codes/hwloc-1.11.7/install-cdt/lib/:$LD_LIBRARY_PATH
        fi
    fi
fi

#Print linked libraries to stderr
if [ "$dry" != "dry" ]; then
    ldd $BIN 1>&2
    echo $LD_LIBRARY_PATH  1>&2
fi

#Parameters for HiCMA
_factor=4 #TODO #ratio of mb/maxrank
_acc=8          #fixed accuracy


#Sizes of matrices for experimental cases
if [ ! -f "$sizefile" ]; then
    echo "Size file does not exist: $sizefile"
    exit
fi

#DEFAULT PARAMETERS
_wavek=0
_decay=0
#CUSTOM PARAMETERS
. $sizefile


#cham gives this error for prio:
#[starpu][_starpu_priority_push_task][assert failure] task priority 1464 is not between minimum -5 and maximum 5

#Loop over scheduling algorithms of StarPU
for sched in $scheds;do
    export STARPU_SCHED=$sched
    #Loop over experimental cases
    for iexp in $exps;do 
        echo Experiment case:$iexp nrows:${nrows[iexp]} mb:${nb[iexp]}
        _m=${nrows[iexp]} 
        _b=${nb[iexp]}
        _maxrank=${maxrank[iexp]}
        _appdata=${appdata[iexp]}
        _nrhs=${nrhs[iexp]}
        _rad=${rad[iexp]}
        _reg=${reg[iexp]}
        _numobj=${numobj[iexp]}
        _order=${order[iexp]}
        _mesh_file=${mesh_file[iexp]}
        _rbf_kernel=${rbf_kernel[iexp]}
        for _acc in ${acc[iexp]}; do
          for _ker in $_rbf_kernel; do
            if [ ! -z "${decay[iexp]}" ]; then 
                _decay=${decay[iexp]}
            fi
            if [ ! -z "${_compmaxrank}" ]; then 
                #echo $_compmaxrank
                :
            elif [ ! -z "${compmaxrank[iexp]}" ]; then 
                _compmaxrank=${compmaxrank[iexp]}
            else
                ## calculate maxrank used for buffers during computation
                scaledb=$((_b/10))
                scaledmaxrank=$((_maxrank*4))
                val=$scaledmaxrank
                if [ $scaledb -lt $scaledmaxrank ]; then
                    val=$scaledb
                fi
                if [ $val -le $_wavek ]; then
                    val=$((_wavek+50))
                fi
                _compmaxrank=$val
                _compmaxrank=$((_b/2))
            fi
            if [ -z "$_m" ]; then 
                continue
            fi
            if [ -z "$_b" ]; then 
                continue
            fi
            #if [$(_m %_mb) != 0]
           #   _mb=$_b;
             # _ts=$(_m/_mb)
             # _n=$(ts*_maxrank)
            #  _nb=$_maxrank
            #else
            #  _mb=$_b;
            #  _n=$(((_m/_mb)*_maxrank))
            #  _nb=$_maxrank
            #fi
               _mb=$_b;
if [ $(( _m % _mb )) -eq 0 ]; then
               _n=$(((_m/_mb)*_maxrank))
             _nb=$_maxrank
else
    printf '%s\n' "Number $n is Odd"
              _ts=$(((_m/_mb)+1))
               _n=$((_ts*_maxrank))
              _nb=$_maxrank
fi
            echo $_n
            echo $_nb
            echo $_mb
        if [ $sqrt_numnodes -eq $sqrt_numnodesQ ]; then
            pdims=$sqrt_numnodes 
        else
            pdims="$sqrt_numnodes $sqrt_numnodesQ" 
        fi
        for pdim in $pdims; do
            if [ "$prog" == "hic" ]; then
                rankfile=$BINDIR/../../exp/ranks/$prog-$sched-$_m-$_mb-$_nb-$numnodes-$nummpi-$numthreads-$_acc-$_rad-$_order-$_ker-$SLURM_JOBID
                if [ "$op" == "potrf" ]; then
                    cmd="$BIN \
                        --m=$_m \
                        --n_range=$_n:$_n \
                        --k=$_m \
                        --mb=$_mb \
                        --nb=$_maxrank \
                        --nowarmup \
                        --threads=$numthreads \
                        --p=$pdim \
                        $tracearg \
                        --rk=0 \
                        --acc=$_acc \
                        $_appdata \
                        --starshwavek=$_wavek \
                        --starshdecay=$_decay \
                        --starshmaxrank=$_compmaxrank \
                        --rankfile=$rankfile \
                        --rbf_kernel=$_ker\
                        --rad=$_rad \
                        --numobj=$_numobj\
                        --order=$_order\
                        --mesh_file=$_mesh_file \
                        --csolve \
                        --solve \
                        #--check \
                        "
                elif [ "$op" == "posv" ]; then
                    cmd="$BIN \
                        $numthreads \
                        0 \
                        $op \
                        $_m \
                        $_m \
                        $_nrhs \
                        $_m \
                        $_mb \
                        1e-$_acc \
                        0 \
                        $_maxrank \
                        $_compmaxrank \
                        $pdim \
                        $((numnodes/pdim)) \
                        "
                elif [ "$op" == "genc" ]; then
                    cmd="$BIN \
                        --m=$_m \
                        --n_range=$_n:$_n \
                        --k=$_m \
                        --mb=$_mb \
                        --nb=$_maxrank \
                        --nowarmup \
                        --threads=$numthreads \
                        --p=$pdim \
                        $tracearg \
                        --rk=0 \
                        --acc=$_acc \
                        $_appdata \
                        --starshwavek=$_wavek \
                        --starshdecay=$_decay \
                        --starshmaxrank=$_compmaxrank \
                        --rankfile=$rankfile \
                        #--check \
                       "
                fi
            elif [ "$prog" == "cham" ]; then
                export STARPU_SCHED=$sched
                export STARPU_SCHED=eager ## prio does not work with chameleon potrf
                if [ "$trace" != "-" ]; then
                    cmd="/lustre/project/k1205/akbudak/hicma/chameleon/build-cdt-fxt-s13/timing/time_dpotrf_tile --nowarmup --P=$pdim --m=$_m --n_range=$_m:$_m:$_m --nb=$_b --threads=$numthreads $tracearg"
                else
                    cmd="$CHAMDIR/timing/time_dpotrf_tile --P=$pdim --m=$_m --n_range=$_m:$_m:$_m $cham_block_size_arg --threads=$numthreads --nowarmup "
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
                #export STARPU_SILENT=1

                msg="M:$_m N:$_n MB:$_mb NB:$_nb MAXRANK:$_maxrank DATE:`date` SCHED:$STARPU_SCHED CMD:$cmd $minmaxsubinfo CASE:$sizefile" 
                echo "!BEGIN:" $msg 
                if [ "$trace" != "-" ]; then
                    traceprefix=`pwd`/exp/trace/$prog-$sched-$_m-$_mb-$_nb-$numnodes-$nummpi-$numthreads/$SLURM_JOBID
                    mkdir -p $traceprefix
                    export STARPU_FXT_PREFIX=$traceprefix
                fi
                if [ "$dry" == "dry" ]; then
                    echo $cmd
                    #echo $sruncmd $cmd
                else
                    echo "!BEGIN:" $msg 1>&2 
                    tstart=$SECONDS
                    $sruncmd $cmd
                    tend=$SECONDS
                    time_sec=$((tend-tstart))
                    time_min=$((time_sec/60))
                    time_hour=$((time_min/60))
                    echo
                    echo "!END:" $msg SECOND:$time_sec MINUTE:$time_min HOUR:$time_hour
                    echo "!END:" $msg 1>&2 
                    
                    if [ "$trace" != "-" ]; then
                        combinedtrace=${traceprefix}trace
                        /project/k1205/akbudak/codes/starpu13-install-cdt-fxt/bin/starpu_fxt_tool -i ${traceprefix}prof_file_akbudak_* -o ${combinedtrace}
                        mv activity.data  dag.dot  data.rec  distrib.data  tasks.rec  trace.html paje.trace    $traceprefix/ #
                        echo $SLURM_JOBID > $traceprefix/0_jobid 
                        #echo "Dot is starting:"
                        #dot -Tpng $traceprefix/dag.dot -o $traceprefix/$prog.png
                        #echo "Dot ended"
                    fi
                fi
                date
                echo
            done
         done
        done ## for each accuracy
    done ## for each case
done ## for each scheduler

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

