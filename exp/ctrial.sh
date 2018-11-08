#!/bin/bash
#SBATCH --job-name=cham
#SBATCH --output=/project/k1205/akbudak/hicma/exp/cout/%j.o
#SBATCH --error=/project/k1205/akbudak/hicma/exp/cerr/%j.e
#SBATCH --partition=workq
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time 12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=akbudak


#3teraflops
export STARPU_SCHED=eager
#export STARPU_SCHED=prio
#export STARPU_MIN_PRIO=-1
#export STARPU_MAX_PRIO=100

#These do not affect perf for dpotrf M=36000 mb=1000
export STARPU_CALIBRATE=0
export STARPU_LIMIT_CPU_MEM=60000
export STARPU_LIMIT_MAX_SUBMITTED_TASKS=20000
export STARPU_LIMIT_MIN_SUBMITTED_TASKS=18000

export LD_LIBRARY_PATH=/opt/intel/composer_xe_2015.2.164/mkl/lib/intel64/:/project/k1205/akbudak/codes/chameleon/install/lib/:/project/k1205/akbudak/codes/hwloc-1.11.7/install/lib/:/project/k1205/akbudak/codes/starpu-1.2.1/install/lib/:$LD_LIBRARY_PATH
echo STARPU_SCHED                     $STARPU_SCHED
echo STARPU_CALIBRATE                 $STARPU_CALIBRATE
echo STARPU_LIMIT_CPU_MEM             $STARPU_LIMIT_CPU_MEM
echo STARPU_LIMIT_MAX_SUBMITTED_TASKS $STARPU_LIMIT_MAX_SUBMITTED_TASKS
echo STARPU_LIMIT_MIN_SUBMITTED_TASKS $STARPU_LIMIT_MIN_SUBMITTED_TASKS
kernel=time_dgemm_tile
kernel=time_dpotrf_tile
#try with these
nb=460
_nb=320
_nb=320
#_nb=1156
_b=1156
for _i in `seq 12 20`; do
for _nb in 320 460 1156;do
    _is=$((_i*_i))
    _m=$((_is*_b))
    echo "# $i $_m $nb"
    echo -n "#" date
    STARPU_SILENT=1  numactl --interleave=all  srun --hint=nomultithread  /project/k1205/akbudak/hicma/chameleon/build/timing/$kernel --n_range=$_m:$_m:$_m --m=$_m --k=$_m --nb=$_nb --p=2 --threads=31
    echo -n "#" date
done
done
