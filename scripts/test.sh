#!/bin/bash -le
module purge
if [ "$HOSTNAME" == "thana" ]; then
	. ./scripts/power8.modules
else
	. ./scripts/intel.modules
fi

module list

set -x

# TEST
PRE=""
factor=1
debug=0
mpi=0


if [ $# -ne 5 ]; then
    echo "Usage: factor[1,2,4...] debug[d,-] mpi[m,-] nthreads problem"
    ii=0
    for i in $*; do 
       echo "$ii: |$i|"
       ii=$((ii+1))
    done
    exit -1
fi

factor=$1
debug=$2
mpi=$3
nthreads=$4
problem=$5

nmpi=8
if [ "$mpi" != "-" ]; then
    nmpi=$mpi
    if [ $debug == "d" ]; then
        CMD="mpirun -n $nmpi xterm -hold -e gdb -ex run --args "
    else 
        CMD="mpirun -n $nmpi  "
    fi
else
    if [ $debug == "d" ]; then
        CMD="gdb -ex run --args"
    fi
fi
echo $CMD

export STARPU_SILENT=1
irange="3 3";  nta=$nthreads;_b=400; acc=8; check="--check"
for nt in $nta;do 
    n=$((m/factor))
    maxrank=$((nb/factor))
    #echo BASH-MAXRANK: $maxrank
    #echo BASH-DEBUG: $debug
    #_b=1600
    #_b=324
    nb=$_b;
    for _i in `seq $irange`;do
        _is=$((_i*_i))
        m=$((_is*_b));
        n=$((_is*_b/factor));
        maxrank=$((_b/factor))
        run="./build/timing/time_zpotrf_tile \
            --m=$m \
            --n_range=$n:$n \
            --k=$m \
            --mb=$nb \
            --nb=$maxrank \
            --nowarmup \
            --threads=$nt \
            --rk=0 \
            --acc=$acc \
            $check \
            $problem \
            --starshwavek=40 \
            --starshdecay=2 \
            --starshmaxrank=$maxrank"
        echo "Executing: $CMD $run"
        $CMD $run
    done
done
