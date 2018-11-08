PRE=""
factor=1
debug=0
mpi=0
if [ $# -ne 4 ]; then
    echo "Usage: factor[1,2,4...] debug[d,-] mpi[m,-] nthreads"
    exit -1
fi
factor=$1
debug=$2
mpi=$3
nthreads=$4
nmpi=8
if [ "$mpi" != "-" ]; then
    nmpi=$mpi
    if [ $debug == "d" ]; then
        CMD="mpirun -n $nmpi xterm -hold -e gdb -ex run --args "
    else 
        #CMD="mpirun -n $nmpi -tag-output  "
        #CMD="mpirun -n $nmpi -quiet  "
        CMD="mpirun -n $nmpi  "
    fi
else
    if [ $debug == "d" ]; then
        CMD="gdb -ex run --args"
    fi
fi
echo $CMD

export STARPU_SILENT=1
nb=2;m=4
nb=3;m=6
#nb=3;m=9
nb=4;m=4
nb=4;m=8
nb=9;m=36
acc=3;nb=16;m=64
acc=4;nb=128;m=1024 
#acc=7;nb=128;m=1024 
#acc=6;nb=128;m=4096 
#nb=7;m=21    #works for rndusr potrf
nta="1 4 8" 
nta="27 28 29 30 31 32" #shihab has 36 physical cores 
#irange="2 2";nta=$nthreads;_b=1156; acc=8; check=""
irange="3 3";  nta=$nthreads;_b=400; acc=8; check="--check"
#irange="4 4";nta=$nthreads;_b=1156; acc=8; check=""
#irange="3 3"; nta="2"; _b=324; acc=8; check="--check"  #ALEKS'S SUGGESTION
#irange="3 3"; nta="2"; _b=256; acc=6; check="--check"  #FAST
#nta="1"
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
        run="./build/main \
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
            --ss \
            --starshdecay=2 \
            --starshmaxrank=$maxrank"
        runcham="/home/akbudak/hicma/chameleon/build/timing/time_dpotrf_tile \
            --m=$m \
            --n_range=$n:$n \
            --nowarmup \
            --threads=$nt \
            "
        $CMD $run
    done
done
exit
    --check \
    --printmat \
    --printindex \

    #--starshmaxrank=$nb

    #--rndusr \
    #--rnd \
    #--ss \

#not used in code
#--maxrank=$((nb/2)) \
