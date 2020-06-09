PRE=""
factor=1
debug=0
mpi=0
if [ $# -ne 7 ]; then
    echo "Usage: factor[1,2,4...] debug[d,-] mpi[m,-] nthreads [cham,hic] [trace,-] [potrf,gemm]"
    exit -1
fi
factor=$1
debug=$2
mpi=$3
nthreads=$4
prog=$5
trace=$6
op=$7
nmpi=8
procP=
procQ=
if [ "$mpi" != "-" ]; then
    nmpi=$mpi
    if [ $nmpi -gt 1 ]; then
        nP=2
        nQ=$((nmpi/2))
    else
        nP=1
        nQ=$nmpi
    fi
    procP="--p=$nP"
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
    if [ $debug == "v" ]; then
        CMD="valgrind --gen-suppressions=all --suppressions=/home/kadir/hicma/exp/valgrind.supp "
        #CMD="valgrind "
    fi
fi
echo "|$CMD|"
tracearg=
if [ "$trace" != "-" ]; then
    tracearg="--trace"
    export STARPU_FXT_TRACE=1
    traceprefix=`pwd`/exp/trace/$prog
    export STARPU_FXT_PREFIX=$traceprefix

fi
app="--ss"
app="--edsin"
app="--geostat"
if [ "$op" == "gemm" ]; then
    app="--rnd"
fi
#app="--rnd" # for fixed rank without any rank descriptors

#irange="2 2";nta=$nthreads;_b=1156; acc=8; check=""
irange="2 2";  nta=$nthreads;_b=400; acc=5; check="--check"
#irange="20 100 20";  nta=$nthreads;_b=400; acc=5;  check="--check"
#irange="4 4";  nta=$nthreads;_b=400; acc=8; check=""
irange="4 4";  nta=$nthreads;_b=400; acc=8; check="--check"
#irange="3 3";  nta=$nthreads;_b=400; acc=8; check="--check"
#irange="6 6";  nta=$nthreads;_b=400; acc=6; wavek=100;check=""
#irange="9 9";  nta=$nthreads;_b=121; acc=4; check=""
#irange="3 3";  nta=$nthreads;_b=256; acc=6; check="--check"
#irange="2 2";nta=$nthreads;_b=400; acc=8; check=""
#irange="4 4";nta=$nthreads;_b=1156; acc=8; check=""
#irange="3 3"; nta="2"; _b=324; acc=8; check="--check"  #ALEKS'S SUGGESTION
#irange="3 3"; nta="2"; _b=256; acc=6; check="--check"  #FAST
#irange="24 24";  nta=$nthreads;_b=100; acc=8; check=""
customsize=0
for nt in $nta;do 
    n=$((m/factor))
    maxrank=$((nb/factor))
    mb=$_b;
    for _i in `seq $irange`;do
        _is=$((_i*_i))
        m=$((_is*_b));
        n=$((_is*_b/factor));
        maxrank=$((_b/factor))
        nrhs=$((_b/4))

        #maxrank is used for input matrix
        #maxrank=91 
        #maxrank=60 
        n=$((m/mb*maxrank))
        nb=$maxrank


        #m=1937; mb=401; n=1000; nb=200;
        #m=2005; mb=401; n=1000; nb=200;
        #compmaxrank is used for buffers used during computations
        compmaxrank=$((mb/2)) 
        if [ $customsize -eq 1 ]; then
            #acc=4;m=$((38+36)); mb=38; n=$((31+31)); nb=31; compmaxrank=$((mb)) 
            #acc=4;m=$((38+38)); mb=38; n=$((31+31)); nb=31; compmaxrank=$((mb)) 
            #acc=2;m=$((17+16)); mb=17; n=$((31+31)); nb=31; compmaxrank=$((mb)) 
            #acc=2;m=$((17+17)); mb=17; n=$((31+31)); nb=31; compmaxrank=$((mb)) 
            acc=4;nt=1;m=$((36*nt)); mb=36; n=$((31*nt)); nb=31; compmaxrank=$((mb)) 
            acc=4;nt=10;m=$((33*nt)); mb=33; n=$((27*nt)); nb=27; compmaxrank=$((mb)) 
            acc=8;    #potrf
            nt=4; mb=401; m=$((mb*nt)); nb=309; n=$((nb*nt)); compmaxrank=$((mb)); nrhs=10 
            nt=4; nb=400; n=1999; nrhs=1; acc=10; #trsmd
            nt=4; nb=400; n=1871; nrhs=1; acc=10; #trsmd
            #nt=4; nb=400; n=1600; nrhs=1; acc=10; #trsmd
            #nt=4; nb=400; n=800; nrhs=1; acc=8; #trsmd
            #nt=4; nb=400; n=671; nrhs=1; acc=8; #trsmd
            #nt=4; nb=40; n=79; nrhs=1; acc=3; #trsmd
            #nt=4; nb=40; n=80; nrhs=1; acc=4; #trsmd
            #nt=4; nb=14; n=28; nrhs=1; acc=3 #trsmd
            #nt=4; nb=14; n=27; nrhs=1; acc=3 #trsmd
            #m=$((401+390)); mb=401; n=$((321+321)); nb=321; compmaxrank=$((mb)) 
            #m=$((309)); mb=309; n=$((321+321)); nb=321; compmaxrank=$((mb)) 
            #m=$((401)); mb=401; n=$((321)); nb=321; compmaxrank=$((mb)) 
            #m=$((401+401)); mb=401; n=$((321+321)); nb=321; compmaxrank=$((mb)) 
        fi
        compmaxrank=$((mb)) 
        #compmaxrank=160

        if [ "$op" == "gemm" -o "$op" == "potrf" ]; then
            run="$HOME/hicma/build/timing/time_z${op}_tile \
                $procP \
                --m=$m \
                --n_range=$n:$n \
                --k=$m \
                --mb=$mb \
                --nb=$nb \
                --nowarmup \
                --threads=$nt \
                --rk=0 \
                --acc=$acc \
                $check \
                $app \
                $tracearg \
                --starshwavek=$wavek \
                --starshdecay=0.41 \
                --starshmaxrank=$compmaxrank \
                --rankfile=RANKS \
            "
                #--printindexend \
                #--printindexall \
                #--printmat \
            runcham="$HOME/hicma/chameleon/buildfxt/timing/time_d${op}_tile \
                $procP \
                --m=$m \
                --n_range=$m:$m \
                --nb=$nb \
                $tracearg \
                --nowarmup \
                --threads=$nt \
                "
        elif [ "$op" == "posv" ]; then
            run="$HOME/hicma/build/testing/testing_z${op} \
                1 0 $op $n $n $nrhs $n $nb 1e-$acc 0 $((nb/2)) $((nb/2)) 1 1 1 \
                "
            runcham="$HOME/hicma/chameleon/build/testing/dtesting \
                1 0 POSV 200 200 40 200 \
                "
        elif [ "$op" == "trsmd" ]; then
            run="$HOME/hicma/build/testing/testing_z${op} \
                1 0 $op $n $n $nrhs $n $nb 1e-$acc 0  \
                "
        fi
        combinedtrace=${traceprefix}trace
        if [ "$trace" != "-" ]; then
            rm ${traceprefix}prof_file_kadir_* 
            rm -f ${traceprefix}/*
        fi
        export STARPU_SCHED=ws
        #export STARPU_SCHED=prio
        export MORSE_TESTING_VERBOSE=1
        if [ "$prog" == "cham" ]; then
            $CMD $runcham
        elif [ "$prog" == "hic" ]; then
            echo $CMD $run
            $CMD $run
        fi
        if [ "$trace" != "-" ]; then
            starpu_fxt_tool -i ${traceprefix}prof_file_kadir_* -o ${combinedtrace}
            mkdir -p $traceprefix
            mv activity.data  dag.dot  data.rec  distrib.data  tasks.rec  trace.html $traceprefix/ #paje.trace    
        dot -Tpng $traceprefix/dag.dot -o $traceprefix/$prog.png

        fi
    done
done
exit 0



            --rankfile=RANKS \
    --alwaysfixedrank \
    
    --geostat \
    --ss \


    --check \
    --printmat \
    --printindexall \

    #--starshmaxrank=$nb

    #--rndusr \
    #--rnd \
    #--ss \

#not used in code
#--maxrank=$((nb/2)) \
