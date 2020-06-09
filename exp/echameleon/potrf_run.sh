ROOT=`pwd`
irange="3 3";nt=2;_bs="289 578 1156"; 
for _i in `seq $irange`;do
    for _b in $_bs;do
        _is=$((_i*_i))
        m=$((_is*_b));
        n=$m
        echo $m $_b
        mpirun -np 4 $ROOT/chameleon/build/timing/time_dpotrf_tile 

        exit
            --m=$m \
            --n_range=$n:$n \
            --k=$m \
            --nb=$_b \
            --threads=$nt \
            
    done;
done
exit
            --nowarmup \
