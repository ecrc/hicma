#for _m_ in 12000 24000 36000 48000; do
for _m_ in 24000; do
    for _mb_ in 600 1200; do
        _maxrank_=1200
        for _acc_ in 3; do
            for reorderinnerproducts in 0 1; do
                cmd="./timing/time_zgemm_tile --m=$_m_ --n_range=$_m_:$_m_ --k=$_m_ --mb=$_mb_ --nb=$_mb_ --nowarmup --threads=55 --starshmaxrank=$_maxrank_ --st-3D-exp --acc=$_acc_ --reorderinnerproducts=$reorderinnerproducts"
                echo $cmd
                eval $cmd
            done
        done
    done
done
