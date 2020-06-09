for m in 21960 36162 43920 67110 73728 86028 96960 109761 119388 131781 147456 ; do
    echo $m;
    $HOME/hicma-dev/chameleon/build/timing/time_zgemm_tile --m=$m --k=$m   --n_range=$m:$m:$m
done
