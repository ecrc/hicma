nbs=(1120 1260 1440 1680 2520 1120 1260 1440 1680 2240 2520 2880 1080 1120 1260 1440 1680 1890 2160 2520 1040 1120 1560 1680 1820 2080 2730 1040 1080 1170 1440 1560 2080 2160 2340 1040 1120 1560 1680 1820 2080 2240 2730 1080 1170 1260 1560 1820 1890 2340 2520 2730 1040 1080 1170 1440 1560 2080 2160 2340 2880 1040 1120 1170 1260 1440 1560 1680 1820 2080 2340 2520 2730 1040 1080 1170 1260 1560 1680 1820 1890 2160 2340 2520 2730 )
nts=(9 8 7 6 4 18 16 14 12 9 8 7 28 27 24 21 18 16 14 12 42 39 28 26 24 21 16 54 52 48 39 36 27 26 24 84 78 56 52 48 42 39 32 91 84 78 63 54 52 42 39 36 108 104 96 78 72 54 52 48 39 126 117 112 104 91 84 78 72 63 56 52 48 189 182 168 156 126 117 108 104 91 84 78 72 )
lennbs=${#nbs[@]}
ncases=$((lennbs))
echo "Number of nbs:$lennbs Number of cases:$ncases"
for i in "${!nbs[@]}"; do
    __nb=${nbs[i]}
    __nt=${nts[i]}
    __m=$((__nb*__nt))
    __halfnb=$((__nb/2))
    __maxrank=$((__nb/2))
    #__maxrank=1000
    if [[ $__halfnb -lt $__maxrank ]]; then
        __maxrank=$__halfnb;
    fi
    __compmaxrank=$((__nb/3*2))
    _ci=$((i+1)); nrows[$_ci]=$__m;     nb[$_ci]=$__nb;     acc[$_ci]="1e-8"; maxrank[$_ci]=$__maxrank;   compmaxrank[$_ci]=$__compmaxrank; appdata[$_ci]="--st-3D-sqexp"; rbf_kernel[$_ci]="NA"; denst[$_ci]="NA"; rad[$_ci]="NA";mesh_file[$_ci]="NA";numobj[$_ci]="NA";order[$_ci]="NA"
done

nprocs="1"
allcaseids[1]="`seq 1 $ncases`"  
prog="hic"
#prog="mkl"
step=1
timelimit="02:00:00"
note="Hicma $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "

