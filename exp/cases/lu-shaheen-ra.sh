ic=1;nrows[$ic]=21960;  nb[$ic]=1464;   acc[$ic]=3;   maxrank[$ic]=1464;  compmaxrank[ $ic]=2928; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=3
ic=2;nrows[$ic]=21960;  nb[$ic]=1830;   acc[$ic]=3;   maxrank[$ic]=1830;  compmaxrank[ $ic]=3660; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=3
ic=3;nrows[$ic]=43920;  nb[$ic]=1830;   acc[$ic]=3;   maxrank[$ic]=1281;  compmaxrank[ $ic]=2562; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=6
ic=4;nrows[$ic]=67110;  nb[$ic]=6711;   acc[$ic]=3;   maxrank[$ic]=5000;  compmaxrank[ $ic]=10000; appdata[ $ic]="--ss"; ntrian[$ic]=22370; nipp[$ic]=3
ic=5;nrows[$ic]=73728;  nb[$ic]=2304;   acc[$ic]=3;   maxrank[$ic]=1728;  compmaxrank[ $ic]=3456; appdata[ $ic]="--ss"; ntrian[$ic]=12288; nipp[$ic]=6
ic=6;nrows[$ic]=86028;  nb[$ic]=1284;   acc[$ic]=3;   maxrank[$ic]=1284;  compmaxrank[ $ic]=2568; appdata[ $ic]="--ss"; ntrian[$ic]=14338; nipp[$ic]=6
nprocs="1"

allcaseids[1]="4" #19 20 21 36

prog="hic"
op="getrf"


step=1
_appdata="--ss"; timelimit="10:00:00"
note="Hicma beta=0.01 $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "
