ic=1;nrows[$ic]=21960;  nb[$ic]=1464;   acc[$ic]=3;   maxrank[$ic]=1024;  compmaxrank[ $ic]=2048; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=3
ic=2;nrows[$ic]=21960;  nb[$ic]=2196;   acc[$ic]=3;   maxrank[$ic]=1500;  compmaxrank[ $ic]=3000; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=3
nprocs="1"
ic=3;nrows[$ic]=43920;  nb[$ic]=1830;   acc[$ic]=3;   maxrank[$ic]=1281;  compmaxrank[ $ic]=2562; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=6

ic=4;nrows[$ic]=73728;  nb[$ic]=2304;   acc[$ic]=3;   maxrank[$ic]=900;  compmaxrank[ $ic]=1800; appdata[ $ic]="--ss"; ntrian[$ic]=12288; nipp[$ic]=6

ic=5;nrows[$ic]=103422;  nb[$ic]=4701;   acc[$ic]=2;   maxrank[$ic]=900;  compmaxrank[ $ic]=1800; appdata[ $ic]="--ss"; ntrian[$ic]=34474; nipp[$ic]=3

allcaseids[1]="5" #19 20 21 36

prog="hic"
op="getrf"


step=1
_appdata="--ss"; timelimit="10:00:00"
note="Hicma beta=0.01 $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "


