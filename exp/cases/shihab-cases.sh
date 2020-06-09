ic=1;nrows[$ic]=21960;  nb[$ic]=1830;   acc[$ic]=3;   maxrank[$ic]=600;  compmaxrank[ $ic]=1200; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=3
ic=2;nrows[$ic]=43920;  nb[$ic]=1830;   acc[$ic]=3;   maxrank[$ic]=800;  compmaxrank[ $ic]=1600; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=6
ic=3;nrows[$ic]=73728;  nb[$ic]=2304;   acc[$ic]=3;   maxrank[$ic]=600;  compmaxrank[ $ic]=1200; appdata[ $ic]="--ss"; ntrian[$ic]=12288; nipp[$ic]=6
ic=4;nrows[$ic]=86028;  nb[$ic]=1284;   acc[$ic]=3;   maxrank[$ic]=600;  compmaxrank[ $ic]=1200; appdata[ $ic]="--ss"; ntrian[$ic]=14338; nipp[$ic]=6
ic=5;nrows[$ic]=96960;  nb[$ic]=2424;   acc[$ic]=3;   maxrank[$ic]=1200;  compmaxrank[ $ic]=2400; appdata[ $ic]="--ss"; ntrian[$ic]=16160; nipp[$ic]=6
ic=6;nrows[$ic]=147456;  nb[$ic]=2304;   acc[$ic]=3;   maxrank[$ic]=800;  compmaxrank[ $ic]=1600; appdata[ $ic]="--ss"; ntrian[$ic]=12288; nipp[$ic]=12
nprocs="1"

allcaseids[1]="1 2 3 4 5 6" #19 20 21 36

prog="hic"
op="getrf"


step=1
_appdata="--ss"; timelimit="10:00:00"
note="Hicma beta=0.01 $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "
