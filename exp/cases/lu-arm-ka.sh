ic=1;nrows[$ic]=21960; nb[$ic]=1464; acc[$ic]=3; maxrank[$ic]=600; compmaxrank[ $ic]=1200; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=3
ic=3;nrows[$ic]=43920; nb[$ic]=1830; acc[$ic]=3; maxrank[$ic]=760; compmaxrank[ $ic]=1250; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=6
ic=4;nrows[$ic]=73728; nb[$ic]=2304; acc[$ic]=3; maxrank[$ic]=700; compmaxrank[ $ic]=1400; appdata[ $ic]="--ss"; ntrian[$ic]=12288; nipp[$ic]=6
ic=5;nrows[$ic]=86028; nb[$ic]=1284; acc[$ic]=3; maxrank[$ic]=600; compmaxrank[ $ic]=1200; appdata[ $ic]="--ss"; ntrian[$ic]=14338; nipp[$ic]=6
ic=6;nrows[$ic]=97080; nb[$ic]=4854; acc[$ic]=3; maxrank[$ic]=3000; compmaxrank[ $ic]=6000; appdata[ $ic]="--ss"; ntrian[$ic]=16180; nipp[$ic]=6
ic=7;nrows[$ic]=147456; nb[$ic]=3072; acc[$ic]=3; maxrank[$ic]=1000; compmaxrank[ $ic]=2000; appdata[ $ic]="--ss"; ntrian[$ic]=12288; nipp[$ic]=12
nprocs="1"

allcaseids[1]="`seq 1 6`"
#allcaseids[1]="1"

prog="hic"
op="getrf"


step=1
_appdata="--ss"; timelimit="10:00:00"
note="HiCMA-LU $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "


