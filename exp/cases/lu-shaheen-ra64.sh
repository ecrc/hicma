ic=1;nrows[$ic]=21960;  nb[$ic]=1464;   acc[$ic]=3;   maxrank[$ic]=1464;  compmaxrank[ $ic]=2928; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=3
ic=2;nrows[$ic]=21960;  nb[$ic]=1830;   acc[$ic]=3;   maxrank[$ic]=1830;  compmaxrank[ $ic]=3660; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=3
ic=3;nrows[$ic]=43920;  nb[$ic]=1830;   acc[$ic]=3;   maxrank[$ic]=1281;  compmaxrank[ $ic]=2562; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=6
ic=4;nrows[$ic]=67110;  nb[$ic]=6711;   acc[$ic]=3;   maxrank[$ic]=5000;  compmaxrank[ $ic]=10000; appdata[ $ic]="--ss"; ntrian[$ic]=22370; nipp[$ic]=3
ic=5;nrows[$ic]=73728;  nb[$ic]=2304;   acc[$ic]=3;   maxrank[$ic]=1200;  compmaxrank[ $ic]=2400; appdata[ $ic]="--ss"; ntrian[$ic]=12288; nipp[$ic]=6
ic=6;nrows[$ic]=73728;  nb[$ic]=2304;   acc[$ic]=3;   maxrank[$ic]=1400;  compmaxrank[ $ic]=2700; appdata[ $ic]="--ss"; ntrian[$ic]=12288; nipp[$ic]=6

ic=7;nrows[$ic]=86028;  nb[$ic]=1284;   acc[$ic]=3;   maxrank[$ic]=1200;  compmaxrank[ $ic]=2400; appdata[ $ic]="--ss"; ntrian[$ic]=14338; nipp[$ic]=6

ic=8;nrows[$ic]=86028;  nb[$ic]=1284;   acc[$ic]=3;   maxrank[$ic]=1400;  compmaxrank[ $ic]=2700; appdata[ $ic]="--ss"; ntrian[$ic]=14338; nipp[$ic]=6

ic=9;nrows[$ic]=147456;  nb[$ic]=2304;   acc[$ic]=2;   maxrank[$ic]=1200;  compmaxrank[ $ic]=2400; appdata[ $ic]="--ss"; ntrian[$ic]=12288; nipp[$ic]=12
ic=10;nrows[$ic]=147456;  nb[$ic]=2304;   acc[$ic]=2;   maxrank[$ic]=1400;  compmaxrank[ $ic]=2700; appdata[ $ic]="--ss"; ntrian[$ic]=12288; nipp[$ic]=12

ic=11;nrows[$ic]=187452;  nb[$ic]=1476;   acc[$ic]=2;   maxrank[$ic]=1200;  compmaxrank[ $ic]=2400; appdata[ $ic]="--ss"; ntrian[$ic]=15621; nipp[$ic]=12

ic=12;nrows[$ic]=187452;  nb[$ic]=1476;   acc[$ic]=2;   maxrank[$ic]=1200;  compmaxrank[ $ic]=2400; appdata[ $ic]="--ss"; ntrian[$ic]=15621; nipp[$ic]=12

ic=13;nrows[$ic]=190230;  nb[$ic]=2238;   acc[$ic]=2;   maxrank[$ic]=1800;  compmaxrank[ $ic]=3600; appdata[ $ic]="--ss"; ntrian[$ic]=63410; nipp[$ic]=3
ic=14;nrows[$ic]=190230;  nb[$ic]=2238;   acc[$ic]=2;   maxrank[$ic]=2000;  compmaxrank[ $ic]=4000; appdata[ $ic]="--ss"; ntrian[$ic]=63410; nipp[$ic]=3

ic=15;nrows[$ic]=168744;  nb[$ic]=1896;   acc[$ic]=2;   maxrank[$ic]=600;  compmaxrank[ $ic]=1200; appdata[ $ic]="--ss"; ntrian[$ic]=12288; nipp[$ic]=12

nprocs="64"

#allcaseids[64]="1 3 5 6 7 8 9 10 11 12 13 14" #19 20 21 36
allcaseids[64]="13 14" #19 20 21 36

prog="hic"
op="getrf"


step=1
_appdata="--ss"; timelimit="10:00:00"
note="Hicma beta=0.01 $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "
