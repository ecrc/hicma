ic=1; nrows[$ic]=21960;  nb[$ic]=1464; acc[$ic]=3; maxrank[$ic]=1464; compmaxrank[ $ic]=2928;  appdata[ $ic]="--ss"; ntrian[$ic]=7320;  nipp[$ic]=3
ic=2; nrows[$ic]=21960;  nb[$ic]=1830; acc[$ic]=3; maxrank[$ic]=1830; compmaxrank[ $ic]=3660;  appdata[ $ic]="--ss"; ntrian[$ic]=7320;  nipp[$ic]=3
ic=3; nrows[$ic]=36162;  nb[$ic]=882;  acc[$ic]=3; maxrank[$ic]=500;  compmaxrank[ $ic]=1000;  appdata[ $ic]="--ss"; ntrian[$ic]=6027;  nipp[$ic]=6
ic=4; nrows[$ic]=43920;  nb[$ic]=1830; acc[$ic]=3; maxrank[$ic]=800;  compmaxrank[ $ic]=1600;  appdata[ $ic]="--ss"; ntrian[$ic]=7320;  nipp[$ic]=6
ic=5; nrows[$ic]=66216;  nb[$ic]=2136; acc[$ic]=3; maxrank[$ic]=1000; compmaxrank[ $ic]=2000;  appdata[ $ic]="--ss"; ntrian[$ic]=5518;  nipp[$ic]=12
ic=6; nrows[$ic]=67110;  nb[$ic]=6711; acc[$ic]=3; maxrank[$ic]=5000; compmaxrank[ $ic]=10000; appdata[ $ic]="--ss"; ntrian[$ic]=22370; nipp[$ic]=3
ic=7; nrows[$ic]=67116;  nb[$ic]=1428; acc[$ic]=3; maxrank[$ic]=1000; compmaxrank[ $ic]=2000;  appdata[ $ic]="--ss"; ntrian[$ic]=5593;  nipp[$ic]=12
ic=8; nrows[$ic]=73728;  nb[$ic]=2304; acc[$ic]=3; maxrank[$ic]=1500; compmaxrank[ $ic]=3000;  appdata[ $ic]="--ss"; ntrian[$ic]=12288; nipp[$ic]=6
ic=9; nrows[$ic]=86028;  nb[$ic]=1284; acc[$ic]=3; maxrank[$ic]=900;  compmaxrank[ $ic]=1800;  appdata[ $ic]="--ss"; ntrian[$ic]=14338; nipp[$ic]=6
ic=10;nrows[$ic]=96960;  nb[$ic]=2424; acc[$ic]=3; maxrank[$ic]=1500; compmaxrank[ $ic]=3000;  appdata[ $ic]="--ss"; ntrian[$ic]=16160; nipp[$ic]=6
ic=11;nrows[$ic]=97080;  nb[$ic]=4854; acc[$ic]=3; maxrank[$ic]=3000; compmaxrank[ $ic]=6000;  appdata[ $ic]="--ss"; ntrian[$ic]=16180; nipp[$ic]=6
ic=12;nrows[$ic]=110316; nb[$ic]=3804; acc[$ic]=3; maxrank[$ic]=1500; compmaxrank[ $ic]=3000;  appdata[ $ic]="--ss"; ntrian[$ic]=9193;  nipp[$ic]=12
ic=13;nrows[$ic]=110424; nb[$ic]=2568; acc[$ic]=3; maxrank[$ic]=1500; compmaxrank[ $ic]=3000;  appdata[ $ic]="--ss"; ntrian[$ic]=9202;  nipp[$ic]=12
ic=14;nrows[$ic]=110424; nb[$ic]=1284; acc[$ic]=3; maxrank[$ic]=1000; compmaxrank[ $ic]=2000;  appdata[ $ic]="--ss"; ntrian[$ic]=9202;  nipp[$ic]=12
ic=15;nrows[$ic]=118272; nb[$ic]=2112; acc[$ic]=3; maxrank[$ic]=1000; compmaxrank[ $ic]=2000;  appdata[ $ic]="--ss"; ntrian[$ic]=9856;  nipp[$ic]=12
ic=16;nrows[$ic]=131781; nb[$ic]=1209; acc[$ic]=3; maxrank[$ic]=1200; compmaxrank[ $ic]=2400;  appdata[ $ic]="--ss"; ntrian[$ic]=43927; nipp[$ic]=3
ic=17;nrows[$ic]=131781; nb[$ic]=4251; acc[$ic]=3; maxrank[$ic]=1200; compmaxrank[ $ic]=2400;  appdata[ $ic]="--ss"; ntrian[$ic]=43927; nipp[$ic]=3 
ic=18;nrows[$ic]=147546; nb[$ic]=2304; acc[$ic]=3; maxrank[$ic]=600;  compmaxrank[ $ic]=1200;  appdata[ $ic]="--ss"; ntrian[$ic]=49152; nipp[$ic]=3 

nprocs="1"

##compaxrank was 1600 but it was not enough for acc=6 so I set it to 2000
#ic=7;nrows[$ic]=147456;  nb[$ic]=2304;   acc[$ic]="6";   maxrank[$ic]=800;  compmaxrank[ $ic]=2000; appdata[ $ic]="--ss"; ntrian[$ic]=12288; nipp[$ic]=12
allcaseids[1]="7" #different accuracy thresholds

allcaseids[1]="`seq 1 18`" 

prog="hic"
op="getrf"


step=1
_appdata="--ss"; timelimit="10:00:00"
note="Hicma beta=0.01 $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "


