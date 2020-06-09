nrows[1]=24;  nb[1]=12; acc[1]=1;   maxrank[1]=6;
nrows[2]=36;  nb[2]=12; acc[2]=1;   maxrank[2]=6;
nrows[3]=1024;  nb[3]=64; acc[3]=5;   maxrank[3]=32;
nrows[4]=12;  nb[4]=12; acc[4]=1;   maxrank[4]=6;

nrows[1]=24;   nb[1]=12;     acc[1]=1;   maxrank[1]=6;   compmaxrank[ 1]=6; appdata[ 1]="--ss"
nrows[2]=72;   nb[2]=36;     acc[2]=4;   maxrank[2]=18;   compmaxrank[ 2]=18; appdata[ 2]="--ss"
#nrows[3]=1024;   nb[3]=512;     acc[3]=5;   maxrank[3]=256;   compmaxrank[3]=256; appdata[ 3]="--ss"
nrows[4]=36;   nb[4]=12;     acc[4]=1;   maxrank[4]=6;   compmaxrank[ 4]=6; appdata[ 4]="--ss"
ic=5;nrows[$ic]=6;      nb[$ic]=6;     acc[$ic]=2;   maxrank[$ic]=3;   compmaxrank[ $ic]=3; appdata[ $ic]="--ss"
ic=6;nrows[$ic]=12;     nb[$ic]=12;    acc[$ic]=2;   maxrank[$ic]=12;  compmaxrank[ $ic]=12; appdata[ $ic]="--ss"
ic=7;nrows[$ic]=2048;   nb[$ic]=1024;  acc[$ic]=8;   maxrank[$ic]=512; compmaxrank[ $ic]=512; appdata[ $ic]="--ss"
ic=8;nrows[$ic]=3;      nb[$ic]=3;     acc[$ic]=2;   maxrank[$ic]=3;   compmaxrank[ $ic]=3; appdata[ $ic]="--ss"
ic=9;nrows[$ic]=6;      nb[$ic]=3;     acc[$ic]=2;   maxrank[$ic]=3;   compmaxrank[ $ic]=3; appdata[ $ic]="--ss"
ic=10;nrows[$ic]=9;     nb[$ic]=3;     acc[$ic]=2;   maxrank[$ic]=3;   compmaxrank[ $ic]=3; appdata[ $ic]="--ss" ntrian[$ic]=3 nipp[$ic]=3
ic=11;nrows[$ic]=3;     nb[$ic]=3;     acc[$ic]=2;   maxrank[$ic]=3;   compmaxrank[ $ic]=3; appdata[ $ic]="--ss" ntrian[$ic]=3 nipp[$ic]=3
ic=12;nrows[$ic]=18;    nb[$ic]=6;     acc[$ic]=2;   maxrank[$ic]=6;   compmaxrank[ $ic]=10; appdata[ $ic]="--ss" ntrian[$ic]=6 nipp[$ic]=3

#one tile
ic=13;nrows[$ic]=3;   nb[$ic]=3;     acc[$ic]=2;   maxrank[$ic]=3;   compmaxrank[ $ic]=3;  appdata[ $ic]="--ss"; ntrian[$ic]=1;  nipp[$ic]=3
ic=14;nrows[$ic]=90;  nb[$ic]=90;    acc[$ic]=2;   maxrank[$ic]=90;  compmaxrank[ $ic]=90; appdata[ $ic]="--ss"; ntrian[$ic]=30; nipp[$ic]=3
ic=17;nrows[$ic]=27;  nb[$ic]=27;    acc[$ic]=4;   maxrank[$ic]=9;  compmaxrank[ $ic]=15; appdata[ $ic]="--ss"; ntrian[$ic]=9;  nipp[$ic]=3

#two tiles
ic=15;nrows[$ic]=6;   nb[$ic]=3;     acc[$ic]=2;   maxrank[$ic]=3;   compmaxrank[ $ic]=3;  appdata[ $ic]="--ss"; ntrian[$ic]=2;  nipp[$ic]=3
ic=16;nrows[$ic]=90;  nb[$ic]=45;    acc[$ic]=6;   maxrank[$ic]=45;  compmaxrank[ $ic]=50; appdata[ $ic]="--ss"; ntrian[$ic]=30; nipp[$ic]=3
ic=18;nrows[$ic]=18;  nb[$ic]=9;     acc[$ic]=5;   maxrank[$ic]=45;  compmaxrank[ $ic]=50; appdata[ $ic]="--ss"; ntrian[$ic]=6;  nipp[$ic]=3

#three tiles
ic=19;nrows[$ic]=18;  nb[$ic]=6;    acc[$ic]=5;   maxrank[$ic]=6;   compmaxrank[ $ic]=10;  appdata[ $ic]="--ss"; ntrian[$ic]=6; nipp[$ic]=3
ic=20;nrows[$ic]=27;  nb[$ic]=9;    acc[$ic]=4;   maxrank[$ic]=9;   compmaxrank[ $ic]=20;  appdata[ $ic]="--ss"; ntrian[$ic]=9; nipp[$ic]=3
ic=39;nrows[$ic]=36;  nb[$ic]=12;    acc[$ic]=4;   maxrank[$ic]=12;   compmaxrank[ $ic]=30;  appdata[ $ic]="--ss"; ntrian[$ic]=12; nipp[$ic]=3
ic=40;nrows[$ic]=45;  nb[$ic]=15;    acc[$ic]=3;   maxrank[$ic]=15;   compmaxrank[ $ic]=50;  appdata[ $ic]="--ss"; ntrian[$ic]=15; nipp[$ic]=3
ic=41;nrows[$ic]=63;  nb[$ic]=21;    acc[$ic]=4;   maxrank[$ic]=21;   compmaxrank[ $ic]=50;  appdata[ $ic]="--ss"; ntrian[$ic]=21; nipp[$ic]=3
ic=42;nrows[$ic]=72;  nb[$ic]=27;    acc[$ic]=4;   maxrank[$ic]=27;   compmaxrank[ $ic]=50;  appdata[ $ic]="--ss"; ntrian[$ic]=27;nipp[$ic]=3
ic=21;nrows[$ic]=90;  nb[$ic]=30;   acc[$ic]=5;   maxrank[$ic]=30;  compmaxrank[ $ic]=80; appdata[ $ic]="--ss"; ntrian[$ic]=30; nipp[$ic]=3
#connot be run on Matlab since these are supported for now 2,3,9,12,120,396,622 and 1126
ic=22;nrows[$ic]=120;  nb[$ic]=30;   acc[$ic]=3;   maxrank[$ic]=30;  compmaxrank[ $ic]=90; appdata[ $ic]="--ss"; ntrian[$ic]=40; nipp[$ic]=3;

#four tiles
ic=43;nrows[$ic]=36;  nb[$ic]=9;    acc[$ic]=3;   maxrank[$ic]=9;   compmaxrank[ $ic]=20;  appdata[ $ic]="--ss"; ntrian[$ic]=12; nipp[$ic]=3
ic=44;nrows[$ic]=72;  nb[$ic]=18;    acc[$ic]=3;   maxrank[$ic]=18;   compmaxrank[ $ic]=50;  appdata[ $ic]="--ss"; ntrian[$ic]=24; nipp[$ic]=3
ic=45;nrows[$ic]=720;  nb[$ic]=180;    acc[$ic]=4;   maxrank[$ic]=180;   compmaxrank[ $ic]=500;  appdata[ $ic]="--ss"; ntrian[$ic]=120; nipp[$ic]=6

ic=24;nrows[$ic]=1440;  nb[$ic]=24;   acc[$ic]=3;   maxrank[$ic]=24;  compmaxrank[ $ic]=50; appdata[ $ic]="--ss"; ntrian[$ic]=120; nipp[$ic]=12
ic=60;nrows[$ic]=1440;  nb[$ic]=240;   acc[$ic]=3;   maxrank[$ic]=240;  compmaxrank[ $ic]=500; appdata[ $ic]="--ss"; ntrian[$ic]=120; nipp[$ic]=12
ic=25;nrows[$ic]=9;     nb[$ic]=3;     acc[$ic]=2;   maxrank[$ic]=3;   compmaxrank[ $ic]=3; appdata[ $ic]="--ss" ntrian[$ic]=3 nipp[$ic]=3
ic=26;nrows[$ic]=360;     nb[$ic]=24;     acc[$ic]=3;   maxrank[$ic]=24;   compmaxrank[ $ic]=48; appdata[ $ic]="--ss" ntrian[$ic]=120 nipp[$ic]=3
ic=53;nrows[$ic]=360;     nb[$ic]=60;     acc[$ic]=4;   maxrank[$ic]=60;   compmaxrank[ $ic]=150; appdata[ $ic]="--ss" ntrian[$ic]=120 nipp[$ic]=3
ic=27;nrows[$ic]=720;     nb[$ic]=120;     acc[$ic]=4;   maxrank[$ic]=120;   compmaxrank[ $ic]=240; appdata[ $ic]="--ss" ntrian[$ic]=120 nipp[$ic]=6
ic=28;nrows[$ic]=1440;     nb[$ic]=240;     acc[$ic]=3;   maxrank[$ic]=120;   compmaxrank[ $ic]=192; appdata[ $ic]="--ss" ntrian[$ic]=120 nipp[$ic]=12
ic=29;nrows[$ic]=1188;     nb[$ic]=24;     acc[$ic]=2;   maxrank[$ic]=6;   compmaxrank[ $ic]=6; appdata[ $ic]="--ss" ntrian[$ic]=396 nipp[$ic]=3
ic=30;nrows[$ic]=2376;     nb[$ic]=24;     acc[$ic]=2;   maxrank[$ic]=6;   compmaxrank[ $ic]=6; appdata[ $ic]="--ss" ntrian[$ic]=396 nipp[$ic]=6
ic=31;nrows[$ic]=4752;     nb[$ic]=24;     acc[$ic]=2;   maxrank[$ic]=6;   compmaxrank[ $ic]=6; appdata[ $ic]="--ss" ntrian[$ic]=396 nipp[$ic]=12
ic=32;nrows[$ic]=1866;     nb[$ic]=24;     acc[$ic]=2;   maxrank[$ic]=6;   compmaxrank[ $ic]=6; appdata[ $ic]="--ss" ntrian[$ic]=622 nipp[$ic]=3
ic=33;nrows[$ic]=3732;     nb[$ic]=24;     acc[$ic]=2;   maxrank[$ic]=6;   compmaxrank[ $ic]=6; appdata[ $ic]="--ss" ntrian[$ic]=622 nipp[$ic]=6
ic=34;nrows[$ic]=7464;     nb[$ic]=24;     acc[$ic]=2;   maxrank[$ic]=6;   compmaxrank[ $ic]=6; appdata[ $ic]="--ss" ntrian[$ic]=622 nipp[$ic]=12
ic=35;nrows[$ic]=2400;     nb[$ic]=36;     acc[$ic]=2;   maxrank[$ic]=36;   compmaxrank[ $ic]=36; appdata[ $ic]="--ss" ntrian[$ic]=400 nipp[$ic]=6

ic=36;nrows[$ic]=1440;  nb[$ic]=240;   acc[$ic]=3;   maxrank[$ic]=240;  compmaxrank[ $ic]=500; appdata[ $ic]="--ss"; ntrian[$ic]=120; nipp[$ic]=12
ic=37;nrows[$ic]=1440;  nb[$ic]=240;   acc[$ic]=3;   maxrank[$ic]=240;  compmaxrank[ $ic]=240; appdata[ $ic]="--ss"; ntrian[$ic]=120; nipp[$ic]=12
ic=38;nrows[$ic]=7680;  nb[$ic]=384;   acc[$ic]="2 4 6 8";   maxrank[$ic]=384;  compmaxrank[ $ic]=1000; appdata[ $ic]="--ss"; ntrian[$ic]=2560; nipp[$ic]=3
ic=46;nrows[$ic]=21960;  nb[$ic]=549;   acc[$ic]=2;   maxrank[$ic]=549;  compmaxrank[ $ic]=1098; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=3
ic=47;nrows[$ic]=21960;  nb[$ic]=549;   acc[$ic]=2;   maxrank[$ic]=360;  compmaxrank[ $ic]=720; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=3
ic=48;nrows[$ic]=21960;  nb[$ic]=732;   acc[$ic]=3;   maxrank[$ic]=732;  compmaxrank[ $ic]=1464; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=3
ic=49;nrows[$ic]=21960;  nb[$ic]=915;   acc[$ic]=3;   maxrank[$ic]=570;  compmaxrank[ $ic]=1240; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=3
ic=50;nrows[$ic]=21960;  nb[$ic]=1098;   acc[$ic]=3;   maxrank[$ic]=500;  compmaxrank[ $ic]=1000; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=3
ic=51;nrows[$ic]=21960;  nb[$ic]=1464;   acc[$ic]=2;   maxrank[$ic]=500;  compmaxrank[ $ic]=1000; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=3
ic=52;nrows[$ic]=43920;  nb[$ic]=1830;   acc[$ic]=3;   maxrank[$ic]=1281;  compmaxrank[ $ic]=2562; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=6
ic=53;nrows[$ic]=43920;  nb[$ic]=1830;   acc[$ic]=3;   maxrank[$ic]=800;  compmaxrank[ $ic]=1600; appdata[ $ic]="--ss"; ntrian[$ic]=7320; nipp[$ic]=6
ic=54;nrows[$ic]=1866;  nb[$ic]=48;   acc[$ic]="2 4 6 8";   maxrank[$ic]=48;  compmaxrank[ $ic]=100; appdata[ $ic]="--ss"; ntrian[$ic]=622; nipp[$ic]=3
ic=55;nrows[$ic]=1866;  nb[$ic]=24;   acc[$ic]="2";   maxrank[$ic]=48;  compmaxrank[ $ic]=100; appdata[ $ic]="--ss"; ntrian[$ic]=622; nipp[$ic]=3
ic=57;nrows[$ic]=31992;  nb[$ic]=1032;   acc[$ic]="2";   maxrank[$ic]=1032;  compmaxrank[ $ic]=2000; appdata[ $ic]="--ss"; ntrian[$ic]=5332; nipp[$ic]=6
ic=58;nrows[$ic]=7680; nb[$ic]=768; acc[$ic]="3"; maxrank[$ic]=768; compmaxrank[ $ic]=1536; appdata[ $ic]="--ss"; ntrian[$ic]=2560; nipp[$ic]=3
nprocs="1"

allcaseids[1]="22" #19 36
allcaseids[1]="27" #19 20 21 36
#nprocs="4"
allcaseids[1]="58" 
#allcaseids[2]="60" 
#allcaseids[4]="52" #46 
#allcaseids[8]="52"  #53

BINNAME=time_zgetrf_tile

prog="hic"
op="getrf"


step=1
_appdata="--ss"; 
timelimit="01:00:00";que="workq"
#timelimit="00:10:00";que="debug"
#timelimit="00:20:00";que="workq"
note="Hicma beta=0.01 $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "
