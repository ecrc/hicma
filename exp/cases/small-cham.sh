nrows[1]=24;  nb[1]=12; acc[1]=1;   maxrank[1]=6;
nrows[2]=36;  nb[2]=12; acc[2]=1;   maxrank[2]=6;
nrows[3]=1024;  nb[3]=64; acc[3]=5;   maxrank[3]=32;
nrows[4]=12;  nb[4]=12; acc[4]=1;   maxrank[4]=6;

nrows[1]=24;   nb[1]=12;     acc[1]=1;   maxrank[1]=6;   compmaxrank[ 1]=6; appdata[ 1]="--ss" 
nrows[2]=72;   nb[2]=36;     acc[2]=4;   maxrank[2]=18;   compmaxrank[ 2]=18; appdata[ 2]="--ss" 
#nrows[3]=1024;   nb[3]=512;     acc[3]=5;   maxrank[3]=256;   compmaxrank[3]=256; appdata[ 3]="--ss"
nrows[4]=36;   nb[4]=12;     acc[4]=1;   maxrank[4]=6;   compmaxrank[ 4]=6; appdata[ 4]="--ss" 
ic=5;nrows[$ic]=60;   nb[$ic]=20;     acc[$ic]=1;   maxrank[$ic]=10;   compmaxrank[ $ic]=10; appdata[ $ic]="--ss" 
ic=6;nrows[$ic]=192;   nb[$ic]=64;     acc[$ic]=2;   maxrank[$ic]=32;   compmaxrank[ $ic]=40; appdata[ $ic]="--ss" 
ic=7;nrows[$ic]=2048;   nb[$ic]=256;     acc[$ic]=8;   maxrank[$ic]=128;   compmaxrank[ $ic]=200; appdata[ $ic]="--ss" 
nprocs="1"
allcaseids[1]="7"

prog="hic"
prog="cham"
op="getrf"


step=1
_appdata="--ss"; timelimit="10:00:00"
note="Hicma beta=0.01 $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "


