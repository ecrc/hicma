#compmaxrank is IPARAM_HICMA_STARSH_MAXRANK, is used
#in Starsh matrix generation as a limit
#to allocated temporary buffer for hcore_gemm
#in checking the width of concatenated matrices in QR of HCORE_GEMM

#maxrank is used to determine the number of columns corresponding to rank of a tile.
#it is used in this formula: number_of_tiles * maxrank.
#this formula gives the total number of columns of the matrix storing low rank tiles
#this value is passed as --nb to program.

#nb is the number of ROWS in a tile.

_ci=0; nrows[$_ci]=10370;   nb[$_ci]=1037;    acc[$_ci]="1e-5";   maxrank[$_ci]=200;  compmaxrank[$_ci]=400; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.000460;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus10370.txt";  rbf_kernel[$_ci]="0"; numobj[$_ci]=1; denst[$_ci]=-1;    
_ci=1; nrows[$_ci]=10370;   nb[$_ci]=1037;    acc[$_ci]="1e-5";   maxrank[$_ci]=400;  compmaxrank[$_ci]=800; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.000460;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus10370.txt";  rbf_kernel[$_ci]="1"; numobj[$_ci]=1; denst[$_ci]=-1;
_ci=2; nrows[$_ci]=10370;   nb[$_ci]=1037;    acc[$_ci]="1e-5";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1200; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.000460;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus10370.txt";  rbf_kernel[$_ci]="2"; numobj[$_ci]=1; denst[$_ci]=-1;
_ci=3; nrows[$_ci]=10370;   nb[$_ci]=1037;    acc[$_ci]="1e-5";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1200; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.000460;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus10370.txt";  rbf_kernel[$_ci]="3"; numobj[$_ci]=1; denst[$_ci]=-1;
_ci=4; nrows[$_ci]=10370;   nb[$_ci]=1037;    acc[$_ci]="1e-5";   maxrank[$_ci]=500;  compmaxrank[$_ci]=1000; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.6;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus10370.txt";  rbf_kernel[$_ci]="7"; numobj[$_ci]=1; denst[$_ci]=-1;
_ci=5; nrows[$_ci]=10370;   nb[$_ci]=1037;    acc[$_ci]="1e-5";   maxrank[$_ci]=300;  compmaxrank[$_ci]=600; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.6;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus10370.txt";  rbf_kernel[$_ci]="9"; numobj[$_ci]=1; denst[$_ci]=-1;

_ci=6; nrows[$_ci]=44932;   nb[$_ci]=1000;    acc[$_ci]="1e-5";   maxrank[$_ci]=200;  compmaxrank[$_ci]=400; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.000370;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus44932.txt";  rbf_kernel[$_ci]="0"; numobj[$_ci]=1; denst[$_ci]=-1; 
_ci=7; nrows[$_ci]=44932;   nb[$_ci]=1000;    acc[$_ci]="1e-6 ";   maxrank[$_ci]=500;  compmaxrank[$_ci]=1000; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.000370;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus44932.txt";  rbf_kernel[$_ci]="1"; numobj[$_ci]=1; denst[$_ci]=-1; 
_ci=8; nrows[$_ci]=44932;   nb[$_ci]=1000;    acc[$_ci]="1e-7";   maxrank[$_ci]=800;  compmaxrank[$_ci]=1600; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.000370;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus44932.txt";  rbf_kernel[$_ci]="2"; numobj[$_ci]=1; denst[$_ci]=-1; 
_ci=9; nrows[$_ci]=44932;   nb[$_ci]=1000;    acc[$_ci]="1e-5";   maxrank[$_ci]=800;  compmaxrank[$_ci]=1600; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.000370;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus44932.txt";  rbf_kernel[$_ci]="3"; numobj[$_ci]=1; denst[$_ci]=-1;
_ci=10; nrows[$_ci]=44932;   nb[$_ci]=1000;    acc[$_ci]="1e-5";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1200; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.6;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus44932.txt";  rbf_kernel[$_ci]="7"; numobj[$_ci]=1; denst[$_ci]=-1;
_ci=11; nrows[$_ci]=44932;   nb[$_ci]=1000;    acc[$_ci]="1e-5";   maxrank[$_ci]=300;  compmaxrank[$_ci]=600; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.6;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus44932.txt";  rbf_kernel[$_ci]="9"; numobj[$_ci]=1; denst[$_ci]=-1;

_ci=12; nrows[$_ci]=117715;   nb[$_ci]=1811;    acc[$_ci]="1e-5";   maxrank[$_ci]=300;  compmaxrank[$_ci]=600; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.000200;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus117715.txt"; rbf_kernel[$_ci]="0";  numobj[$_ci]=1; denst[$_ci]=-1; 
_ci=13; nrows[$_ci]=117715;   nb[$_ci]=1811;    acc[$_ci]="1e-6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1200; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.000200;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus117715.txt"; rbf_kernel[$_ci]="1";  numobj[$_ci]=1;denst[$_ci]=-1; 
_ci=14; nrows[$_ci]=117715;   nb[$_ci]=1811;    acc[$_ci]="1e-7";   maxrank[$_ci]=900;  compmaxrank[$_ci]=1800; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.000200;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus117715.txt"; rbf_kernel[$_ci]="2";  numobj[$_ci]=1; denst[$_ci]=-1; 
_ci=15; nrows[$_ci]=117715;   nb[$_ci]=1811;    acc[$_ci]="1e-5";   maxrank[$_ci]=900;  compmaxrank[$_ci]=1800; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.000200;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus117715.txt"; rbf_kernel[$_ci]="3";  numobj[$_ci]=1; denst[$_ci]=-1;
_ci=16; nrows[$_ci]=117715;   nb[$_ci]=1811;    acc[$_ci]="1e-5";   maxrank[$_ci]=800;  compmaxrank[$_ci]=1600; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.6;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus117715.txt"; rbf_kernel[$_ci]="7";  numobj[$_ci]=1; denst[$_ci]=-1;
_ci=17; nrows[$_ci]=117715;   nb[$_ci]=1811;    acc[$_ci]="1e-5";   maxrank[$_ci]=400;  compmaxrank[$_ci]=800; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.6;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus117715.txt"; rbf_kernel[$_ci]="9";  numobj[$_ci]=1; denst[$_ci]=-1;

_ci=18; nrows[$_ci]=142416;   nb[$_ci]=2400;    acc[$_ci]="1e-5";   maxrank[$_ci]=500;  compmaxrank[$_ci]=1000; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.000036;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus142418.txt"; rbf_kernel[$_ci]="0"; numobj[$_ci]=1; denst[$_ci]=-1; 
_ci=19; nrows[$_ci]=142416;   nb[$_ci]=2400;    acc[$_ci]="1e-6";   maxrank[$_ci]=800;  compmaxrank[$_ci]=1600; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.000036;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus142418.txt"; rbf_kernel[$_ci]="1"; numobj[$_ci]=1; denst[$_ci]=-1; 
_ci=20; nrows[$_ci]=142416;   nb[$_ci]=2400;    acc[$_ci]="1e-7";   maxrank[$_ci]=1000;  compmaxrank[$_ci]=2000; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.000036;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus142418.txt"; rbf_kernel[$_ci]="2"; numobj[$_ci]=1; denst[$_ci]=-1; 
_ci=21; nrows[$_ci]=142416;   nb[$_ci]=2400;    acc[$_ci]="1e-7";   maxrank[$_ci]=1000;  compmaxrank[$_ci]=2000; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.000036;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus142418.txt"; rbf_kernel[$_ci]="3"; numobj[$_ci]=1; denst[$_ci]=-1;
_ci=22; nrows[$_ci]=142416;   nb[$_ci]=2400;    acc[$_ci]="1e-7";   maxrank[$_ci]=800;  compmaxrank[$_ci]=1600; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.6;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus142418.txt"; rbf_kernel[$_ci]="7"; numobj[$_ci]=1; denst[$_ci]=-1;
_ci=23; nrows[$_ci]=142416;   nb[$_ci]=2400;    acc[$_ci]="1e-7";   maxrank[$_ci]=400;  compmaxrank[$_ci]=800; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=0.6;  order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus142418.txt"; rbf_kernel[$_ci]="9"; numobj[$_ci]=1; denst[$_ci]=-1;



norocs="1 2 4 8 16"
nprocs="1"

allcaseids[1]="`seq 0 23`"

prog="hic"

allcaseids[2]="`seq 1 4` `seq 9 12`"

step=1
timelimit="00:30:00"
#_compmaxrank=150 #for 54000 maxrank=100
note="Hicma beta=0.01 $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "

#prog="hic"
#prog="cham"

