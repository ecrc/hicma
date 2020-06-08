#compmaxrank is IPARAM_HICMA_STARSH_MAXRANK, is used
#in Starsh matrix generation as a limit
#to allocated temporary buffer for hcore_gemm
#in checking the width of concatenated matrices in QR of HCORE_GEMM

#maxrank is used to determine the number of columns corresponding to rank of a tile.
#it is used in this formula: number_of_tiles * maxrank.
#this formula gives the total number of columns of the matrix storing low rank tiles
#this value is passed as --nb to program.

#nb is the number of ROWS in a tile.



_ci=40; nrows[$_ci]=235430;   nb[$_ci]=3622;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=2; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file";  rbf_kernel[$_ci]="9"; 
_ci=41; nrows[$_ci]=470860;   nb[$_ci]=400;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=4; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file";  rbf_kernel[$_ci]="9"; 
_ci=42; nrows[$_ci]=706290;   nb[$_ci]=5120;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=6; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file";  rbf_kernel[$_ci]="9"; 
_ci=43; nrows[$_ci]=941720;   nb[$_ci]=5120;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=8; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file";  rbf_kernel[$_ci]="9"; 
_ci=44; nrows[$_ci]=1177150;   nb[$_ci]=5120;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=10; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file";  rbf_kernel[$_ci]="9"; 
_ci=45; nrows[$_ci]=1412580;   nb[$_ci]=5120;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=12; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file";  rbf_kernel[$_ci]="9"; 
_ci=46; nrows[$_ci]=1648010;   nb[$_ci]=5120;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=14; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file";  rbf_kernel[$_ci]="9"; 
_ci=47; nrows[$_ci]=1883440;   nb[$_ci]=5120;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=16; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file";  rbf_kernel[$_ci]="9"; 
_ci=48; nrows[$_ci]=2118870;   nb[$_ci]=5120;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=18; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file";  rbf_kernel[$_ci]="9"; 
_ci=49; nrows[$_ci]=2354300;   nb[$_ci]=5120;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=20; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file";  rbf_kernel[$_ci]="9"; 
_ci=50; nrows[$_ci]=2825160;   nb[$_ci]=5120;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=24; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file";  rbf_kernel[$_ci]="9"; 
_ci=51; nrows[$_ci]=3296020;   nb[$_ci]=5120;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=28; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file";  rbf_kernel[$_ci]="9"; 
_ci=52; nrows[$_ci]=4002310;   nb[$_ci]=5120;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=34; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file"; rbf_kernel[$_ci]="9"; 
_ci=53; nrows[$_ci]=4237740;   nb[$_ci]=5120;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=36; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file"; rbf_kernel[$_ci]="9"; 
_ci=54; nrows[$_ci]=4237740;   nb[$_ci]=7244;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=36; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file"; rbf_kernel[$_ci]="9"; 
_ci=55; nrows[$_ci]=4473170;   nb[$_ci]=5120;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=38; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file"; rbf_kernel[$_ci]="9"; 
_ci=56; nrows[$_ci]=4473170;   nb[$_ci]=7244;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=38; order[$_ci]=2; mesh_file[$_ci]="path/to/mesh/file"; rbf_kernel[$_ci]="9"; 



norocs="1 2 4 8 16"
nprocs="32"
allcaseids[32]="`seq 40 51`" 

prog="hic"
step=1
timelimit="04:00:00"
note="Hicma beta=0.01 $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "


