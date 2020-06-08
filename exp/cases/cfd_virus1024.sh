#compmaxrank is IPARAM_HICMA_STARSH_MAXRANK, is used
#in Starsh matrix generation as a limit
#to allocated temporary buffer for hcore_gemm
#in checking the width of concatenated matrices in QR of HCORE_GEMM

#maxrank is used to determine the number of columns corresponding to rank of a tile.
#it is used in this formula: number_of_tiles * maxrank.
#this formula gives the total number of columns of the matrix storing low rank tiles
#this value is passed as --nb to program.

#nb is the number of ROWS in a tile.





_ci=48; nrows[$_ci]=4002310;   nb[$_ci]=5120;    acc[$_ci]="1e-6 1e-8";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=34; order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/4002310.txt"; rbf_kernel[$_ci]="9"; reg[$_ci]=1.1;

_ci=49; nrows[$_ci]=4237740;   nb[$_ci]=5120;    acc[$_ci]="1e-6 1e-8";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=36; order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/4002310.txt"; rbf_kernel[$_ci]="9"; reg[$_ci]=1.1;

_ci=50; nrows[$_ci]=4237740;   nb[$_ci]=7244;    acc[$_ci]="1e-6 1e-8";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=36; order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/4237740.txt"; rbf_kernel[$_ci]="9"; reg[$_ci]=1.1;

_ci=51; nrows[$_ci]=4473170;   nb[$_ci]=5120;    acc[$_ci]="1e-6 1e-8";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=38; order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/4473170.txt"; rbf_kernel[$_ci]="9"; reg[$_ci]=1.1;

_ci=52; nrows[$_ci]=4473170;   nb[$_ci]=7244;    acc[$_ci]="1e-6 1e-8";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=38; order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/4473170.txt"; rbf_kernel[$_ci]="9"; reg[$_ci]=1.1;


nprocs="1024"
allcaseids[1024]="48 49 50 51 52" 


prog="hic"
step=1
timelimit="04:00:00"
note="Hicma beta=0.01 $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "


