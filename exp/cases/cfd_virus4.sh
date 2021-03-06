#compmaxrank is IPARAM_HICMA_STARSH_MAXRANK, is used
#in Starsh matrix generation as a limit
#to allocated temporary buffer for hcore_gemm
#in checking the width of concatenated matrices in QR of HCORE_GEMM

#maxrank is used to determine the number of columns corresponding to rank of a tile.
#it is used in this formula: number_of_tiles * maxrank.
#this formula gives the total number of columns of the matrix storing low rank tiles
#this value is passed as --nb to program.

#nb is the number of ROWS in a tile.



_ci=40; nrows[$_ci]=103700;   nb[$_ci]=2074;    acc[$_ci]="1e-5 1e-6 1e-7";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=10; order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/S1data/SortVirus103700.txt";  rbf_kernel[$_ci]="9"; 
_ci=41; nrows[$_ci]=207400   nb[$_ci]=3050;    acc[$_ci]="1e-5 1e-6 1e-7";   maxrank[$_ci]=100;  compmaxrank[$_ci]=200; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=20; order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/S1data/SortVirus207400.txt"; rbf_kernel[$_ci]="9"; 


nprocs="4"
allcaseids[4]="40 41" 


prog="hic"
step=1
timelimit="04:00:00"
note="Hicma beta=0.01 $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "


