#compmaxrank is IPARAM_HICMA_STARSH_MAXRANK, is used
#in Starsh matrix generation as a limit
#to allocated temporary buffer for hcore_gemm
#in checking the width of concatenated matrices in QR of HCORE_GEMM

#maxrank is used to determine the number of columns corresponding to rank of a tile.
#it is used in this formula: number_of_tiles * maxrank.
#this formula gives the total number of columns of the matrix storing low rank tiles
#this value is passed as --nb to program.

#nb is the number of ROWS in a tile.


_ci=40; nrows[$_ci]=20740;   nb[$_ci]=1037;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=8; numsubobj[$_ci]=2; order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/S2data_20k/"; rbf_kernel[$_ci]="9"; denst[$_ci]=-1;

_ci=41; nrows[$_ci]=207400   nb[$_ci]=3050;    acc[$_ci]="1e-6";   maxrank[$_ci]=50;  compmaxrank[$_ci]=100; appdata[$_ci]="--m-3D-rbf"; rad[$_ci]=-1; numobj[$_ci]=60; numsubobj[$_ci]=240; order[$_ci]=2; mesh_file[$_ci]="stars-h/SARS-CoV-2-meshes/S2data_200k/"; rbf_kernel[$_ci]="9"; denst[$_ci]=-1;


nprocs="1"
allcaseids[1]="40" 


prog="hic"
step=1
timelimit="04:00:00"
note="Hicma beta=0.01 $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "


