sizes="big"
sizes="small"
if [ "$sizes" == "small"  ];then
#compmaxrank is IPARAM_HICMA_STARSH_MAXRANK, is used
#in Starsh matrix generation as a limit
#to allocated temporary buffer for hcore_gemm
#in checking the width of concatenated matrices in QR of HCORE_GEMM

#maxrank is used to determine the number of columns corresponding to rank of a tile.
#it is used in this formula: number_of_tiles * maxrank.
#this formula gives the total number of columns of the matrix storing low rank tiles
#this value is passed as --nb to program.

#nb is the number of ROWS in a tile.

_ci=101; nrows[$_ci]=54000;    nb[$_ci]=1350;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss"; 
_ci=101; nrows[$_ci]=5400;     nb[$_ci]=135;     acc[$_ci]="6";    maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss"; 
_ci=101; nrows[$_ci]=1350;     nb[$_ci]=135;     acc[$_ci]="3e-8"; maxrank[$_ci]=145;   compmaxrank[$_ci]=250; appdata[$_ci]="--ss"; 
_ci=102; nrows[$_ci]=54000;    nb[$_ci]=1500;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss"; 
_ci=103; nrows[$_ci]=54000;    nb[$_ci]=1800;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss"; 
_ci=104; nrows[$_ci]=54000;    nb[$_ci]=2250;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=105; nrows[$_ci]=81000;    nb[$_ci]=1350;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=106; nrows[$_ci]=81000;    nb[$_ci]=1500;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=107; nrows[$_ci]=81000;    nb[$_ci]=1800;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=108; nrows[$_ci]=81000;    nb[$_ci]=2250;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=109; nrows[$_ci]=108000;   nb[$_ci]=1350;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=110; nrows[$_ci]=108000;   nb[$_ci]=1500;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=111; nrows[$_ci]=108000;   nb[$_ci]=1800;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=112; nrows[$_ci]=108000;   nb[$_ci]=2250;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=113; nrows[$_ci]=135000;   nb[$_ci]=1350;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=114; nrows[$_ci]=135000;   nb[$_ci]=1500;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=115; nrows[$_ci]=135000;   nb[$_ci]=1800;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=116; nrows[$_ci]=135000;   nb[$_ci]=2250;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=117; nrows[$_ci]=162000;   nb[$_ci]=1350;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=118; nrows[$_ci]=162000;   nb[$_ci]=1500;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=119; nrows[$_ci]=162000;   nb[$_ci]=1800;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=120; nrows[$_ci]=162000;   nb[$_ci]=2250;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=121; nrows[$_ci]=189000;   nb[$_ci]=1350;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=122; nrows[$_ci]=189000;   nb[$_ci]=1500;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=123; nrows[$_ci]=189000;   nb[$_ci]=1800;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";
_ci=124; nrows[$_ci]=189000;   nb[$_ci]=2250;    acc[$_ci]="4 6";  maxrank[$_ci]=100;   compmaxrank[$_ci]=150; appdata[$_ci]="--ss";

_ci=201; nrows[$_ci]=54000;    nb[$_ci]=1350;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp"; 
_ci=202; nrows[$_ci]=54000;    nb[$_ci]=1500;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp"; 
_ci=203; nrows[$_ci]=54000;    nb[$_ci]=1800;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp"; 
_ci=204; nrows[$_ci]=54000;    nb[$_ci]=2250;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=205; nrows[$_ci]=81000;    nb[$_ci]=1350;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=206; nrows[$_ci]=81000;    nb[$_ci]=1500;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=207; nrows[$_ci]=81000;    nb[$_ci]=1800;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=208; nrows[$_ci]=81000;    nb[$_ci]=2250;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=209; nrows[$_ci]=108000;   nb[$_ci]=1350;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=210; nrows[$_ci]=108000;   nb[$_ci]=1500;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=211; nrows[$_ci]=108000;   nb[$_ci]=1800;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=212; nrows[$_ci]=108000;   nb[$_ci]=2250;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=213; nrows[$_ci]=135000;   nb[$_ci]=1350;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=214; nrows[$_ci]=135000;   nb[$_ci]=1500;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=215; nrows[$_ci]=135000;   nb[$_ci]=1800;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=216; nrows[$_ci]=135000;   nb[$_ci]=2250;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=217; nrows[$_ci]=162000;   nb[$_ci]=1350;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=218; nrows[$_ci]=162000;   nb[$_ci]=1500;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=219; nrows[$_ci]=162000;   nb[$_ci]=1800;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=220; nrows[$_ci]=162000;   nb[$_ci]=2250;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=221; nrows[$_ci]=189000;   nb[$_ci]=1350;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=222; nrows[$_ci]=189000;   nb[$_ci]=1500;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=223; nrows[$_ci]=189000;   nb[$_ci]=1800;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
_ci=224; nrows[$_ci]=189000;   nb[$_ci]=2250;    acc[$_ci]="4 6";   maxrank[$_ci]=600;  compmaxrank[$_ci]=1100; appdata[$_ci]="--st-3D-sqexp";
                                                                 
nrows[25]=216000;   nb[25]=1350;    acc[25]=8;  maxrank[25]=100; 
nrows[26]=216000;   nb[26]=1500;    acc[26]=8;  maxrank[26]=100;
nrows[27]=216000;   nb[27]=1800;    acc[27]=8;  maxrank[27]=100;
nrows[28]=216000;   nb[28]=2250;    acc[28]=8;  maxrank[28]=100;
nrows[29]=243000;   nb[29]=1350;    acc[29]=8;  maxrank[29]=100;
nrows[30]=243000;   nb[30]=1500;    acc[30]=8;  maxrank[30]=100;
nrows[31]=243000;   nb[31]=1800;    acc[31]=8;  maxrank[31]=100;
nrows[32]=243000;   nb[32]=2250;    acc[32]=8;  maxrank[32]=100;
nrows[33]=270000;   nb[33]=1350;    acc[33]=8;  maxrank[33]=100;
nrows[34]=270000;   nb[34]=1500;    acc[34]=8;  maxrank[34]=100;
nrows[35]=270000;   nb[35]=1800;    acc[35]=8;  maxrank[35]=100;
nrows[36]=270000;   nb[36]=2250;    acc[36]=8;  maxrank[36]=100;
nrows[37]=297000;   nb[37]=1350;    acc[37]=8;  maxrank[37]=100;
nrows[38]=297000;   nb[38]=1500;    acc[38]=8;  maxrank[38]=100;
nrows[39]=297000;   nb[39]=1800;    acc[39]=8;  maxrank[39]=100;
nrows[40]=297000;   nb[40]=2250;    acc[40]=8;  maxrank[40]=100;
nrows[41]=324000;   nb[41]=1350;    acc[41]=8;  maxrank[41]=100;
nrows[42]=324000;   nb[42]=1500;    acc[42]=8;  maxrank[42]=100;
nrows[43]=324000;   nb[43]=1800;    acc[43]=8;  maxrank[43]=100;
nrows[44]=324000;   nb[44]=2250;    acc[44]=8;  maxrank[44]=100;
nrows[45]=351000;   nb[45]=1350;    acc[45]=8;  maxrank[45]=100;
nrows[46]=351000;   nb[46]=1500;    acc[46]=8;  maxrank[46]=100;
nrows[47]=351000;   nb[47]=1800;    acc[47]=8;  maxrank[47]=100;
nrows[48]=351000;   nb[48]=2250;    acc[48]=8;  maxrank[48]=100;
nrows[49]=378000;   nb[49]=1350;    acc[49]=8;  maxrank[49]=100;
nrows[50]=378000;   nb[50]=1500;    acc[50]=8;  maxrank[50]=100;
nrows[51]=378000;   nb[51]=1800;    acc[51]=8;  maxrank[51]=100;
nrows[52]=378000;   nb[52]=2250;    acc[52]=8;  maxrank[52]=100;
nrows[53]=405000;   nb[53]=1350;    acc[53]=8;  maxrank[53]=100;
nrows[54]=405000;   nb[54]=1500;    acc[54]=8;  maxrank[54]=100;
nrows[55]=405000;   nb[55]=1800;    acc[55]=8;  maxrank[55]=100;
nrows[56]=405000;   nb[56]=2250;    acc[56]=8;  maxrank[56]=100;
nrows[57]=432000;   nb[57]=1350;    acc[57]=8;  maxrank[57]=100;
nrows[58]=432000;   nb[58]=1500;    acc[58]=8;  maxrank[58]=100;
nrows[59]=432000;   nb[59]=1800;    acc[59]=8;  maxrank[59]=100;
nrows[60]=432000;   nb[60]=2250;    acc[60]=8;  maxrank[60]=100;
nrows[61]=459000;   nb[61]=1350;    acc[61]=8;  maxrank[61]=100;
nrows[62]=459000;   nb[62]=1500;    acc[62]=8;  maxrank[62]=100;
nrows[63]=459000;   nb[63]=1800;    acc[63]=8;  maxrank[63]=100;
nrows[64]=459000;   nb[64]=2250;    acc[64]=8;  maxrank[64]=100;
nrows[65]=486000;   nb[65]=1350;    acc[65]=8;  maxrank[65]=100;
nrows[66]=486000;   nb[66]=1500;    acc[66]=8;  maxrank[66]=100;
nrows[67]=486000;   nb[67]=1800;    acc[67]=8;  maxrank[67]=100;
nrows[68]=486000;   nb[68]=2250;    acc[68]=8;  maxrank[68]=100;
nrows[69]=513000;   nb[69]=1350;    acc[69]=8;  maxrank[69]=100;
nrows[70]=513000;   nb[70]=1500;    acc[70]=8;  maxrank[70]=100;
nrows[71]=513000;   nb[71]=1800;    acc[71]=8;  maxrank[71]=100;
nrows[72]=513000;   nb[72]=2250;    acc[72]=8;  maxrank[72]=100;
nrows[73]=540000;   nb[73]=1350;    acc[73]=8;  maxrank[73]=100;
nrows[74]=540000;   nb[74]=1500;    acc[74]=8;  maxrank[74]=100;
nrows[75]=540000;   nb[75]=1800;    acc[75]=8;  maxrank[75]=100;
nrows[76]=540000;   nb[76]=2250;    acc[76]=8;  maxrank[76]=100;
nrows[77]=567000;   nb[77]=1350;    acc[77]=8;  maxrank[77]=100;
nrows[78]=567000;   nb[78]=1500;    acc[78]=8;  maxrank[78]=100;
nrows[79]=567000;   nb[79]=1800;    acc[79]=8;  maxrank[79]=100;
nrows[80]=567000;   nb[80]=2250;    acc[80]=8;  maxrank[80]=100;
nrows[81]=594000;   nb[81]=1350;    acc[81]=8;  maxrank[81]=100;
nrows[82]=594000;   nb[82]=1500;    acc[82]=8;  maxrank[82]=100;
nrows[83]=594000;   nb[83]=1800;    acc[83]=8;  maxrank[83]=100;
nrows[84]=594000;   nb[84]=2250;    acc[84]=8;  maxrank[84]=100;
_ci=555; nrows[$_ci]=10370;     nb[$_ci]=1037;     acc[$_ci]="1e-8"; maxrank[$_ci]=500;   compmaxrank[$_ci]=1000; appdata[$_ci]="--ss"; rbf_kernel[$_ci]="NA"; denst[$_ci]="NA"; rad[$_ci]="NA";mesh_file[$_ci]="NA";numobj[$_ci]="NA";order[$_ci]="NA"
_ci=556; nrows[$_ci]=10370;     nb[$_ci]=1037;     acc[$_ci]="1e-8"; maxrank[$_ci]=500;   compmaxrank[$_ci]=1000; appdata[$_ci]="--st-2D-exp"; rbf_kernel[$_ci]="NA"; denst[$_ci]="NA"; rad[$_ci]="NA";mesh_file[$_ci]="NA";numobj[$_ci]="NA";order[$_ci]="NA"
nbs=(792 858 936 1144 1287 784 1274 1456 1568 2548 896 952 1088 1792 1904 2176 3808 900 1125 1350 1500 1620 2025 2250 2700 3375 4050 4500 1296 1404 1872 1944 2106 2808 3159 3888 4212 5616 1377 1683 1782 1836 2244 2754 3366 3564 5049 5508 1683 2079 2142 2618 3213 3366 3927 4158 5049 1836 2244 2376 2448 2992 3366 3672 4488 4752 5049 2673 2754 3366 4131 5049 5346 2244 2295 2805 2970 3060 3366 3740 4590 5049 5610 5940 2646 3087 3969 4116 5292 2754 3366 3564 3672 4488 5049 5508 3366 3861 3978
4862 5049 5967 3213 3366 3927 4158 4284 5049 5236 3366 4455 4590 5049 5610 3672 4488 4752 4896 5049 5984 3900 4290 4400 5200 5720 4200 4400 4620 5280 5600 5775 4368 4680 4914 5040 5460 5616 )
nts=(13 12 11 9 8 26 16 14 13 8 34 32 28 17 16 14 8 45 36 30 27 25 20 18 15 12 10 9 39 36 27 26 24 18 16 13 12 9 44 36 34 33 27 22 18 17 12 11 42 34 33 27 22 21 18 17 14 44 36 34 33 27 24 22 18 17 16 34 33 27 22 18 17 45 44 36 34 33 30 27 22 20 18 17 42 36 28 27 21 44 36 34 33 27 24 22 39 34 33 27 26 22 44 42 36 34 33 28 27 45 34 33 30 27 44 36 34 33 32 27 44 40 39 33 30 44 42 40 35 33 32 45 42 40 39 36 35 )
nbs=(360 390 420 520 540 560 630 720 780 840 910 1040 1080 1170 1260 1560 1680 1820 1890 2160 2340 2520 2730 3120 3510 3640 3780 4680 5040 5460 )
nts=(546 504 468 378 364 351 312 273 252 234 216 189 182 168 156 126 117 108 104 91 84 78 72 63 56 54 52 42 39 36 )
lennbs=${#nbs[@]}
ncases=$((lennbs*2))
echo "Number of nbs:$lennbs Number of cases:$ncases"
for i in "${!nbs[@]}"; do
    __nb=${nbs[i]}
    __nt=${nts[i]}
    __m=$((__nb*__nt))
    __halfnb=$((__nb/2))
    __maxrank=250
    if [[ $__halfnb -lt $__maxrank ]]; then
        __maxrank=$__halfnb;
    fi
    _ci=$((i+1)); nrows[$_ci]=$__m;     nb[$_ci]=$__nb;     acc[$_ci]="1e-8"; maxrank[$_ci]=$__maxrank;   compmaxrank[$_ci]=$__maxrank; appdata[$_ci]="--ss"; rbf_kernel[$_ci]="NA"; denst[$_ci]="NA"; rad[$_ci]="NA";mesh_file[$_ci]="NA";numobj[$_ci]="NA";order[$_ci]="NA"
    #_ci=$((i+1)); nrows[$_ci]=$__m;     nb[$_ci]=$__nb;     acc[$_ci]="1e-8"; maxrank[$_ci]=500;   compmaxrank[$_ci]=500; appdata[$_ci]="--ss"; rbf_kernel[$_ci]="NA"; denst[$_ci]="NA"; rad[$_ci]="NA";mesh_file[$_ci]="NA";numobj[$_ci]="NA";order[$_ci]="NA"
    __maxrank=800
    if [[ $__halfnb -lt $__maxrank ]]; then
        __maxrank=$__halfnb;
    fi
    _ci=$((i+1+lennbs)); nrows[$_ci]=$__m;     nb[$_ci]=$__nb;     acc[$_ci]="1e-8"; maxrank[$_ci]=$__maxrank;   compmaxrank[$_ci]=$__maxrank; appdata[$_ci]="--st-2D-exp"; rbf_kernel[$_ci]="NA"; denst[$_ci]="NA"; rad[$_ci]="NA";mesh_file[$_ci]="NA";numobj[$_ci]="NA";order[$_ci]="NA"
done

allcaseids[16]="`seq 1 84`"
allcaseids[16]="`seq 1 24`"
allcaseids[16]="`seq 1 4`"
#allcaseids[16]="1"
nprocs="1 2 4 8 16"
nprocs="1"
allcaseids[1]="`seq 101 124` `seq 201 220`"
allcaseids[1]="555 556"  
allcaseids[1]="`seq $lennbs $ncases`"  
allcaseids[1]="`seq 1 $ncases`"  
#allcaseids[1]="`seq 201 220`" # @HATEM st-3d-sqexp two accuracies, 6 matrix sizes
prog="hic"
#prog="mkl"

allcaseids[2]="`seq 1 4` `seq 9 12`"
allcaseids[4]="`seq 1 4` `seq 9 12`"
allcaseids[8]="`seq 1 4` `seq 9 12`"
allcaseids[16]="`seq 14 16` `seq 37 40`"
else
nrows[1]=1080000;   nb[1]=2700;     acc[1]=8;   maxrank[1]=50;   compmaxrank[ 1]=150; appdata[ 1]="--ss" 
nrows[2]=1080000;   nb[2]=3000;     acc[2]=8;   maxrank[2]=50;   compmaxrank[ 2]=150; appdata[ 2]="--ss"
nrows[3]=1080000;   nb[3]=3375;     acc[3]=8;   maxrank[3]=50;   compmaxrank[ 3]=150; appdata[ 3]="--ss"
nrows[4]=1080000;   nb[4]=4500;     acc[4]=8;   maxrank[4]=50;   compmaxrank[ 4]=150; appdata[ 4]="--ss"
nrows[5]=2295000;   nb[5]=2700;     acc[5]=8;   maxrank[5]=50;   compmaxrank[ 5]=150; appdata[ 5]="--ss"
nrows[6]=2295000;   nb[6]=3000;     acc[6]=8;   maxrank[6]=50;   compmaxrank[ 6]=150; appdata[ 6]="--ss"
nrows[7]=2295000;   nb[7]=3375;     acc[7]=8;   maxrank[7]=50;   compmaxrank[ 7]=150; appdata[ 7]="--ss"
nrows[8]=2295000;   nb[8]=4500;     acc[8]=8;   maxrank[8]=50;   compmaxrank[ 8]=150; appdata[ 8]="--ss"
nrows[9]=3510000;   nb[9]=2700;     acc[9]=8;   maxrank[9]=50;   compmaxrank[ 9]=150; appdata[ 9]="--ss"
nrows[10]=3510000;  nb[10]=3000;    acc[10]=8;  maxrank[10]=50;  compmaxrank[10]=150; appdata[10]="--ss"
nrows[11]=3510000;  nb[11]=3375;    acc[11]=8;  maxrank[11]=50;  compmaxrank[11]=150; appdata[11]="--ss"
nrows[12]=3510000;  nb[12]=4500;    acc[12]=8;  maxrank[12]=50;  compmaxrank[12]=150; appdata[12]="--ss"
nrows[13]=4725000;  nb[13]=2700;    acc[13]=8;  maxrank[13]=50;  compmaxrank[13]=150; appdata[13]="--ss"
nrows[14]=4725000;  nb[14]=3000;    acc[14]=8;  maxrank[14]=50;  compmaxrank[14]=150; appdata[14]="--ss"
nrows[15]=4725000;  nb[15]=3375;    acc[15]=8;  maxrank[15]=50;  compmaxrank[15]=150; appdata[15]="--ss"
nrows[16]=4725000;  nb[16]=4500;    acc[16]=8;  maxrank[16]=50;  compmaxrank[16]=150; appdata[16]="--ss"
nrows[17]=5940000;  nb[17]=2700;    acc[17]=8;  maxrank[17]=50;  compmaxrank[17]=150; appdata[17]="--ss"
nrows[18]=5940000;  nb[18]=3000;    acc[18]=8;  maxrank[18]=50;  compmaxrank[18]=150; appdata[18]="--ss"
nrows[19]=5940000;  nb[19]=3375;    acc[19]=8;  maxrank[19]=50;  compmaxrank[19]=150; appdata[19]="--ss"
nrows[20]=5940000;  nb[20]=4500;    acc[20]=8;  maxrank[20]=50;  compmaxrank[20]=150; appdata[20]="--ss"
nrows[21]=8100000;  nb[21]=2700;    acc[21]=8;  maxrank[21]=50;  compmaxrank[21]=150; appdata[21]="--ss"
nrows[22]=8100000;  nb[22]=3000;    acc[22]=8;  maxrank[22]=50;  compmaxrank[22]=150; appdata[22]="--ss"
nrows[23]=8100000;  nb[23]=3375;    acc[23]=8;  maxrank[23]=50;  compmaxrank[23]=150; appdata[23]="--ss"
nrows[24]=8100000;  nb[24]=4500;    acc[24]=8;  maxrank[24]=50;  compmaxrank[24]=150; appdata[24]="--ss"
nrows[25]=10800000; nb[25]=2700;    acc[25]=8;  maxrank[25]=50;  compmaxrank[25]=150; appdata[25]="--ss"
nrows[26]=10800000; nb[26]=3000;    acc[26]=8;  maxrank[26]=50;  compmaxrank[26]=150; appdata[26]="--ss"
nrows[27]=10800000; nb[27]=3375;    acc[27]=8;  maxrank[27]=50;  compmaxrank[27]=150; appdata[27]="--ss"
nrows[28]=10800000; nb[28]=4500;    acc[28]=8;  maxrank[28]=50;  compmaxrank[28]=150; appdata[28]="--ss"


#nrows[1]=270000;   nb[1]=1350;     acc[1]=8;   maxrank[1]=50;   compmaxrank[ 1]=150; appdata[ 1]="--ss" 
#nrows[2]=270000;   nb[2]=2700;     acc[2]=8;   maxrank[2]=50;   compmaxrank[ 2]=150; appdata[ 2]="--ss"
nrows[5]=540000;   nb[5]=1350;     acc[5]=8;   maxrank[5]=50;   compmaxrank[ 5]=150; appdata[ 5]="--ss"
nrows[6]=540000;   nb[6]=2700;     acc[6]=8;   maxrank[6]=50;   compmaxrank[ 6]=150; appdata[ 6]="--ss"
allcaseids[2]="1 2 5 6"

allcaseids[4]="1 2 5 6"
allcaseids[4]="1 2"

allcaseids[8]="1 2 5 6"

allcaseids[16]="`seq 1 12`"
allcaseids[16]="5 6"

allcaseids[32]="`seq 1 16`"
#allcaseids[64]="`seq 1 15`"
#allcaseids[128]="`seq 1 20`"
#allcaseids[256]="`seq 17 24`"
#allcaseids[512]="`seq 21 28`"
#allcaseids[16]="`seq 1 28`"; 
#allcaseids[32]="`seq 1 28`"; 
#allcaseids[64]="`seq 1 28`"; 
#allcaseids[128]="`seq 1 28`"; 
#allcaseids[256]="`seq 1 28`"; 
#allcaseids[512]="`seq 1 28`"; 
#allcaseids[1024]="`seq 1 28`"; 
nprocs="2 4 8 16 32"
nprocs="8 16 32"
nprocs="2 4 8"
nprocs="4 16"
prog="hic"
fi



step=1
timelimit="02:00:00"
#_compmaxrank=70 # for 216000<= <=594000 maxrank=50
#_compmaxrank=100 #for 81000<= <216000 maxrank=100
#_compmaxrank=150 #for 54000 maxrank=100
note="Hicma beta=0.01 $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "

#prog="hic"
#prog="cham"

