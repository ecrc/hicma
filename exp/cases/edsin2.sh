sizes="big"
sizes="small"
if [ "$sizes" == "small"  ];then
nrows[1]=54000; nb[1]=2700; acc[1]=8;   maxrank[1]=106;
nrows[2]=81000; nb[2]=2700; acc[2]=8;   maxrank[2]=106;
nrows[3]=108000;    nb[3]=2700; acc[3]=8;   maxrank[3]=106;
nrows[4]=135000;    nb[4]=2700; acc[4]=8;   maxrank[4]=106;
nrows[5]=162000;    nb[5]=2700; acc[5]=8;   maxrank[5]=106;
nrows[6]=189000;    nb[6]=2700; acc[6]=8;   maxrank[6]=106;
nrows[7]=216000;    nb[7]=2700; acc[7]=8;   maxrank[7]=106;
nrows[8]=243000;    nb[8]=2700; acc[8]=8;   maxrank[8]=106;
nrows[9]=270000;    nb[9]=2700; acc[9]=8;   maxrank[9]=106;
nrows[10]=297000;   nb[10]=2700;    acc[10]=8;  maxrank[10]=106;
nrows[11]=324000;   nb[11]=2700;    acc[11]=8;  maxrank[11]=106;
nrows[12]=351000;   nb[12]=2700;    acc[12]=8;  maxrank[12]=106;
nrows[13]=378000;   nb[13]=2700;    acc[13]=8;  maxrank[13]=106;
nrows[14]=405000;   nb[14]=2700;    acc[14]=8;  maxrank[14]=106;
nrows[15]=432000;   nb[15]=2700;    acc[15]=8;  maxrank[15]=106;
nrows[16]=459000;   nb[16]=2700;    acc[16]=8;  maxrank[16]=106;
nrows[17]=486000;   nb[17]=2700;    acc[17]=8;  maxrank[17]=106;
nrows[18]=513000;   nb[18]=2700;    acc[18]=8;  maxrank[18]=106;
nrows[19]=540000;   nb[19]=2700;    acc[19]=8;  maxrank[19]=106;
nrows[20]=567000;   nb[20]=2700;    acc[20]=8;  maxrank[20]=106;
nrows[21]=594000;   nb[21]=2700;    acc[21]=8;  maxrank[21]=106;
allcaseids[16]="`seq 1 21`"
allcaseids[32]="`seq 1 21`"
allcaseids[64]="`seq 1 21`"
allcaseids[128]="`seq 1 21`"
allcaseids[256]="`seq 1 21`"
nprocs="16 32 64 128 256"
else
nrows[1]=1080000;   nb[1]=2700; acc[1]=8;   maxrank[1]=100;
nrows[2]=1080000;   nb[2]=3000; acc[2]=8;   maxrank[2]=100;
nrows[3]=1080000;   nb[3]=3375; acc[3]=8;   maxrank[3]=100;
nrows[4]=1080000;   nb[4]=4500; acc[4]=8;   maxrank[4]=100;
nrows[5]=2295000;   nb[5]=2700; acc[5]=8;   maxrank[5]=100;
nrows[6]=2295000;   nb[6]=3000; acc[6]=8;   maxrank[6]=100;
nrows[7]=2295000;   nb[7]=3375; acc[7]=8;   maxrank[7]=100;
nrows[8]=2295000;   nb[8]=4500; acc[8]=8;   maxrank[8]=100;
nrows[9]=3510000;   nb[9]=2700; acc[9]=8;   maxrank[9]=100;
nrows[10]=3510000;  nb[10]=3000;    acc[10]=8;  maxrank[10]=100;
nrows[11]=3510000;  nb[11]=3375;    acc[11]=8;  maxrank[11]=100;
nrows[12]=3510000;  nb[12]=4500;    acc[12]=8;  maxrank[12]=100;
nrows[13]=4725000;  nb[13]=2700;    acc[13]=8;  maxrank[13]=100;
nrows[14]=4725000;  nb[14]=3000;    acc[14]=8;  maxrank[14]=100;
nrows[15]=4725000;  nb[15]=3375;    acc[15]=8;  maxrank[15]=100;
nrows[16]=4725000;  nb[16]=4500;    acc[16]=8;  maxrank[16]=100;
nrows[17]=5940000;  nb[17]=2700;    acc[17]=8;  maxrank[17]=100;
nrows[18]=5940000;  nb[18]=3000;    acc[18]=8;  maxrank[18]=100;
nrows[19]=5940000;  nb[19]=3375;    acc[19]=8;  maxrank[19]=100;
nrows[20]=5940000;  nb[20]=4500;    acc[20]=8;  maxrank[20]=100;
nrows[21]=8100000;  nb[21]=2700;    acc[21]=8;  maxrank[21]=100;
nrows[22]=8100000;  nb[22]=3000;    acc[22]=8;  maxrank[22]=100;
nrows[23]=8100000;  nb[23]=3375;    acc[23]=8;  maxrank[23]=100;
nrows[24]=8100000;  nb[24]=4500;    acc[24]=8;  maxrank[24]=100;
nrows[25]=10800000; nb[25]=2700;    acc[25]=8;  maxrank[25]=100;
nrows[26]=10800000; nb[26]=3000;    acc[26]=8;  maxrank[26]=100;
nrows[27]=10800000; nb[27]=3375;    acc[27]=8;  maxrank[27]=100;
nrows[28]=10800000; nb[28]=4500;    acc[28]=8;  maxrank[28]=100;
allcaseids[16]="`seq 1 8`"
allcaseids[32]="`seq 1 12`"
allcaseids[64]="`seq 1 15`"
allcaseids[128]="`seq 1 20`"
allcaseids[256]="`seq 17 24`"
allcaseids[512]="`seq 21 28`"
nprocs="16 32 64 128 256 512"
fi



step=1
_appdata="--edsin"; timelimit="00:25:00"
_appdata="--edsin"; _wavek=100; _compmaxrank=250;
note="Hicma only dense matrix $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "


