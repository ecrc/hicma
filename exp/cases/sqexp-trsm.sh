timelimit="03:00:00"
sizes="big"
sizes="small"
if [ "$sizes" == "small"  ];then
#nrows[1]=1600;  nb[1]=1350; acc[1]=6;   maxrank[1]=200; nrhs[1]=400
nrows[1]=27000; nb[1]=2700; acc[1]=6;   maxrank[1]=1350;    nrhs[1]=2700;
nrows[2]=54000; nb[2]=2700; acc[2]=6;   maxrank[2]=1350;    nrhs[2]=2700;
nrows[3]=81000; nb[3]=2700; acc[3]=6;   maxrank[3]=1350;    nrhs[3]=2700;
nrows[4]=108000;    nb[4]=2700; acc[4]=6;   maxrank[4]=1350;    nrhs[4]=2700;
nrows[5]=135000;    nb[5]=2700; acc[5]=6;   maxrank[5]=1350;    nrhs[5]=2700;
nrows[6]=162000;    nb[6]=2700; acc[6]=6;   maxrank[6]=1350;    nrhs[6]=2700;
nrows[7]=189000;    nb[7]=2700; acc[7]=6;   maxrank[7]=1350;    nrhs[7]=2700;
nrows[8]=216000;    nb[8]=2700; acc[8]=6;   maxrank[8]=1350;    nrhs[8]=2700;
nrows[9]=243000;    nb[9]=2700; acc[9]=6;   maxrank[9]=1350;    nrhs[9]=2700;
nrows[10]=270000;   nb[10]=2700;    acc[10]=6;  maxrank[10]=1350;   nrhs[10]=2700;
nrows[11]=297000;   nb[11]=2700;    acc[11]=6;  maxrank[11]=1350;   nrhs[11]=2700;
nrows[12]=324000;   nb[12]=2700;    acc[12]=6;  maxrank[12]=1350;   nrhs[12]=2700;
nrows[13]=351000;   nb[13]=2700;    acc[13]=6;  maxrank[13]=1350;   nrhs[13]=2700;
nrows[14]=378000;   nb[14]=2700;    acc[14]=6;  maxrank[14]=1350;   nrhs[14]=2700;
nrows[15]=405000;   nb[15]=2700;    acc[15]=6;  maxrank[15]=1350;   nrhs[15]=2700;
nrows[16]=432000;   nb[16]=2700;    acc[16]=6;  maxrank[16]=1350;   nrhs[16]=2700;
nrows[17]=459000;   nb[17]=2700;    acc[17]=6;  maxrank[17]=1350;   nrhs[17]=2700;
nrows[18]=486000;   nb[18]=2700;    acc[18]=6;  maxrank[18]=1350;   nrhs[18]=2700;
nrows[19]=513000;   nb[19]=2700;    acc[19]=6;  maxrank[19]=1350;   nrhs[19]=2700;
nrows[20]=540000;   nb[20]=2700;    acc[20]=6;  maxrank[20]=1350;   nrhs[20]=2700;
nrows[21]=567000;   nb[21]=2700;    acc[21]=6;  maxrank[21]=1350;   nrhs[21]=2700;
nrows[22]=594000;   nb[22]=2700;    acc[22]=6;  maxrank[22]=1350;   nrhs[22]=2700;
allcaseids[16]="`seq 1 22`"
allcaseids[1]="1"
nprocs="1"
#allcaseids[1]="1";nprocs="1";que=debug
#allcaseids[2]="1";nprocs="2";que=debug; timelimit="00:30:00" #TODO
#allcaseids[2]="1";nprocs="2";
else
    :
fi


op=posv
step=1
_appdata="--ss"; 
_compmaxrank=1350 
note="Hicma trsm $_appdata - $sizes - $_wavek - $timelimit - $_compmaxrank "


