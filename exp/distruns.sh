#!/bin/bash
dry="-"
nprocs="16 32 64 128 256"
nprocs="64"
cases=
trace="-"
op="potrf"
que="workq"; timelimit="06:00:00"
#que="debug"
if [ $# -eq 1 -o $# -eq 2 ]; then
    cases=$1
    if [ ! -f "$cases" ];then
        echo "Cases file does not exist: $cases" 
        exit
    else
        . $cases
#allcaseids[16]="`seq 1 28`"
#allcaseids[32]="`seq 1 28`"
#allcaseids[64]="`seq 1 28`"
#allcaseids[128]="`seq 1 28`"
#allcaseids[256]="`seq 1 28`"
#allcaseids[512]="`seq 1 28`"
#nprocs="16 32 64 128 256 512"
##_appdata="--edsin"; _wavek=25; _compmaxrank=100; tlmxrk=30
#_appdata="--edsin"; _wavek=50; _compmaxrank=150; tlmxrk=50
#step=1
##_appdata="--edsin"; _wavek=100; _compmaxrank=300
#timelimit="03:00:00"
#note="Hicma $_appdata - $sizes - wavek:$_wavek - timelimit:$timelimit - compmaxrank:$_compmaxrank - tile max rank:$tlmxrk "
#nrows[1]=1080000;   nb[1]=2700; acc[1]=8;   maxrank[1]=$tlmxrk;
#nrows[2]=1080000;   nb[2]=3000; acc[2]=8;   maxrank[2]=$tlmxrk;
#nrows[3]=1080000;   nb[3]=3375; acc[3]=8;   maxrank[3]=$tlmxrk;
#nrows[4]=1080000;   nb[4]=4500; acc[4]=8;   maxrank[4]=$tlmxrk;
#nrows[5]=2295000;   nb[5]=2700; acc[5]=8;   maxrank[5]=$tlmxrk;
#nrows[6]=2295000;   nb[6]=3000; acc[6]=8;   maxrank[6]=$tlmxrk;
#nrows[7]=2295000;   nb[7]=3375; acc[7]=8;   maxrank[7]=$tlmxrk;
#nrows[8]=2295000;   nb[8]=4500; acc[8]=8;   maxrank[8]=$tlmxrk;
#nrows[9]=3510000;   nb[9]=2700; acc[9]=8;   maxrank[9]=$tlmxrk;
#nrows[10]=3510000;  nb[10]=3000;    acc[10]=8;  maxrank[10]=$tlmxrk;
#nrows[11]=3510000;  nb[11]=3375;    acc[11]=8;  maxrank[11]=$tlmxrk;
#nrows[12]=3510000;  nb[12]=4500;    acc[12]=8;  maxrank[12]=$tlmxrk;
#nrows[13]=4725000;  nb[13]=2700;    acc[13]=8;  maxrank[13]=$tlmxrk;
#nrows[14]=4725000;  nb[14]=3000;    acc[14]=8;  maxrank[14]=$tlmxrk;
#nrows[15]=4725000;  nb[15]=3375;    acc[15]=8;  maxrank[15]=$tlmxrk;
#nrows[16]=4725000;  nb[16]=4500;    acc[16]=8;  maxrank[16]=$tlmxrk;
#nrows[17]=5940000;  nb[17]=2700;    acc[17]=8;  maxrank[17]=$tlmxrk;
#nrows[18]=5940000;  nb[18]=3000;    acc[18]=8;  maxrank[18]=$tlmxrk;
#nrows[19]=5940000;  nb[19]=3375;    acc[19]=8;  maxrank[19]=$tlmxrk;
#nrows[20]=5940000;  nb[20]=4500;    acc[20]=8;  maxrank[20]=$tlmxrk;
#nrows[21]=8100000;  nb[21]=2700;    acc[21]=8;  maxrank[21]=$tlmxrk;
#nrows[22]=8100000;  nb[22]=3000;    acc[22]=8;  maxrank[22]=$tlmxrk;
#nrows[23]=8100000;  nb[23]=3375;    acc[23]=8;  maxrank[23]=$tlmxrk;
#nrows[24]=8100000;  nb[24]=4500;    acc[24]=8;  maxrank[24]=$tlmxrk;
#nrows[25]=10800000; nb[25]=2700;    acc[25]=8;  maxrank[25]=$tlmxrk;
#nrows[26]=10800000; nb[26]=3000;    acc[26]=8;  maxrank[26]=$tlmxrk;
#nrows[27]=10800000; nb[27]=3375;    acc[27]=8;  maxrank[27]=$tlmxrk;
#nrows[28]=10800000; nb[28]=4500;    acc[28]=8;  maxrank[28]=$tlmxrk;
    fi
    if [ $# -eq 2 ]; then
        dry=$2
    fi
else 
    echo "Usage: $0 casesFile [dry]" 1>&2
    exit
fi
#corespernode=40; numthreads=39 #skylake
corespernode=30; numthreads=31 #shaheen
corespernode=30; numthreads="1 2 4 8 16 31" #shaheen

#que="debug"; timelimit="00:30:00"; echo "DEBUG DEBUG DEBUG DEBUG DEBUG "
for nodes in $nprocs; do
    echo "#Number of nodes: $nodes ============================="
    if [ "$dry" == "dry" ]; then
        cmdbatch=""
    else
        cmdbatch="sbatch --parsable --nodes=$nodes --partition=$que --time=$timelimit --account=k1205 --job-name=hicma-$nodes "
    fi
    maxsubs[1]=500; minsubs[1]=250
    maxsubs[2]=1000; minsubs[2]=500
    maxsubs[3]=2000; minsubs[3]=1000
    maxsubs[4]=4000; minsubs[4]=2000
    maxsubs[5]=6000; minsubs[5]=3000
    #for isub in 1 2 3 4 5;do
    for isub in 0;do
        #maxsub=${maxsubs[isub]}
        #minsub=${minsubs[isub]}
        sched="prio eager random ws lws dm dmda dmdar dmdas dmdasd heft pheft peager" 
        sched="prio"

        #cases="exp/cases/cham32.txt";       startt=107; endt=108; step=200
        #ids="`seq 1 5`"
        #$cmdbatch exp/trial.sh $nodes $nodes 31 - cham - "$ids" $dry $cases "13:00:00" $que $sched $maxsub $minsub;continue


        prog="hic"
        #prog="cham"
        caseids=${allcaseids[$nodes]}      
        ncases=`echo "$caseids" | wc -w`    
        startt=0;   endt=$((ncases-1)); 
        echo -n "#`date` on $nodes nodes. $note - $prog - $cases "
        echo -n \"$caseids\" 
        echo " $ncases ($startt-$endt-$step)"

        ct=$startt
        while [ $ct -le $endt ]; do
            et=$((ct+step-1))
            if [ $et -gt $endt ]; then
                et=$endt
            fi
            #ids="`seq $ct $et`"
            arrcaseids=($caseids)
            ids=
            for icase in `seq $ct $et`; do
                ids="$ids ${arrcaseids[icase]}"
            done
            ct=$((ct+step))
            echo "#case ids: $ids"
            for nt in $numthreads; do
                corespernode=$nt
                maxsub=$((nodes*(corespernode+10)))
                minsub=$((nodes*corespernode))
                $cmdbatch exp/distmem.sh $nodes $nodes $nt $trace $prog - "$ids" $dry $cases $timelimit $que $sched $maxsub $minsub $op
            done
        done
    done
done

