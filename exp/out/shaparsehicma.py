#!/usr/bin/python
import sys
import traceback
import socket
import glob
hn=socket.gethostname()
hostname=None
if hn.startswith('cdl'):
    hostname="shaheen"
if hn.startswith('xci'):
    hostname="isambard"
#import platform
#print(platform.node(), hn)
#print(hn)
#verbose=True
verbose=True
verbose=False

if hostname == "isambard":
    ROOT="/home/ri-kakbudak/kexphicma/ipy/hicma-shaheen/"
    ROOT=""
    od="/home/ri-kakbudak/hicma-dev/exp/out"
    jobsfile=[
            #'/lustre/home/ri-kakbudak/hicma-dev/exp/jobids/2019-09-05-hicma-stat-isambard-1.txt'
            '/lustre/home/ri-kakbudak/hicma-dev/exp/jobids/2019-09-10-hicma-stat-isambard-1.txt'
            ]
elif hostname == "shaheen":
    ROOT="/lustre/project/k1205/akbudak/kexphicma/ipy/hicma-shaheen/"
    od="/lustre/project/k1205/akbudak/hicma-dev/exp/out"
    jobsfile="jobshicma.txt"
    ROOT=""
    #IPDPS, 2017
    jobsfile="/lustre/project/k1205/akbudak/hicma-dev/jobs-10-17-1.txt"
    jobsfile="/lustre/project/k1205/akbudak/hicma-dev/exp/jobids/2017-10-18-hicma-10M-1.txt"
    jobsfile='/lustre/project/k1205/akbudak/hicma-dev/exp/jobids/2017-10-18-hicma-10M-2.txt'
    #jobsfile='/lustre/project/k1205/akbudak/hicma-dev/exp/jobids/2017-10-19-hicma-50K500K-1.txt'
    jobsfile='/lustre/project/k1205/akbudak/hicma-dev/exp/jobids/all-hicma-10M.txt'
    #jobsfile='/lustre/project/k1205/akbudak/hicma-dev/exp/jobids/2017-10-21-hicma-acc-1.txt'
    #jobsfile='/lustre/project/k1205/akbudak/hicma-dev/exp/jobids/2017-10-21-hicma-acc-nocomp-1.txt'
    jobsfile='/lustre/project/k1205/akbudak/hicma-dev/exp/jobids/2017-10-21-hicma-compare-1.txt' # to compare with KNL and Skylake
    #jobsfile='/lustre/project/k1205/akbudak/hicma-dev/exp/jobids/2017-10-22-hicma-acc-mem-1.txt' # different acc for 16 nodes
    jobsfile='/lustre/project/k1205/akbudak/hicma-dev/exp/jobids/2017-11-02-hicma-more-points-BIG-1.txt'
    jobsfile='/lustre/project/k1205/akbudak/hicma-dev/exp/jobids/2017-11-02-hicma-more-points-1.txt'

    #EuroPar runs were started in January 22, 2018
    jobsfile='/lustre/project/k1205/akbudak/hicma-dev/exp/jobids/2018-01-22-hicma-edsin-acc-1.txt'
    jobsfile=['/lustre/project/k1205/akbudak/hicma-dev/exp/jobids/2018-01-23-hicma-edsin-acc-custommaxrk-1.txt']
    ROOT='/lustre/project/k1205/akbudak/hicma-dev/'
    jobsfile=['exp/jobids/2018-01-23-hicma-edsin-acc-custommaxrk-1.txt']
    Xjobsfile=['exp/jobids/2018-02-10-hicma-geostat-16-1.txt'
            ,'exp/jobids/2018-02-10-hicma-edsin-16-1.txt'
            ,'exp/jobids/2018-02-10-hicma-rnd-16-1.txt'
            ]
    Xjobsfile=[
            'exp/jobids/2018-02-12-hicma-edsin-1.txt'
            ,'exp/jobids/2018-02-12-hicma-sqexp-1.txt'
            ,'exp/jobids/2018-02-12-hicma-matern-1.txt'
            ,'exp/jobids/2018-02-10-hicma-rnd-16-1.txt'
            ]
    jobsfile=['exp/jobids/2018-02-12-hicma-edsin-scale-1.txt','exp/jobids/2018-02-16-hicma-edsin-scale-1.txt','exp/jobids/2018-02-16-hicma-edsin-scale-2.txt']
    jobsfile=['exp/jobids/2018-02-17-hicma-only-prob-1.txt']
    Xjobsfile=[
            #'exp/jobids/ji-2018-02-12-turbo_off.txt'
            'exp/jobids/ji-2018-02-12-turbo_on_small_block_sizes.txt',
            'exp/jobids/ji-2018-02-12-turbo_on.txt'
            ]
    jobsfile=['exp/jobids/2018-02-17-hicma-trsm-4.txt']
    jobsfile=['exp/jobids/2018-03-19-hicma-edsin-big-256-afterupgrade-1.txt']
    jobsfile=['exp/jobids/2018-04-15-hicma-threads-2.txt']
    #jobsfile=['exp/jobids/2018-04-15-hicma-threads-smalltile-1.txt']
    jobsfile=['exp/jobids/2018-04-17-hicma-threads-updatecham-1.txt']
    jobsfile=['exp/jobids/2018-05-07-hicma-threads-scheds-1.txt']
else:
    ROOT="/Users/akbudak/kexphicma/ipy/hicma-shaheen/"
    ROOT="/Users/kadirakbudak/hicma-dev/exp/out/"
    ROOT="./"
    od="/Users/akbudak/kexphicma/ipy/hicma-shaheen"
    od=None
    jobsfile=["jobshicma.txt"]
    if verbose is True:
        print(jobsfile)
jobids=[]
for fn in jobsfile:
    if verbose is True:
        print("Jobsfile:", fn)
    with open(ROOT+fn) as f:
        jobidsAll = f.readlines()
        for x in jobidsAll:
            if x.startswith('#STOP'):
                break
            if not x.startswith('#'):
                jobids.append(x)
if verbose is True:
    print("Jobids:", jobids)
#jobids=['4315255'] #test file to get time compress
#jobids=[jobids[0]]
jobids = [x.strip() for x in jobids]
for id in jobids:
    if "STOP" in id:
        break
    try:
        if od is not None:
            files = [f for f in glob.glob( od+'/*'+id+"*")]
            resfile=files[0]
        else:
            resfile=id
        if verbose is True:
            print(resfile)
        with open(resfile) as f:
            reslines = f.readlines();
            #print(resfile)
            # print(reslines)
        for line in list(reslines):
            if line.startswith('<'):
                reslines.remove(line)
    except:
        sys.stdout.flush()
        if verbose is True:
            sys.stderr.write("File not found: " + resfile + "\n")
        sys.stderr.flush()
        continue
    nlines=len(reslines)
    iline = 0;
    if verbose is True:
        print("Starting to parse file: ",resfile)
        #print(reslines)
    while iline < nlines:
        try:
            line = reslines[iline]
            #if verbose is True:
            #    print("Processing:", line)
            if "SCHED:" in line and "!BEGIN" in line:
                if verbose is True:
                    print("Found SCHED and !BEGIN in line:", line)
                i = iline
                seek_next_begin=False
                while i < iline+6 and i < nlines:
                    if "!END" in reslines[i]:
                        iline = i
                        seek_next_begin=True
                        if verbose is True:
                            print("Will start next result because found !END in line:", line)
                        break
                    i += 1
                if seek_next_begin is True:
                    continue
                schedline = line
                iline += 2
                i = iline
                while iline < nlines:
                    if "hicma-dev" in reslines[iline]:
                        break
                    iline += 1
                if iline >= nlines:
                    break

                found_hicmadev=False
                i=iline
                __line = reslines[i]
                while i < nlines:
                    if verbose is True:
                        print("Does this line |"+__line+"| contains hicma-dev/build? ", "hicma-dev/build" in __line)
                    if "hicma-dev/build" in __line:
                        found_hicmadev=True
                        break
                    i += 1
                    __line = reslines[i]
                if found_hicmadev is True:
                    iline = i
                code="hic"
                if "chameleon" in reslines[iline]:
                    code="cham"
                    iavgrk, iminrk, imaxrk, favgrk, fminrk, fmaxrk, tproblem, tcompress, fixedrank, fixedacc, wavek, shmaxrk, shprob, shdecay = [0] * 14
                    if verbose is True:
                        print('code switched to cham')

                if "hicma-dev" not in reslines[iline]:
                    continue
                try:
                    if verbose is True:
                        print("Splitted binary:", reslines[iline].split(' ') )
                    binary = reslines[iline].split(' ')[2].strip()
                    if code is "cham":
                        binary = reslines[iline].split(' ')[2].strip()
                except Exception as e:
                    print('Error white getting binary at line', iline, binary, resfile)
                    print(e)
                    sys.exit()
                if verbose is True:
                    print("Binary:",binary)
                iline += 1
                if iline >= nlines:
                    break
                threadline = reslines[iline]
                try:
                    nthreads = int(threadline.split(':')[1])
                except:
                    print('Error in line', iline, threadline, resfile)
                    sys.exit()
                iline += 3
                if iline >= nlines:
                    break
                mpiline = reslines[iline]
                if verbose is True:
                    print("MPI line:", iline, resfile, mpiline)
                try:
                    nmpiP, nmpiQ = mpiline.split(':')[1].split('x')
                except:
                    #print "Error occured in reading", resfile
                    #print "SCHEDLINE:", schedline
                    #print "MPILINE:", mpiline
                    #sys.exit
                    nmpiP = 1
                    nmpiQ = 1
                    iline -= 2


                nmpiP = int(nmpiP)
                nmpiQ = int(nmpiQ)
                nmpi = nmpiP * nmpiQ
                if verbose is True:
                    print(str(nmpiP)+"x"+str(nmpiQ)+"="+str(nmpi));

                if code is "hic":
                    iline += 5
                    fixedrank=int(reslines[iline].split(':')[1])
                    iline += 1
                    fixedacc=float(reslines[iline].split(':')[1])
                    iline += 2
                    wavek=float(reslines[iline].split(':')[1])
                    if verbose is True:
                        print("wavek: ", wavek)
                    iline += 1
                    shmaxrk=int(reslines[iline].split(':')[1])
                    iline += 1
                    shprob=int(reslines[iline].split(':')[1])
                    iline += 1
                    shdecay=float(reslines[iline].split(':')[1])
                    iline += 3
                    if iline >= nlines:
                        break
                    resline = reslines[iline]
                    tproblem = float('nan')
                    tcompress = float('nan')
                    if verbose is True:
                        print(resline)
                    found_tproblem=False
                    i=iline
                    __line = reslines[i]
                    abort=False
                    while i < nlines:
                        if verbose is True:
                            print("Is Tgeneratecompress1 in:", __line)
                        if "Tgeneratecompress1:" in __line:
                            found_tproblem=True
                            break
                        i += 1
                        if "!END" in __line:
                            abort = True
                            break
                        __line = reslines[i]
                    if abort is True:
                        print("==========================================ABORTING SINCE !END IS ENCOUNTERED EARLY", schedline, file=sys.stderr)
                        continue
                    if found_tproblem is True:
                        iline = i
                        resline = reslines[iline]
                    if "Tgeneratecompress1:" in resline:
                        tproblem=float(resline.split(':')[1])
                        iline+=1
                        if iline >= nlines:
                            break
                        resline = reslines[iline]
                        found_tcompress=False
                        i=iline
                        while i < nlines:
                            __line = reslines[i]
                            if verbose is True:
                                print("Is Tgeneratecompress2 in:", __line)
                            if "Tgeneratecompress2:" in __line:
                                found_tcompress=True
                                break
                            i += 1
                        if found_tcompress is True:
                            iline = i
                            resline = reslines[iline]
                        #print resline
                        while "Tgeneratecompress2:" not in resline:
                            iline+=1
                            if iline >= nlines:
                                break
                            resline = reslines[iline]
                            #print resfile
                            #print iline
                            #print "loop:", resline
                        if "Tgeneratecompress2:" not in resline:
                            break
                        linetcompress=resline.split(':')
                        if len(linetcompress) < 2:
                            break
                        tcompress=float(linetcompress[1])
                        iline+=1
                        if iline >= nlines:
                            break
                        resline = reslines[iline]
                        found_initialranks=False
                        i=iline
                        abort = False
                        while i < nlines:
                            __line = reslines[i]
                            if verbose is True:
                                print("Is initial_ranks in:", __line)
                            if "!END" in __line:
                                abort = True
                                break
                            if "initial_ranks:" in __line:
                                found_initialranks=True
                                break
                            i += 1
                        if abort is True:
                            print("==========================================ABORTING SINCE !END IS ENCOUNTERED EARLY", schedline, file=sys.stderr)
                            continue
                        if found_initialranks is True:
                            iline = i
                            resline = reslines[iline]
                        while "initial_ranks" not in resline and "+-" not in resline and iline < nlines:
                            iline+=1
                            if iline >= nlines:
                                break
                            resline = reslines[iline]
                            if "!END" in resline:
                                abort = True
                                break
                        if abort is True:
                            print("==========================================ABORTING SINCE !END IS ENCOUNTERED EARLY", schedline, file=sys.stderr)
                            continue

                    rank_stats_printed = False
                    iavgrk = float('nan')
                    iminrk = float('nan')
                    imaxrk = float('nan')
                    favgrk = float('nan')
                    fminrk = float('nan')
                    fmaxrk = float('nan')
                    if "initial_ranks:" in resline:
                        iavgrk, iminrk, imaxrk = [float(x.split(':')[-1]) for x in resline.split(' ')]
                        if verbose is True:
                            print("initial ranks:", iavgrk, iminrk, imaxrk)
                        rank_stats_printed = True
                    if rank_stats_printed is True:
                        iline += 1
                        if iline >= nlines:
                            break
                        resline = reslines[iline]
                        found_finalranks=False
                        i=iline
                        while i < nlines:
                            __line = reslines[i]
                            if verbose is True:
                                print("Is final_ranks in:", __line)
                            if "final_ranks:" in __line:
                                found_finalranks=True
                                break
                            if "!END" in __line:
                                abort = True
                                break
                            i += 1
                        if abort is True:
                            if verbose is True:
                                print("==========================================ABORTING SINCE !END IS ENCOUNTERED EARLY", schedline, resfile, file=sys.stderr)
                            continue
                        if found_finalranks is True:
                            iline = i
                            resline = reslines[iline]

                        if "final_ranks:" in resline:
                            favgrk, fminrk, fmaxrk = [float(x.split(':')[-1]) for x in resline.split(' ')]
                            if verbose is True:
                                print("final ranks:", favgrk, fminrk, fmaxrk)
                        iline += 1
                        if iline >= nlines:
                            break
                        resline = reslines[iline]
                if code is "cham":
                    abort=False
                    found_tproblem=False
                    i=iline
                    while i < nlines:
                        if verbose is True:
                            print("Is Tgenerate in:", __line)
                        if "Tgenerate:" in __line:
                            found_tproblem=True
                            break
                        i += 1
                        if "!END" in __line:
                            abort = True
                            break
                        __line = reslines[i]
                    if abort is True:
                        print("==========================================ABORTING SINCE !END IS ENCOUNTERED EARLY", schedline, file=sys.stderr)
                        continue
                    if found_tproblem is True:
                        iline = i
                        resline = reslines[iline]
                        linetcompress=resline.split(':')
                        if len(linetcompress) < 2:
                            break
                        tproblem=float(linetcompress[1])
                        iline+=1
                found_chamres=False
                i=iline
                while i < nlines:
                    __line = reslines[i]
                    if verbose is True:
                        print("Is +- in:", __line)
                    if "+-" in __line:
                        found_chamres=True
                        break
                    i += 1
                if found_chamres is True:
                    iline = i
                    resline = reslines[iline]
                if "+-" in resline:
                    #print schedline
                    dsched = {}
                    for x in schedline.split(' '):
                        if ":" in x:
                            kvs = x.split(':')
                            dsched[kvs[0].strip()] = kvs[1].strip()
                    if 'MAXSUB' not in dsched:
                        dsched["MAXSUB"]=0
                        dsched['MINSUB']=0
                    #print dres['SCHED']
                    ares = " ".join(resline.split()).split(' ')
                    hn=hostname
                    if hostname is None: #if hostname is not set, try to grep it from name of the RAW result file
                        for name in ['shihab','jasmine','flamingo']:
                            if name in resfile:
                                hn=name
                                break
                    print(code, dsched["SCHED"], ares[0], dsched['MB'], nmpi, nmpiP, nmpiQ, ares[3], ares[4], dsched['NB'], dsched["MAXSUB"], dsched['MINSUB'], nthreads, id, iavgrk, iminrk, imaxrk, favgrk, fminrk, fmaxrk, tproblem, tcompress, fixedrank, fixedacc, wavek, hn, binary, shmaxrk, shprob, shdecay)
        except Exception as e:
            #pass
            print('File:', resfile)
            print(traceback.format_exc())
            raise(e)
        iline+=1
    #break #TODO
sys.exit()

# ids=`cat $ROOT/jobshicma.txt|grep -v "#"`
# echo Output of these jobs will be processed: $ids 1>&2
# od=/lustre/project/k1205/akbudak/hicma-dev/exp/out
# for id in $ids; do
    # f=$od/$id
    # echo processing $f  1>&2
    # np=`cat $f|grep "# Nb mpi:"|head -1|awk -F':' '{print $2}'`
    # cat $f|grep "+-"|sed -e "s/^/hicma $np/"
    # continue
    # while IFS='' read -r line || [[ -n "$line" ]]; do
        # #echo "Text read from file: $line"
        # sched=`echo $line|grep  "SCHED:"|awk -F':| ' '{print $21}'`
        # if [ ! -z "$sched" ]; then
            # echo $sched
        # fi
        # mb=`echo $line|grep  "MB:"|awk -F':| ' '{print $6}'`
        # if [ ! -z "$mb" ]; then
            # echo $mb
        # fi
    # done < "$f"
    # #sched=`cat $f|grep -2 "SCHED:"|head -1|awk -F':| ' '{print $21}'`
    # #res=`cat $f|grep -2 "SCHED:"`
    # #echo "Sched:"$sched
    # #echo $res
    # break
# done
