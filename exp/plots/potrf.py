#!/usr/bin/env python3

#/Users/akbudak/anaconda3/bin/python

import re
import sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
if len(sys.argv) != 4:
    print("Usage: potrf.py hicma-cpu hicma-gpu mkl")
    print("Sample run:")
    print("\t./potrf.py ../out/cpu.txt ../out/gpu.txt ../out/mkl.txt")
    sys.exit(-1);

print("This is the name of the script: ", sys.argv[0])
print("Number of arguments: ", len(sys.argv))
print("The arguments are: " , str(sys.argv))

listfile_cpu=sys.argv[1]
listfile_gpu=sys.argv[2]
listfile_mkl=sys.argv[3]
listfiles={'cpu':listfile_cpu, 'gpu':listfile_gpu, 'mkl':listfile_mkl}

def process_raw_file(key, filename):
    print(filename)
    if key == 'cpu':
        df = process_cpu_raw_file(filename)
    elif key == 'gpu':
        df = process_gpu_raw_file(filename)
    elif key == 'mkl':
        df = process_mkl_raw_file(filename)
    return df

def process_cpu_raw_file(filename):
    print('processing cpu', filename)
    with open(filename) as f:
        lines = f.readlines()
    allres=[]
    iline = 0
    while iline < len(lines):
        line=lines[iline]
        if line.startswith("# MB:"):
            res={}
            res["mb"]=int(line.split(':')[1])
        if line.startswith("# fixed acc:"):
            res["acc"]=float(line.split(':')[1])
        if line.startswith("# shprob:"):
            intprob=int(line.split(':')[1])
            if intprob == 2:
                prob = 'sqexp'
                probndim = '2'
            elif intprob == 15:
                prob = 'exp'
                probndim = '2'
            res["KERNEL"]=prob
            res["PROB"]=probndim
        if line.startswith("# reorder inner products:"):
            res["rip"]=int(line.split(':')[1])
        for (pre,txt) in [
                ("i","Ainitial_ranks:"),
                ("f", "Cfinal_ranks:")]:
            if line.startswith(txt):
                res[pre+"avg"] = float(re.split(':| ',line)[2])
                res[pre+"min"] = int(re.split(':| ',line)[4])
                res[pre+"max"] = int(re.split(':| ',line)[6])
        if line.startswith("ReShg"):
            arr=re.split('\t', line);
            res["flop"]    =   int(arr[1])
            res["gflop"]   = float(arr[2])
            res["GF/s"] = float(arr[3])
            res["TCHOL"]    = float(arr[4])
            line=lines[iline+1]
            endoffile = False
            newresult = False
            while "+-" not in line:
                iline += 1
                if iline > len(lines):
                    endoffile = True
                    break
                if line.startswith("# MB:"):
                    newresult = True
                    break
                line = lines[iline+1]
            if endoffile is True:
                break
            if newresult is True:
                continue

            arr=line.strip().split(' ');
            #print(arr)
            res["N"]     =  int(arr[0])
            #res["N"]     = int(arr[1])
            #res["K"]     = int(arr[2])
            #res["time2"] = float(arr[3])
            #print(res)
            allres.append(res)
        iline+=1


    df=pd.DataFrame.from_dict(allres)
    df['Label']='HiCMA-CPU'
    df['HOST']='Vulture'
    df['NTHR']='0'
    df['STRG']='H'
    df['NB']=df['mb']
    return df


def process_gpu_raw_file(filename):
    print('processing gpu', filename)
    with open(filename) as f:
        lines = f.readlines()
        lines2=[]
        for l in lines:
            if 'R-STATIC' in l or 'MAGMA' in l or 'CUSOLVER' in l:
                ncols=len(l.split())
                lines2.append(l.split())
    print('ncols:', ncols, 'nresults:', len(lines2))
    colsv1=['Label','PROB','KERNEL','ARCH','N','NB','NTHR','ACC','STMXRK','GF/s','TCOMP','TCHOL','IAVG','IMAX','FAVG','FMAX']
    colsv2=['Label','PROB','KERNEL','STRG','ARCH','N','NB','NTHR','ACC','STMXRK','GF/s','TCOMP','TCHOL','IAVG','IMAX','FAVG','FMAX']
    try:
        df = pd.DataFrame(lines2, columns=colsv1)
        df['STRG']='H'
    except:
        df = pd.DataFrame(lines2, columns=colsv2)
    if 'v100' in filename:
        df['HOST']='V100'
    if 'a100' in filename:
        df['HOST']='A100'
    if 'ibexrome' in filename:
        df['HOST']='IbexRome'
        #print(df.to_string())
    #df = pd.read_csv(filename, delim_whitespace=True) #header=None,
    #print(df.to_string())
    return df

def process_mkl_raw_file(filename):
    print('processing mkl', filename)
    with open(filename) as f:
        lines = f.readlines()
        lines2=[]
        for l in lines:
            if l.startswith('R-MKL-DPOTRF'):
                lines2.append(l.split())
    df = pd.DataFrame(lines2, columns=['Label','N','TCHOL','GF/s'])
    df['Label']='MKL'
    df['HOST']='Vulture'
    df['PROB']='0'
    df['KERNEL']='0'
    df['NTHR']='0'
    df['STRG']='H'
    #print(df)
    return df


dfs=[]
for lf in listfiles:
    fl=listfiles[lf]
    if fl != "NA":
        print("reading list:", fl)
        with open(fl) as f:
            lines = f.readlines()
            for l in lines:
                if l.startswith('#'):
                    continue
                else:
                    df = process_raw_file(lf, '../out/'+l.strip())
                    dfs.append(df)
df = pd.concat(dfs)
df["TCHOL"] = pd.to_numeric(df["TCHOL"])
df["N"] = pd.to_numeric(df["N"])
df["NB"] = pd.to_numeric(df["NB"])

## block size vs total number of operations
#print(df[(df['N']==196560)&(df['Label']=='HiCMA-CPU')][['KERNEL','mb','gflop','TCHOL']])
#print(df[(df['HOST']=='IbexRome')&(df['Label']=='R-STATIC')][['KERNEL','N','NB','TCHOL']].sort_values(['KERNEL','N','NB']).to_string())
prettyval={'gflop':'#ops (gflop)','TCHOL':'Time (s)'}
for val in ['gflop','TCHOL']:
    for kernel in ['exp','sqexp']:
        fig, ax = plt.subplots()
        print(kernel)
        df5=df[(df['HOST']=='Vulture')&(df['KERNEL']==kernel)&(df['Label']=='HiCMA-CPU')][['KERNEL','N','NB',val]].sort_values(['KERNEL','N','NB'])
        #df5=df[(df['HOST']=='IbexRome')&(df['KERNEL']==kernel)&(df['Label']=='R-STATIC')][['KERNEL','N','NB','TCHOL']].sort_values(['KERNEL','N','NB'])
        print(df5.to_string())
        matrixsizes = np.sort(df5['N'].unique())
        for matrixsize in matrixsizes:
            ax = df5[(df5['N']==matrixsize)].plot(ax=ax, x='NB', y=val, label=matrixsize)
        plt.xlabel('NB')
        plt.ylabel(prettyval[val])
        ax.legend(bbox_to_anchor=(1.3, 1.00), title=kernel)
        fig.savefig(kernel+'-'+val+'-NB.pdf', bbox_inches='tight')
        fig.savefig(kernel+'-'+val+'-NB.png', bbox_inches='tight')
sys.exit()

## agg with a dict of functions
## https://stackoverflow.com/questions/47360510/pandas-groupby-and-aggregation-output-should-include-all-the-original-columns-i
df2 = df.groupby(['Label','HOST','N','PROB','KERNEL','NTHR','STRG'], as_index=False).agg([('TCHOL','min')]).reset_index()
df2.columns = df2.columns.droplevel(1)

df3=df2[['Label','HOST','N','NB','gflop','PROB','KERNEL','NTHR','STRG','TCHOL']]

#print(df3.to_string())
labels=df3['Label'].unique()
print(labels)
labels=['MKL', 'HiCMA-CPU', 'R-STATIC'] ## custom order
hosts=df3['HOST'].unique()
probndims=df3['PROB'].unique()
kernels=df3['KERNEL'].unique()
ngpus=df3['NTHR'].unique()
prettylabels={'R-STATIC':'HiCMA-GPU'}


def plot(hosts, labels, probndims, kernels, ngpus, strgs, plotname):
    fig, ax = plt.subplots()
    for host in hosts:
        for label in labels:
            for probndim in probndims:
                for kernel in kernels:
                    for ngpu in ngpus:
                        for strg in strgs:
                            df=df3[(df3['HOST']==host)&(df3['Label']==label)&(df3['PROB']==probndim)&(df3['KERNEL']==kernel)&(df3['NTHR']==ngpu)&(df3['STRG']==strg)]
                            nresults=len(df)
                            if nresults == 0: continue ## skip this result
                            print(label, probndim, kernel, nresults)
                            #print(df.columns)
                            #print(df[['Label', 'HOST', 'N', 'PROB', 'KERNEL', 'NTHR', 'STRG', 'TCHOL']])
                            try:
                                prettylabel=prettylabels[label]
                            except:
                                prettylabel=label
                            if label == 'MKL':
                                case=prettylabel
                            elif label == 'HiCMA-CPU':
                                case=prettylabel+"-"+str(probndim)+"D-"+kernel
                            elif ngpu == '1':
                                case=host+'-'+prettylabel+"-"+ngpu+" GPU-"+str(probndim)+"D-"+kernel+'-'+strg
                                case=host+'-'+prettylabel+"-"+str(probndim)+"D-"+kernel+'-'+strg
                            else:
                                case=host+'-'+prettylabel+"-"+ngpu+" GPUs-"+str(probndim)+"D-"+kernel+'-'+strg
                                case=host+'-'+prettylabel+"-"+str(probndim)+"D-"+kernel+'-'+strg
                            if label in ['MAGMA']:
                                case=host+'-'+prettylabel+"-"+ngpu+" GPUs"
                            if label in ['CUSOLVER']:
                                case=host+'-'+prettylabel
                            ax = df.plot(ax=ax, x='N', y='TCHOL', label=case)
        plt.grid()
        plt.xlabel('Matrix size')
        plt.ylabel('Time (s)')
        ax.legend(bbox_to_anchor=(1.1, 1.05))
        fig.savefig(plotname+'.pdf', bbox_inches='tight')
        fig.savefig(plotname+'.png', bbox_inches='tight')

#plot(['V100'], ['R-STATIC'], probndims, ['0','NA','exp'], ngpus, 'v100-exp')
#sys.exit(0)

#plot(hosts, ['MKL', 'R-STATIC', 'HiCMA-CPU'], probndims, kernels, ngpus, 'with-mkl')
#plot(hosts, ['R-STATIC','HiCMA-CPU',], probndims, kernels, ngpus, 'only-hicma')
#plot(hosts, ['R-STATIC'], probndims, kernels, ngpus, 'only-gpu')
#plot(['Vulture','A100'], ['R-STATIC','HiCMA-CPU','MAGMA','CUSOLVER'], probndims, ['0','NA','sqexp'], ngpus, 'a100-sqexp')
#plot(['Vulture','A100'], ['R-STATIC','HiCMA-CPU','MAGMA','CUSOLVER'], probndims, ['0','NA','exp'], ngpus, 'a100-exp')
plot(['Vulture','V100','IbexRome'], ['R-STATIC','HiCMA-CPU','MAGMA','CUSOLVER'], probndims, ['0','NA','sqexp'], ngpus, ['G','H'], 'v100-sqexp')
plot(['Vulture','V100','IbexRome'], ['R-STATIC','HiCMA-CPU','MAGMA','CUSOLVER'], probndims, ['0','NA','exp'], ngpus, ['G','H'], 'v100-exp')
#plot(['Vulture','A100','V100'], ['R-STATIC','HiCMA-CPU','MAGMA','CUSOLVER'], probndims, ['0','NA','exp'], ['0','NA','1'], 'a100-v100-exp')
#plot(['Vulture','A100','V100'], ['R-STATIC','HiCMA-CPU','MAGMA','CUSOLVER'], probndims, ['0','NA','sqexp'], ['0','NA','1'], 'a100-v100-sqexp')
sys.exit(0)
