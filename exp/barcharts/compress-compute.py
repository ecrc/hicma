import pandas as pd
import math
import numpy
import sys
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib as mpl
sys.path.append('./')
fs = 16
font = {'size'   : fs}
mpl.rc('font', **font)


scal_cols=['code','kernel','m','t1','nb','nmpi','p','q','time','gflops','t2']


hsched_cols2=['sched','m','mb','nnodes','p','q','time','nb','maxsub','minsub','nthreads', 'jobid','iavgrk', 'iminrk', 'imaxrk', 'favgrk', 'fminrk', 'fmaxrk', 'tproblem','tcompress']
colspecv2=['sched','m','mb','nnodes','p','q','time','nb','maxsub','minsub','nthreads', 'jobid','iavgrk', 'iminrk', 'imaxrk', 'favgrk', 'fminrk', 'fmaxrk', 'tproblem','tcompress', 'fixedrank', 'fixedacc','wavek', 'hostname','binary','shmaxrk', 'shprob', 'shdecay']
colspecv3=['code','sched','m','mb','nnodes','p','q','time','gflop/s','nb','maxsub','minsub','nthreads', 'jobid','iavgrk', 'iminrk', 'imaxrk', 'favgrk', 'fminrk', 'fmaxrk', 'tproblem','tcompress', 'fixedrank', 'fixedacc','wavek', 'hostname','binary','shmaxrk', 'shprob', 'shdecay']
fhminmax='res/2017-09-24-hicma-maxrk50-1.txt'
fhminmax='res/2017-10-17-hicma-2M-1.txt'
fhminmax='res/2017-10-17-hicma-2M-2.txt'
fhminmax='res/2017-10-17-hicma-2M-3.txt'
fhicma='res/2017-10-18-hicma-3M-1.txt'




global_nplist=[8, 16, 32, 64,128, 256, 512, 1024]
#df_hicma = df_hicma[df_hicma.m <= 594000]
#df_hicma = df_hicma[df_hicma.m > 594000]
#df_scalm = df_scalm[df_scalm.m > 594000]
#print(df_hicma)
from common import *
#print(katrial)
#print(appmarkers)
#sys.exit()


def fig_finalize(fig, ax, plt, plotname, legncol=2, legloc=4, ylim=None,
                 xpoints=[80,170,260,350,440,600,800],
                 ypoints=[200, 330, 400, 500,  800, 1000, 1200, 1800, 2500, 3000, 8000]
                ):
    print(plotname)
    #ax.set_xlabel("Matrix size",fontweight='bold')
    #ax.set_ylabel("Time(s)",fontweight='bold')

    #proposal
    ax.set_xlabel("Matrix size")
    ax.set_ylabel("Time(s)")
    if ylim is not None: ax.set_ylim(ylim)

    ax.yaxis.grid(True, which='both')
        #ax.set_title('TLR (HiCMA) vs Dense (ScaLAPACK-libsci) POTRF Results on Shaheen')

    if legncol is None:
        legncol=1
    if legloc is None:
        legloc=4
    if legncol == 100:
        plt.legend(handlelength=3,loc=4, ncol=3, columnspacing=1, bbox_to_anchor=(1.1, 1.05))
            #lgd=plt.legend(handlelength=3,loc='lower right',fontsize=12) #, bbox_to_anchor=(1, 0.5))
            #bbox_extra_artists=(lgd,), bbox_inches='tight'
    else:
        plt.legend(handlelength=3,loc=legloc, ncol=legncol, columnspacing=1) #fontsize


    ax.set_yscale('log')
    ax.set_xscale('log')

    if plotname is 'time_hsw_cluster' or plotname is 'time_hsw_cluster_with_tcompress':
        plt.minorticks_off()
        points=[4,6,8,10,14,20,26,34,44]
        ticks=[i*13500 for i in points ]
        print('xxticks', ticks)
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
        fig.canvas.draw()
        labels = [str(int(numpy.rint(i*13500/1000)))+'K' for i in points]
        print("labels:",labels)
        ax.set_xticklabels(labels)
        ax.xaxis.grid(True, which='major')
        ax.yaxis.grid(True, which='major')
        ax.yaxis.grid(False, which='minor')
        ax.grid(color='gray', linestyle='-', linewidth=.5)
        for child in ax.get_children():
            if isinstance(child, mpl.spines.Spine):
                child.set_color('gray')
    elif plotname.startswith('time_hsw_cluster_BIG') or  plotname.startswith('time_skl_cluster_BIG'):
        plt.minorticks_off()
        ax.xaxis.grid(True, which='major')
        ax.yaxis.grid(True, which='major')


        ticks=[i*13500 for i in xpoints ]
        print('xticks', ticks)
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
        #fig.canvas.draw()
        labels = [str(int(numpy.rint(i*13500/1000000)))+'M' for i in xpoints]
        print("labels:",labels)
        ax.set_xticklabels(labels)

        ticks=[i for i in ypoints]
        print('yticks', ticks)
        ax.set_yticks(ticks)
        #fig.canvas.draw()
        labels = ["{0:.0f}".format(i/(60)) for i in ypoints ]
        print("labels:",labels)
        ax.set_yticklabels(labels)
#        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
        ax.set_ylabel("Time(minutes)") #,fontweight='bold')
        ax.grid(color='gray', linestyle='-', linewidth=.5)
        for child in ax.get_children():
            if isinstance(child, mpl.spines.Spine):
                child.set_color('gray')
    else:
        plt.minorticks_on()#ax.minorticks_on()
        #ax.xaxis.grid(True, which='both')
        #return
        plt.minorticks_off()
        points=[4,8,12,16,24,32,40]
        if plotname == 'time_TRSM_hsw_cluster':
            points=[4,8,12,16,24]
        ticks=[i*13500 for i in points ]
        print('ticks', ticks)
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
        fig.canvas.draw()
        labels = [str(int(numpy.rint(i*13500/1000)))+'K' for i in points]
        print("labels:",labels)
        ax.set_xticklabels(labels)
        ax.xaxis.grid(True, which='major')
        ax.yaxis.grid(True, which='major')
        ax.yaxis.grid(False, which='minor')
        ax.grid(color='gray', linestyle='-', linewidth=.5)
        for child in ax.get_children():
            if isinstance(child, mpl.spines.Spine):
                child.set_color('gray')
    fig.savefig(plotname+'_ylog.pdf', bbox_inches = 'tight');
    plt.show()

    #ax.set_yscale('linear')
    #fig.savefig(plotname+'_ylinear.pdf');
    #plt.show()


def plot(df_hicma, df_scalm, plotname, show_data_labels=True, show_data=False, yhicma_display_compress=False,
         nplist=global_nplist,
         scalnplist=global_nplist,
         df_hicmaB=None, legncol=None, legloc=None, ylim=None, applist=None, yval='time',
         xpoints=[80,170,260,350,440,600,800],
         ypoints=[200, 330, 400, 500,  800, 1000, 1200, 1800, 2500, 3000, 8000],
         colnrows='m',colnumnodes='nnodes',showdatacols=['m','mb','nb','shmaxrk','fmaxrk','favgrk','time','jobid','tcompress']
        ):
    #if plotname is not 'time_hsw_cluster':
    #    return
    ax=None
    fig=None
    figsize_width=8
    figsize_height=8

    #proposal
    #figsize_width=9
    #figsize_height=5.8

    if show_data_labels is True:
        figsize_width=15
        figsize_height=100

    fig, ax = plt.subplots(figsize=(figsize_width, figsize_height))
    npcolors={8:'red',16:'green',32:'magenta',64:'blue',128:'pink',256:'cyan',512:'black'}
    goodnpcolors={8:'red',16:'#0ca7ab',32:'#59529c',64:'#fa3c3c',128:'#ecaf4d',256:'#6ea35e',512:'#26df2b'#406c9d'
                  ,1024:'blue'}
    """
    HICMA_STARSH_PROB_RND    1
 HICMA_STARSH_PROB_SS     2
 HICMA_STARSH_PROB_RNDUSR 3
 HICMA_STARSH_PROB_FILE   4
 HICMA_STARSH_PROB_GEOSTAT   5
 HICMA_STARSH_PROB_EDSIN  6
 """
    if True:
        if True:
            if True:
                for inp, np in enumerate(scalnplist):
                    label='ScaLAPACK '+str(np) +' nodes'
                    #nmpi=np*nmpi_per_node
                    df = df_scalm[df_scalm.nnodes==np]
                    if len(df) == 0:
                        continue
                    ax=df.plot(ax=ax,x='m', y='time',ls='-',marker='o',label=label,
                               color=goodnpcolors[np],
                              markersize=markersize,mfc='white', mew=markeredgeweight,
                               lw=linewidth,
                               zorder=100, clip_on=False)
                    if show_data is True:
                        print('ScaLAPACK - '+str(np))
                        print(df[['m','nnodes','p','q','time']])

            if True:
                for inp, np in enumerate(nplist):
                    singleapp = False
                    if len(applist) == 1:
                        singleapp = True
                    for app in applist:
                        if singleapp:
                            label='HiCMA-'+str(np)# +' nodes'
                            linestyle='-'
                            marker='o'
                            color=goodnpcolors[np]
                        else:
                            label='HiCMA-'+appnames[app]+'-'+str(np)# +' nodes'
                            linestyle='--'
                            marker=appmarkers[app]
                            color='black'
                            #for label,color,marker in [[str(minsub)+' '+sched+' - HiCMA - '+str(np)+' nodes',npcolors[inp],"x"]]:
                        #df=df_hsched[(df_hsched['numnodes']==np)&(df_hsched['sched']==sched)&(df_hsched['iminsub']==minsub)]
                        dfnp=df_hicma[(df_hicma[colnumnodes]==np)&(df_hicma['sched']=='prio')&(df_hicma['shprob']==app)]
                        df=dfnp.loc[dfnp.groupby([colnrows])["time"].idxmin()]
                        if len(df) == 0:
                            continue
                        if show_data is True:
                            print(label)
                            #print(dfnp[['nnodes','m','mb','nb','minsub','time','jobid','tcompress']])
                            print(df[showdatacols]) #,'wavek'
                        ax=df.plot(ax=ax,x=colnrows, y=yval,ls=linestyle,marker=marker,label=label,
                                   color=color,
                                  markersize=markersize,mfc='white', mew=markeredgeweight,
                                   lw=linewidth, zorder=100, clip_on=False)
                        if show_data_labels is True:
                            for xy,m,nb,jobid,time in zip(zip(df.m, df.time),df.m,df.mb,df.jobid,df.time):                                       # <--
                                ax.annotate(str(int(m/1000))+' '+str(nb)+'\n'+str(int(time))+'\n'+str(jobid), xy=xy, textcoords='data') # <--
                        if yhicma_display_compress is True:
                            ax=df.plot(ax=ax,x='m', y='tcompress',ls=None,marker='X',label='$T_{compress}$ - '+str(np),
                                       color=npcolors[inp],markersize=markersize)
            if df_hicmaB is not None and colnumnodes != 'PQ':
                for inp, np in enumerate(nplist):
                    label='HiCMA KNL '+str(np)
                #for label,color,marker in [[str(minsub)+' '+sched+' - HiCMA - '+str(np)+' nodes',npcolors[inp],"x"]]:
                    #df=df_hsched[(df_hsched['numnodes']==np)&(df_hsched['sched']==sched)&(df_hsched['iminsub']==minsub)]
                    dfnp=df_hicmaB[(df_hicmaB[colnumnodes]==np)&(df_hicmaB['sched']=='prio')]
                    df=dfnp.loc[dfnp.groupby([colnrows])["time"].idxmin()]
                    if len(df) == 0:
                        continue
                    ax=df.plot(ax=ax,x=colnrows, y='time',ls=':',marker='s', mew=markeredgeweight, lw=linewidth,
                               label=label,color=npcolors[np],markersize=markersize)
                    #print(df[['m','nb','time']])
    fig_finalize(fig, ax, plt, plotname, legncol, legloc, ylim, xpoints, ypoints)

def process(plotname, hicma_files, scal_files,  show_data_labels=False, show_data=False,
            mlowerbound=0, mupperbound=1000000000, yhicma_display_compress=False,
            nplist=global_nplist,
            scalnplist=global_nplist,
            hicma_filesB=None, legncol=None, legloc=None, ylim=None,colspecs=None,applist=None,filterfunc=None,
            yval='time',waveklist=None,
            xpoints=[80,170,260,350,440,600,800],
            ypoints=[200, 330, 400, 500,  800, 1000, 1200, 1800, 2500, 3000, 8000],
            colnrows='m',colnumnodes='nnodes',showdatacols=['m','mb','nb','shmaxrk','fmaxrk','favgrk','time','jobid','tcompress']
           ):
    frames = []
    if colspecs is None:
        print("define colspecs")
        sys.exit()
    if len(colspecs) != len(hicma_files):
        print("colspecs and hicma_files should have same length")
        sys.exit()
    for ifile, f in enumerate(hicma_files):
        #print(f,colspecs[ifile])
        df = pd.read_csv(f, delim_whitespace=True,engine='python',index_col=False,header=None,names=colspecs[ifile])
        #print(df)
        frames.append(df)
    df_hicma = pd.concat(frames, ignore_index=True)
    #print(df_hicma.columns)
    #return
    #sys.exit()
    if 'shprob' in df_hicma.columns:
        for index, row in df_hicma.iterrows():
            if math.isnan(row['shprob']):
                df_hicma.set_value(index, 'shprob', 2)
    else:
        df_hicma['shprob']=2
    if 'shmaxrk' not in df_hicma.columns:
        df_hicma['shmaxrk']=math.nan
    if 'sched' not in df_hicma.columns:
        df_hicma['sched'] = 'prio'
    if 'time' not in df_hicma.columns:
        df_hicma['time'] =  df_hicma['tpotrf']+df_hicma['ttrsm1']+df_hicma['ttrsm2']
    #print(df_hicma)
    #sys.exit()
    framesB = []
    if hicma_filesB is not None:
        for f in hicma_filesB:
            df = pd.read_csv(f, delim_whitespace=True,engine='python',index_col=False,header=None,names=hsched_cols2)
            framesB.append(df)
        df_hicmaB = pd.concat(framesB, ignore_index=True)
        df_hicmaB = df_hicmaB[df_hicmaB[colnrows] > mlowerbound]
        df_hicmaB = df_hicmaB[df_hicmaB[colnrows] < mupperbound]
    else:
        df_hicmaB = pd.DataFrame(columns=hsched_cols2)

    if len(scal_files) > 0:
        df_scal = pd.read_csv(scal_files[0], delim_whitespace=True,engine='python',index_col=False,header=None,names=scal_cols)
        df_scalm=df_scal.loc[df_scal.groupby(['m','nmpi'])["nb"].idxmin()]
        df_scalm['nnodes']=df_scalm.nmpi/nmpi_per_node
        df_scalm['nnodes']=df_scalm['nnodes'].astype(int)
        df_scalm = df_scalm[df_scalm['m'] > mlowerbound]
        df_scalm = df_scalm[df_scalm['m'] < mupperbound]
    else:
        df_scalm = pd.DataFrame(columns=scal_cols)
        df_scalm['nnodes']=0
    #print(df_scal)
    #sys.exit()

    if yval is "tcompcalc":
        df_hicma['tcompcalc'] = df_hicma['time'] + df_hicma['tcompress']

    df_hicma = df_hicma[df_hicma[colnrows] > mlowerbound]

    df_hicma = df_hicma[df_hicma[colnrows] < mupperbound]

    if filterfunc is not None:
        df_hicma=filterfunc(df_hicma)

    if waveklist is not None:
        df_hicma = df_hicma[(df_hicma['wavek'].isin(waveklist))]
    df_hicma = df_hicma[df_hicma[colnumnodes] < 1024]
    #dftmp=df_hicma#[df_hicma['shprob']==6.0]
    #print(dftmp[['shprob','nnodes','m','mb','minsub','time','jobid','tcompress']])
    #sys.exit()
    plot(df_hicma, df_scalm, plotname, show_data_labels, show_data, yhicma_display_compress,
         nplist, scalnplist, df_hicmaB, legncol, legloc, ylim, applist, yval, xpoints, ypoints,colnrows,colnumnodes,
        showdatacols)
    return df_hicma, df_scalm

fscal='res/2017-2017-10-20-KNL-scal-small-1-grepped.txt'
nmpi_per_node=64
files=['res/2017-10-21-hicma-compare-KNL-1.txt']
if False: process(plotname='time_compare_hsw_cluster_MED', hicma_files=files,
        scal_files=[], hicma_filesB=['res/2017-10-20-KNL-hicma-small-1-grepped.txt'], show_data_labels=False,
        show_data=False, yhicma_display_compress=False, mlowerbound=539999, mupperbound=1100000, legncol=100)

files=['res/2017-10-20-KNL-hicma-small-1-grepped.txt']
if False:process(plotname='time_knl_cluster', hicma_files=files, scal_files=[fscal],
        show_data_labels=False, show_data=False, mupperbound=594000, nplist=[16],legloc=2,ylim=[4,10000])
if False:process(plotname='time_knl_cluster_BIG', hicma_files=files, scal_files=[],
        show_data_labels=False, show_data=False, mlowerbound=539999, mupperbound=1100000)


files=['res/2017-10-19-hicma-50K500K-1.txt', 'res/2017-10-18-hicma-3M-1.txt', 'res/2017-10-20-hicma-10M-1.txt']
files=['res/2018-02-11-hicma-16-1.txt'  ,'res/2017-10-22-hicma-1.txt'      ]
files=['res/2018-02-11-hicma-16-1.txt']
colspecs=[colspecv2         ]

fscal='res/2017-06-02-scal-2.txt'
nmpi_per_node=32

shprob=2
#shprob=6
fscal2='res/2018-02-23-scal-1.txt'
df_scal = pd.read_csv(fscal2, delim_whitespace=True,engine='python',index_col=False,header=None,names=scal_cols)
df_scalm=df_scal.loc[df_scal.groupby(['m','nmpi'])["time"].idxmin()]
df_scalm['nnodes']=df_scalm.nmpi/nmpi_per_node
df_scalm['nnodes']=df_scalm['nnodes'].astype(int)
df_scalm=df_scalm[df_scalm['nnodes']==16]
df_scalm['mK']=df_scalm['m']/1000
df_scalm['mK']=df_scalm['mK'].astype(int)

fdense='res/2018-02-17-only-dense-1.txt'
df_dense=pd.read_csv(fdense, delim_whitespace=True,engine='python',index_col=False,header=None,names=colspecv2)
df_dense=df_dense[df_dense['shprob']==shprob]
df_dense=df_dense.loc[df_dense.groupby(['m','nnodes'])["tcompress"].idxmin()]
df_dense=df_dense[df_dense['nnodes']==16]


fhicma='res/2018-02-11-hicma-16-1.txt'
#df_hicma=pd.read_csv(fhicma, delim_whitespace=True,engine='python',index_col=False,header=None,names=colspecv2)
fhicma='/Users/kadirakbudak/hicma-dev/exp/out/vulture-hicma-lu-2019-10-26-custommaxrank-1-parsed.txt' #nipp=3
fhicma='/Users/kadirakbudak/hicma-dev/exp/out/vulture-hicma-lu-2019-10-28-custommaxrank-nipp6-O3-1-parsed.txt'
df_hicma=pd.read_csv(fhicma, delim_whitespace=True,engine='python',index_col=False,header=None,names=colspecv3)

df_cham=df_hicma[df_hicma['code']=='cham']
df_cham['tcompress']=df_cham['tproblem']
df_cham['tcompcalc'] = df_cham['time'] + df_cham['tcompress']
df_cham['mK']=df_cham['m']/1000
df_cham['mK']=df_cham['mK'].astype(int)

df_hicma=df_hicma[df_hicma['code']=='hic']
df_hicma=df_hicma[df_hicma['shprob']==shprob]
df_hicma=df_hicma.loc[df_hicma.groupby(['m','nnodes'])["time"].idxmin()]
df_hicma['tcompcalc'] = df_hicma['time'] + df_hicma['tcompress']
df_hicma=df_hicma[df_hicma['m']<=486000]
df_hicma['mK']=df_hicma['m']/1000
df_hicma['mK']=df_hicma['mK'].astype(int)

#print(df_hicma)
#print(df_dense)
print(df_dense[['m','tcompress']])
for index, row in df_scalm.iterrows():
    m=row['m']
    time=row['time']
    print(m)
    try:
        tdense=df_dense[df_dense['m']==m]['tcompress'].tolist()[0]
    except :
        print("Error occured", file=sys.stderr )
        continue
    df_scalm.set_value(index, 'tcompcalc', tdense+time)
    df_scalm.set_value(index, 'tcompress', tdense)


#scal_cols=['code','kernel','m','t1','nb','nmpi','p','q','time','gflops','t2']
#colspecv2=['sched','m','mb','nnodes','p','q','time','nb','maxsub','minsub',
#'nthreads', 'jobid','iavgrk', 'iminrk', 'imaxrk', 'favgrk', 'fminrk', 'fmaxrk',
#'tproblem','tcompress', 'fixedrank', 'fixedacc','wavek', 'hostname','binary','shmaxrk', 'shprob', 'shdecay']


#print(df_scalm[['nnodes','m','time','tcompcalc']])
#print(df_dense[['nnodes','m','tcompress']])
print('Chameleon')
print(df_cham[['nnodes','m','tcompress','time','tcompcalc']])
print('HiCMA')
print(df_hicma[['nnodes','m','tcompress','time','tcompcalc']])



figsize_width=6
figsize_height=8

fig, ax = plt.subplots(figsize=(figsize_width, figsize_height))
axes = plt.gca()
#axes.set_xlim([xmin,xmax])
axes.set_ylim([0,1600])
df_hicma[['tcompress','time','mK']].plot(kind='bar', stacked=True,x='mK',ax=ax,width=.9,zorder=2)
#df_scalm[['tcompress','tcompcalc','m']].plot(kind='bar', stacked=True,x='m',ax=ax)
ax.set_xlabel("Matrix size (x1000)")
ax.set_ylabel("Time(s)")
ax.legend(['Generation&compression','HiCMA-LU'], loc=2);
ax.grid(zorder=1)
fig.savefig('hicma-lu.pdf', bbox_inches = 'tight');
#plt.show()
#sys.exit()

fig, ax = plt.subplots(figsize=(figsize_width, figsize_height))
df_cham[['tcompress','time','mK']].plot(kind='bar', stacked=True,x='mK',ax=ax,width=.9, zorder=2)
ax.set_xlabel("Matrix size (x1000)")
ax.set_ylabel("Time(s)")
ax.legend(['Dense matrix generation','Chameleon-LU'], loc=2);
ax.grid(zorder=1)
fig.savefig('chameleon-lu.pdf', bbox_inches = 'tight');
#plt.show()


