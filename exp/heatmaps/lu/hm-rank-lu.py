import matplotlib as mpl
#print("matplotlib version:",mpl.__version__)
mpl.use('agg')
#['cairo', 'GTKCairo', 'gdk', 'MacOSX', 'GTK3Agg', 'Qt5Agg', 'GTKAgg', 'template', 'pdf', 'ps', 'WXAgg', 'GTK3Cairo', 'Qt4Agg', 'svg', 'WX', 'TkAgg', 'GTK', 'agg', 'WebAgg', 'pgf', 'nbAgg']
#%matplotlib inline

## divide one file into files with equal number of lines
#cat 2019-10-18-lu-1e2468-1.txt | grep -B 21  "initial_ranks" | grep -v "initial_ranks" |grep -v "\-\-"|split -l 21
#cat 2019-10-19-lu-acc-vul-1.sh | grep -B 25  "initial_ranks" | grep -v "initial_ranks" |grep -v "\-\-"|split -l 25 --additional-suffix=ACC

import os
import pandas as pd
import math
import numpy as np
import sys
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
#print(plt.get_backend())
from matplotlib import cm as CM
from pylab import  text
sys.path.append('./')
fs = 16
font = {'size'   : fs}
mpl.rc('font', **font)
import seaborn as sns;
keys=['rank_tile', 'nop_tile', 'nop_thread']
keys=['init_rank_tile','rank_tile']
filenames=[
    'out/xaa',
    #'out/xab',
    'out/xac',
    #'out/xad',
    'out/xae',
]
process_initial_ranks=True
filenames=[
'out/xaaACC',
'out/xabACC',
'out/xacACC',
'out/xadACC',
'out/xaeACC',
'out/xafACC',
        ]
filenames=[
'out/xaaACC147456',
'out/xabACC147456',
'out/xacACC147456',
'out/xadACC147456',
'out/xaeACC147456',
'out/xafACC147456',
        ]
# process_initial_ranks=False
# filenames=[
# 'out/hic-prio-147456-2304-800-1-1-40-_initialranks'
# ]
process_initial_ranks=True
#line 3217 of hicma-dev/exp/lu-outs/lu-sh4.log
filenames=[
'out/ntrian17598-init.txt',
'out/ntrian17598-final.txt'
]
def parse_ranks(filename, nrows, ncols, mat, display, title, vmin=None, vmax=None, colorbar=True, lttext=None):
    mat2=mat.astype('float')
    mat2[mat2 == 0] = np.NaN
    #std = np.nanstd(mat)
    mean = np.nanmean(mat2)
    minn  = np.nanmin(mat2)
    maxr = np.nanmax(mat2)

    info='Min:{:.0f}'.format(minn)
    #info='Avg:{:.0f}'.format(mean)
    #if display is not 'diff':
    info+='\nMax:{:.0f}'.format(maxr)

    fig, ax = plt.subplots()#figsize=(10,10))

    #text(0.15, 0.9,info,
    text(0.7, 0.9,info,
     horizontalalignment='left',
     verticalalignment='center',
     transform = ax.transAxes,
     color='black',backgroundcolor='white')

    # print(info)
    if lttext is not None:
        text(0.7, 0.4, lttext,
     horizontalalignment='left',
     verticalalignment='center',
     transform = ax.transAxes,
     color='black',backgroundcolor='white')
        # print(lttext)

    if vmax is None:
        vmax = np.max(mat)
    else:
        np.fill_diagonal(mat, vmax)
    cmap = CM.get_cmap('jet') # jet doesn't have white color
    cmap.set_bad('w')
    im=ax.imshow(mat,cmap=cmap, vmin=vmin, vmax=vmax)
    # # inset axes https://matplotlib.org/gallery/subplots_axes_and_figures/zoom_inset_axes.html
    # axins = ax.inset_axes([0.52, 0.55, 0.45, 0.45])
    # axins.imshow(mat, cmap=cmap, vmin=0, vmax=vmax)
    # # sub region of the original image
    # x1, x2, y1, y2 = 190,200,200,190
    # axins.set_xlim(x1, x2)
    # axins.set_ylim(y1, y2)
    # axins.set_xticklabels('')
    # axins.set_yticklabels('')
    # ax.indicate_inset_zoom(axins)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    # print(filename)
    # print(mat)
    if colorbar is True:
        fig.colorbar(im, fraction=0.046, pad=0.04)
    fig.tight_layout()
    if title is not None: plt.title(title)
    plotname=filename+'_'+display+'.pdf'
    print("Plot is generated:",plotname)
    fig.savefig(plotname, bbox_inches = 'tight');
    #plt.show()

def process_file(fname, counts, key, colorbarmax, process_initial_ranks=False, draw_colorbar=False):
    with open(fname) as f:
        lines = f.readlines()
    start=0
    if True:
        if True:
            line = lines[start]
            nrows, ncols = line.split()
            nrows=int(nrows)
            ncols=int(ncols)
            #print(key, st, nrows, ncols)
            start = start + 1
            for i in range(0,nrows):
                line=lines[start+i]
                #print(i,line[0:10],line[-10:])
                row = [int(x) for x in line.split()]
                counts[key].append(row)
            counts[key]=np.array(counts[key])
            # print(fname)
            # print(counts[key])

            if False and (key is 'init_rank_tile' or key is 'rank_tile' or key is 'nop_tile'): # or key is 'op_tile':
                mask=np.tri(counts[key].shape[0],k=0).T
                #print(mask)
                counts[key] = np.ma.array(counts[key], mask=mask)
            #if key is 'rank_tile':
            #    np.fill_diagonal(counts[key], 0)
            mat=counts[key].astype('float')
            mat[mat == 0] = np.NaN
            std = np.nanstd(mat)
            mean = np.nanmean(mat)
            maxr=np.nanmax(mat)
            cov = std/mean
            if std  > 10000: std  /= 1e9
            if mean > 10000: mean /= 1e9
            title = key \
            +'\n mean:'+'{:.2f}'.format(mean) \
            +'\n std:'+'{:.2f}'.format(std)\
            +'\n cov:'+'{:.2f}'.format(cov)\
            +'\n maxr:'+'{:.2f}'.format(maxr)
            # print(title)
            title=None
            if colorbarmax < maxr:
                colorbarmax = maxr
            if process_initial_ranks is True or key is 'rank_tile':
                parse_ranks(fname+'INIT', nrows, ncols, counts['init_rank_tile'], display='init_rank_tile', title=title, vmin=0, vmax=colorbarmax, colorbar=draw_colorbar)
                if process_initial_ranks is False:
                    parse_ranks(fname+'FINAL', nrows, ncols, counts['rank_tile'], display='rank_tile', title=title, vmin=0, vmax=colorbarmax, colorbar=True)
                    mat=counts['rank_tile']
                    mat2=mat.astype('float')
                    mat2[mat2 == 0] = np.NaN
                    mean = np.nanmean(mat2)
                    maxr = np.nanmax(mat2)
                    lttext='Final\navg:{:.0f}'.format(mean)
                    lttext+='\nFinal\nmax:{:.0f}'.format(maxr)
                    parse_ranks(fname+'DIFF', nrows, ncols, counts['rank_tile']-counts['init_rank_tile'], display='diff', title=title, lttext=lttext, vmax=colorbarmax,colorbar=True)
    return colorbarmax

colorbarmax=0
for j in filenames:
    counts={}
    counts['init_rank_tile']=[]
    counts['rank_tile']=[]
    counts['nop_tile']=[]
    counts['nop_thread']=[]
    print("Processing:", j)
    if process_initial_ranks is False:
        colorbarmax=0

    #colorbarmax=2304
    colorbarmax=1200
    colorbarmax=process_file(j, counts, 'init_rank_tile',colorbarmax, process_initial_ranks=process_initial_ranks)
    #colorbarmax=2304
    colorbarmax=1200

    finalrankfile = j+'_finalranks'
    if os.path.isfile(finalrankfile) is True:
        process_file(finalrankfile, counts, 'rank_tile', colorbarmax)
if process_initial_ranks is True:
    print("redrawing with colorbarmax:", colorbarmax)
    draw_colorbar=False
    for ij, j in enumerate(filenames):
        counts={}
        counts['init_rank_tile']=[]
        counts['rank_tile']=[]
        counts['nop_tile']=[]
        counts['nop_thread']=[]
        print("Processing:", j)
        if ij == len(filenames)-1:
            print("Drawing colorbar")
            draw_colorbar=True
        process_file(j, counts, 'init_rank_tile', colorbarmax, process_initial_ranks=process_initial_ranks, draw_colorbar=draw_colorbar)
    #break
