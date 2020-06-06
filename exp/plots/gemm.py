#!/usr/bin/env python3

#/Users/akbudak/anaconda3/bin/python

import re
import sys
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
if len(sys.argv) < 2:
    print("Usage: gemm.py outputFile")
    sys.exit(-1);

print("This is the name of the script: ", sys.argv[0])
print("Number of arguments: ", len(sys.argv))
print("The arguments are: " , str(sys.argv))

rawfile=sys.argv[1]
print("Result file:", rawfile)
with open(rawfile) as f:
    lines = f.readlines()
allres=[]
for (iline,line) in enumerate(lines):
    if line.startswith("# MB:"):
        res={}
        res["mb"]=int(line.split(':')[1])
    if line.startswith("# fixed acc:"):
        res["acc"]=float(line.split(':')[1])
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
        res["gflop/s"] = float(arr[3])
        res["time"]    = float(arr[4])
        line=lines[iline+1]
        arr=line.strip().split(' ');
        #print(arr)
        res["M"]     =  int(arr[0])
        #res["N"]     = int(arr[1])
        #res["K"]     = int(arr[2])
        #res["time2"] = float(arr[3])
        print(res)
        allres.append(res)


df=pd.DataFrame.from_dict(allres)
print(df)

sys.exit(0)
df2 = df.loc[df.groupby(['M','rip'])['time'].idxmin()]
print(df2)
df3=df2[['M','rip','gflop','time']]
print(df3)
#plt.subplots()
#for val in ['time','gflop']:

val='gflop'
what='flops'

for (val, what, ylabel) in [
        ('gflop','flops','GFlop')
        , ('time','time','Time(s)')
        #, ('favg','final rank','Final Average Rank')
        ]:
    table = pd.pivot_table(df3, values=[val], index=['M'], columns=['rip'])
    print("Table: ",table)
    table.plot(kind='bar')

    table.reset_index(inplace=True)


    ax=plt.gca()
    percent_of_goal = ["{}%".format(int(100.*row[1]/row[2]-100)) for name,row in table.iterrows()]
    pairs = len(df3['M'].unique())
    make_pairs = zip(*[ax.get_children()[:pairs],ax.get_children()[pairs:pairs*2]])
    height_shift = ax.get_ylim()[1]*.05
    for i,(left, right) in enumerate(make_pairs):
        if  percent_of_goal[i].startswith("-") is False:
            height=min(left.get_bbox().y1,right.get_bbox().y1)
            t=ax.text(i+.05,right.get_bbox().y1+height_shift,"   "+percent_of_goal[i]+"\nless "+what, horizontalalignment ='left')
            t.set_bbox(dict(facecolor='red', alpha=0.5, edgecolor='red'))
        t=ax.text(i,left.get_bbox().y1-height_shift, "{:.0f}".format(table.iloc[i][1]), horizontalalignment ='right')
        t=ax.text(i,right.get_bbox().y1-height_shift, "{:.0f}".format(table.iloc[i][2]), horizontalalignment ='left')

    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
    ax.patch.set_facecolor('#FFFFFF')
    #ax.spines['bottom'].set_color('#CCCCCC')
    ax.spines['bottom'].set_linewidth(0)
    ax.spines['left'].set_linewidth(1)
    ax.spines['top'].set_linewidth(0)
    ax.spines['right'].set_linewidth(0)

    plt.xlabel('Matrix size')
    plt.ylabel(ylabel)
    plt.legend( ['No reordering','Increasing order of ranks'])
    plt.grid()
    plt.savefig(val+'.pdf', bbox_inches = 'tight')
    plt.savefig(val+'.png', bbox_inches = 'tight')
