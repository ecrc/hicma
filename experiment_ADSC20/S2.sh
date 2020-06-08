#!/bin/bash
#Please set this to available cores in the system
TH=39
export STARPU_SCHED=prio
export STARPU_SILENT=1

#8 viruses of 10370 resolution in which every 2 viurses interact together in one batch using 1e-6
numactl --interleave=all ./build/timing/time_zpotrf_tile_batch --m=20740 --n_range=1000:1000 --k=20740 --mb=1037 --nb=50 --nowarmup --threads=$TH --p=1 --rk=0 --acc=1e-6 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=100 --rbf_kernel=9 --rad=-1 --numobj=8  --numsubobj=2 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/S2data_20k/ --csolve --solve 

#240 viruses of 10370 resolution in which every 20 viurses interact together in one batch using 1e-6

numactl --interleave=all ./build/timing/time_zpotrf_tile_batch --m=207400 --n_range=3400:3400 --k=207400 --mb=3050 --nb=50 --nowarmup --threads=$TH --p=1 --rk=0 --acc=1e-6 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=100 --rbf_kernel=9 --rad=-1 --numobj=240 --numsubobj=20 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/S2data_200k/ 



