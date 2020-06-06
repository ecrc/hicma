#!/bin/bash
#Please set this to available cores in the system
TH=39
export STARPU_SCHED=prio
export STARPU_SILENT=1

#10 viruses of 10370 resolution using 1e-5, 1e-6, and 1e-7
numactl --interleave=all ./build/timing/time_zpotrf_tile --m=103700 --n_range=2500:2500 --k=103700 --mb=2074 --nb=50 --nowarmup --threads=$TH --p=1 --rk=0 --acc=1e-5 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=100 --rbf_kernel=9 --rad=-1 --numobj=10 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/S1data/SortVirus103700.txt --csolve --solve 
numactl --interleave=all  ./build/timing/time_zpotrf_tile --m=103700 --n_range=2500:2500 --k=103700 --mb=2074 --nb=50 --nowarmup --threads=$TH --p=1 --rk=0 --acc=1e-6 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=100 --rbf_kernel=9 --rad=-1 --numobj=10 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/S1data/SortVirus103700.txt --csolve --solve 
numactl --interleave=all  ./build/timing/time_zpotrf_tile --m=103700 --n_range=2500:2500 --k=103700 --mb=2074 --nb=50 --nowarmup --threads=$TH --p=1 --rk=0 --acc=1e-7 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=100 --rbf_kernel=9 --rad=-1 --numobj=10 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/S1data/SortVirus103700.txt --csolve --solve 

#20 viruses of 10370 resolution using 1e-5, 1e-6, and 1e-7

numactl --interleave=all ./build/timing/time_zpotrf_tile --m=207400 --n_range=6800:6800 --k=207400 --mb=3050 --nb=100 --nowarmup --threads=$TH --p=1 --rk=0 --acc=1e-5 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=200 --rbf_kernel=9 --rad=-1 --numobj=20 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/S1data/SortVirus207400.txt --csolve --solve

numactl --interleave=all ./build/timing/time_zpotrf_tile --m=207400 --n_range=6800:6800 --k=207400 --mb=3050 --nb=100 --nowarmup --threads=$TH --p=1 --rk=0 --acc=1e-6 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=200 --rbf_kernel=9 --rad=-1 --numobj=20 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/S1data/SortVirus207400.txt --csolve --solve 

numactl --interleave=all ./build/timing/time_zpotrf_tile --m=207400 --n_range=6800:6800 --k=207400 --mb=3050 --nb=100 --nowarmup --threads=$TH --p=1 --rk=0 --acc=1e-7 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=200 --rbf_kernel=9 --rad=-1 --numobj=20 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/S1data/SortVirus207400.txt --csolve --solve

