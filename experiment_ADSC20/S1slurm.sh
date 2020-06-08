#!/bin/bash
#SBATCH --job-name=hicma
#SBATCH --account=k1205
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --partition=workq
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --hint=nomultithread
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=4
#SBATCH --time=10:00:00


#Please set this to available cores in the system
TH=31
P=2
export STARPU_SCHED=prio
export STARPU_SILENT=1

# you can set  max to number of mpi nodes * number of threads, and min to half the max
export STARPU_LIMIT_MAX_SUBMITTED_TASKS=128
export STARPU_LIMIT_MIN_SUBMITTED_TASKS=64


#10 viruses of 10370 resolution using 1e-5, 1e-6, and 1e-7
srun numactl --interleave=all ./build/timing/time_zpotrf_tile --m=103700 --n_range=2500:2500 --k=103700 --mb=2074 --nb=50 --nowarmup --threads=$TH --p=$P --rk=0 --acc=1e-5 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=100 --rbf_kernel=9 --rad=-1 --numobj=10 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/S1data/SortVirus103700.txt --csolve --solve 
srun numactl --interleave=all  ./build/timing/time_zpotrf_tile --m=103700 --n_range=2500:2500 --k=103700 --mb=2074 --nb=50 --nowarmup --threads=$TH --p=$P --rk=0 --acc=1e-6 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=100 --rbf_kernel=9 --rad=-1 --numobj=10 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/S1data/SortVirus103700.txt --csolve --solve 
srun numactl --interleave=all  ./build/timing/time_zpotrf_tile --m=103700 --n_range=2500:2500 --k=103700 --mb=2074 --nb=50 --nowarmup --threads=$TH --p=$P --rk=0 --acc=1e-7 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=100 --rbf_kernel=9 --rad=-1 --numobj=10 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/S1data/SortVirus103700.txt --csolve --solve 

#20 viruses of 10370 resolution using 1e-5, 1e-6, and 1e-7

srun numactl --interleave=all ./build/timing/time_zpotrf_tile --m=207400 --n_range=6800:6800 --k=207400 --mb=3050 --nb=100 --nowarmup --threads=$TH --p=$P --rk=0 --acc=1e-5 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=200 --rbf_kernel=9 --rad=-1 --numobj=20 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/S1data/SortVirus207400.txt --csolve --solve

srun numactl --interleave=all ./build/timing/time_zpotrf_tile --m=207400 --n_range=6800:6800 --k=207400 --mb=3050 --nb=100 --nowarmup --threads=$TH --p=$P --rk=0 --acc=1e-6 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=200 --rbf_kernel=9 --rad=-1 --numobj=20 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/S1data/SortVirus207400.txt --csolve --solve 

srun numactl --interleave=all ./build/timing/time_zpotrf_tile --m=207400 --n_range=6800:6800 --k=207400 --mb=3050 --nb=100 --nowarmup --threads=$TH --p=$P --rk=0 --acc=1e-7 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=200 --rbf_kernel=9 --rad=-1 --numobj=20 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/S1data/SortVirus207400.txt --csolve --solve

