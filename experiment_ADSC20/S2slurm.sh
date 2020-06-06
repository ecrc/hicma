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


#240 viruses of 10370 resolution in which every 20 viurses interact together in one batch using 1e-6

srun numactl --interleave=all ./build/timing/time_zpotrf_tile_batch --m=207400 --n_range=3400:3400 --k=207400 --mb=3050 --nb=50 --nowarmup --threads=$TH --p=$P --rk=0 --acc=1e-6 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=100 --rbf_kernel=9 --rad=-1 --numobj=240 --numsubobj=20 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/S2data_200k/ 

