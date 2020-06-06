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

#virus of 10370 resolution using 1e-5, 1e-6, and 1e-7
srun numactl --interleave=all ./build/timing/time_zpotrf_tile --m=10370 --n_range=500:500 --k=10370 --mb=1037 --nb=50 --nowarmup --threads=$TH --p=$P --rk=0 --acc=1e-5 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=100  --rbf_kernel=9 --rad=0.6 --numobj=1 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus10370.txt --csolve --solve

srun numactl --interleave=all ./build/timing/time_zpotrf_tile --m=10370 --n_range=600:600 --k=10370 --mb=1037 --nb=60 --nowarmup --threads=$TH  --p=$P --rk=0 --acc=1e-6 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=120  --rbf_kernel=9 --rad=0.6 --numobj=1 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus10370.txt --csolve --solve

srun numactl --interleave=all ./build/timing/time_zpotrf_tile --m=10370 --n_range=1000:1000 --k=10370 --mb=1037 --nb=100 --nowarmup --threads=$TH  --p=$P --rk=0 --acc=1e-7 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=200  --rbf_kernel=9 --rad=0.6 --numobj=1 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus10370.txt --csolve --solve

#virus of 44932 resolution using 1e-5, 1e-6, and 1e-7

srun numactl --interleave=all ./build/timing/time_zpotrf_tile --m=44932 --n_range=2250:2250 --k=44932 --mb=1000 --nb=50 --nowarmup --threads=$TH  --p=$P --rk=0 --acc=1e-5 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=100 --rbf_kernel=9 --rad=0.6 --numobj=1 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus44932.txt --csolve --solve 

srun numactl --interleave=all ./build/timing/time_zpotrf_tile --m=44932 --n_range=2700:2700 --k=44932 --mb=1000 --nb=60 --nowarmup --threads=$TH  --p=$P --rk=0 --acc=1e-6 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=120 --rbf_kernel=9 --rad=0.6 --numobj=1 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus44932.txt --csolve --solve

srun numactl --interleave=all ./build/timing/time_zpotrf_tile --m=44932 --n_range=4500:4500 --k=44932 --mb=1000 --nb=100 --nowarmup --threads=$TH  --p=$P --rk=0 --acc=1e-7 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=200 --rbf_kernel=9 --rad=0.6 --numobj=1 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus44932.txt --csolve --solve

#virus of 117715 resolution using 1e-5, 1e-6, and 1e-7

srun numactl --interleave=all ./build/timing/time_zpotrf_tile --m=117715 --n_range=3250:3250 --k=117715 --mb=1811 --nb=50 --nowarmup --threads=$TH  --p=$P --rk=0 --acc=1e-5 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=100  --rbf_kernel=9 --rad=0.6 --numobj=1 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus117715.txt --csolve --solve 

srun numactl --interleave=all ./build/timing/time_zpotrf_tile --m=117715 --n_range=3900:3900 --k=117715 --mb=1811 --nb=60 --nowarmup --threads=$TH  --p=$P --rk=0 --acc=1e-6 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=120  --rbf_kernel=9 --rad=0.6 --numobj=1 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus117715.txt --csolve --solve

srun numactl --interleave=all ./build/timing/time_zpotrf_tile --m=117715 --n_range=6500:6500 --k=117715 --mb=1811 --nb=100 --nowarmup --threads=$TH  --p=$P --rk=0 --acc=1e-7 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=200 --rbf_kernel=9 --rad=0.6 --numobj=1 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus117715.txt --csolve --solve 


#virus of 142418 resolution using 1e-5, 1e-6, and 1e-7

srun numactl --interleave=all ./build/timing/time_zpotrf_tile --m=142418 --n_range=3000:3000 --k=142418 --mb=2400 --nb=50 --nowarmup --threads=$TH  --p=$P --rk=0 --acc=1e-5 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=50 --rbf_kernel=9 --rad=0.6 --numobj=1 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus142418.txt --csolve --solve

srun numactl --interleave=all ./build/timing/time_zpotrf_tile --m=142418 --n_range=3600:3600 --k=142418 --mb=2400 --nb=60 --nowarmup --threads=$TH  --p=$P --rk=0 --acc=1e-6 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=120  --rbf_kernel=9 --rad=0.6 --numobj=1 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus142418.txt --csolve --solve

srun numactl --interleave=all ./build/timing/time_zpotrf_tile --m=142418 --n_range=6000:6000 --k=142418 --mb=2400 --nb=100 --nowarmup --threads=$TH  --p=$P --rk=0 --acc=1e-7 --m-3D-rbf --starshwavek=0 --starshdecay=0 --starshmaxrank=200 --rbf_kernel=9 --rad=0.6 --numobj=1 --order=2 --mesh_file=stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus142418.txt --csolve --solve
