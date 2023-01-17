#module load mkl/2018-initial
for m in ecrc-extras mkl/2020.0.166 gcc/10.2.0 cmake/3.19.2 hwloc/2.4.0-nocuda-gcc-10.2.0 starpu/1.2.10-gcc-10.2.0-mkl-openmpi-4.1.0 ; do
	module load $m
done
#starpu/1.3.7-gcc-10.2.0-mkl-openmpi-4.1.0 

