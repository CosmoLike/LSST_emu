#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -q mediumq
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=12:00:00
#PBS -J 1-7140
#PBS -N LY10_21
#PBS -e /halo_nobackup/sunglass/teifler/output/
#PBS -o /halo_nobackup/sunglass/teifler/output/

cd $PBS_O_WORKDIR
/home/teifler/CosmoLike/LSSTC_emu/./compute_covariances_fourier21 $PBS_ARRAY_INDEX > /halo_nobackup/sunglass/teifler/job_output.log


