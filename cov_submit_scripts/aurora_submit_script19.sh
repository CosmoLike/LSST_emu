#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -q shortq
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=3:00:00
#PBS -J 1-6670
#PBS -N LY10_19
#PBS -e /aurora_nobackup/sunglass/teifler/output/
#PBS -o /aurora_nobackup/sunglass/teifler/output/

cd $PBS_O_WORKDIR
/home/teifler/CosmoLike/LSSTC_emu/./compute_covariances_fourier19 $PBS_ARRAY_INDEX >& /aurora_nobackup/sunglass/teifler/job_output.log


