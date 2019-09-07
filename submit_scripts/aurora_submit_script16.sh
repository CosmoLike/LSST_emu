#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -q verylongq
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=12:00:00
#PBS -J 150001-160000
#PBS -N LY10_150-160k
#PBS -e /aurora_nobackup/sunglass/teifler/output/
#PBS -o /aurora_nobackup/sunglass/teifler/output/

cd $PBS_O_WORKDIR
/home/teifler/CosmoLike/LSSTC_emu/./compute_covariances_fourier $PBS_ARRAY_INDEX >& /aurora_nobackup/sunglass/teifler/job_output.log


