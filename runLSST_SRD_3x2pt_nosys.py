#!/usr/bin/env python

import sys
sys.path.append('/home/u17/timeifler/LSSTC_emu')

from cosmolike_libs_pessi import * 
from schwimmbad import MPIPool

file_source_z = os.path.join(dirname, "zdistris/SRD_zdistri_model_z0=1.100000e-01_beta=6.800000e-01_Y10_source.txt")
file_lens_z = os.path.join(dirname, "zdistris/SRD_zdistri_model_z0=2.800000e-01_beta=9.000000e-01_Y10_lens.txt")
data_file = os.path.join(dirname, "datav/LSST_3x2pt_fid_pessi")
cov_file = os.path.join(dirname, "cov/LSST_SRD_cov_3x2pt_inv")
chain_file = os.path.join(dirname, "like/like_LSST_3x2pt_nosys")

initcosmo("halofit")
initbins(20,20.0,15000.0,3000.0,21.0,10,10)
initpriors("pessi","none","none","none") 
initsurvey("LSST")
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","LSST_Y10")
initclusters()
initia("none","none")

# test also with
#initpriors("none","none","none","Planck")
#initpriors("none","none","none","random")
initprobes("3x2pt")
initdatainv(cov_file ,data_file)

#sample_params=sample_LCDM_only()
sample_params= sample_cosmology_only()
#sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear())
#sample_params = sample_cosmology_2pt_shear_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_photo_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 

sample_main(sample_params,1000,560,1,chain_file, blind=False, pool=MPIPool())

