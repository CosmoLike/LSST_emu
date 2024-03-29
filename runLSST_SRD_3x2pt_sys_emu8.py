#!/usr/bin/env python

import sys
sys.path.append('/home/u17/timeifler/LSST_emu')

from cosmolike_libs_emu import * 
from schwimmbad import MPIPool

inv=['inv_LSST_Y10_area7.623220e+03_ng1.399100e+01_nl2.297260e+01','inv_LSST_Y10_area1.478630e+04_ng3.339750e+01_nl6.159480e+01','inv_LSST_Y10_area9.931470e+03_ng1.706900e+01_nl2.878090e+01','inv_LSST_Y10_area8.585430e+03_ng2.848750e+01_nl5.143540e+01','inv_LSST_Y10_area1.768180e+04_ng3.536430e+01_nl6.572260e+01','inv_LSST_Y10_area1.512690e+04_ng1.038020e+01_nl1.637770e+01','inv_LSST_Y10_area9.747990e+03_ng1.111050e+01_nl1.768990e+01','inv_LSST_Y10_area8.335080e+03_ng2.959040e+01_nl5.369840e+01','inv_LSST_Y10_area9.533420e+03_ng3.146080e+01_nl5.756200e+01','inv_LSST_Y10_area1.833130e+04_ng1.574630e+01_nl2.626600e+01','inv_LSST_Y10_area1.286780e+04_ng1.330660e+01_nl2.170300e+01','inv_LSST_Y10_area1.741890e+04_ng9.072180e+00_nl1.405870e+01','inv_LSST_Y10_area1.978310e+04_ng9.587190e+00_nl1.496670e+01','inv_LSST_Y10_area1.253880e+04_ng1.632340e+01_nl2.736000e+01','inv_LSST_Y10_area1.526000e+04_ng2.583270e+01_nl4.603640e+01','inv_LSST_Y10_area1.654070e+04_ng1.074150e+01_nl1.702530e+01','inv_LSST_Y10_area1.963680e+04_ng3.804120e+01_nl7.139000e+01','inv_LSST_Y10_area1.111270e+04_ng3.288210e+01_nl6.051840e+01','inv_LSST_Y10_area1.038550e+04_ng1.988930e+01_nl3.422830e+01','inv_LSST_Y10_area1.614020e+04_ng2.703690e+01_nl4.847670e+01','inv_LSST_Y10_area1.892010e+04_ng1.517110e+01_nl2.518120e+01','inv_LSST_Y10_area1.797620e+04_ng1.424180e+01_nl2.344000e+01','inv_LSST_Y10_area1.135200e+04_ng1.918510e+01_nl3.285780e+01','inv_LSST_Y10_area9.214770e+03_ng2.692600e+01_nl4.825130e+01','inv_LSST_Y10_area1.691070e+04_ng2.200120e+01_nl3.837660e+01','inv_LSST_Y10_area1.199560e+04_ng1.265530e+01_nl2.050290e+01','inv_LSST_Y10_area1.619980e+04_ng1.863040e+01_nl3.178300e+01','inv_LSST_Y10_area1.439510e+04_ng1.197870e+01_nl1.926470e+01','inv_LSST_Y10_area8.133860e+03_ng3.687280e+01_nl6.890940e+01','inv_LSST_Y10_area1.351050e+04_ng2.252650e+01_nl3.941680e+01','inv_LSST_Y10_area1.912230e+04_ng1.743810e+01_nl2.948740e+01','inv_LSST_Y10_area1.568450e+04_ng1.234240e+01_nl1.992920e+01','inv_LSST_Y10_area1.201480e+04_ng1.000950e+01_nl1.571620e+01','inv_LSST_Y10_area1.405970e+04_ng2.359790e+01_nl4.154870e+01','inv_LSST_Y10_area1.091930e+04_ng2.087710e+01_nl3.616160e+01','inv_LSST_Y10_area1.321270e+04_ng2.489560e+01_nl4.414800e+01']

data=['LSST_3x2pt_Y10_area7.623220e+03','LSST_3x2pt_Y10_area1.478630e+04', 'LSST_3x2pt_Y10_area9.931470e+03', 'LSST_3x2pt_Y10_area8.585430e+03', 'LSST_3x2pt_Y10_area1.768180e+04', 'LSST_3x2pt_Y10_area1.512690e+04', 'LSST_3x2pt_Y10_area9.747990e+03', 'LSST_3x2pt_Y10_area8.335080e+03', 'LSST_3x2pt_Y10_area9.533420e+03', 'LSST_3x2pt_Y10_area1.833130e+04', 'LSST_3x2pt_Y10_area1.286780e+04', 'LSST_3x2pt_Y10_area1.741890e+04', 'LSST_3x2pt_Y10_area1.978310e+04', 'LSST_3x2pt_Y10_area1.253880e+04', 'LSST_3x2pt_Y10_area1.526000e+04', 'LSST_3x2pt_Y10_area1.654070e+04', 'LSST_3x2pt_Y10_area1.963680e+04', 'LSST_3x2pt_Y10_area1.111270e+04', 'LSST_3x2pt_Y10_area1.038550e+04', 'LSST_3x2pt_Y10_area1.614020e+04', 'LSST_3x2pt_Y10_area1.892010e+04', 'LSST_3x2pt_Y10_area1.797620e+04', 'LSST_3x2pt_Y10_area1.135200e+04', 'LSST_3x2pt_Y10_area9.214770e+03', 'LSST_3x2pt_Y10_area1.691070e+04', 'LSST_3x2pt_Y10_area1.199560e+04', 'LSST_3x2pt_Y10_area1.619980e+04', 'LSST_3x2pt_Y10_area1.439510e+04', 'LSST_3x2pt_Y10_area8.133860e+03', 'LSST_3x2pt_Y10_area1.351050e+04', 'LSST_3x2pt_Y10_area1.912230e+04', 'LSST_3x2pt_Y10_area1.568450e+04', 'LSST_3x2pt_Y10_area1.201480e+04', 'LSST_3x2pt_Y10_area1.405970e+04', 'LSST_3x2pt_Y10_area1.091930e+04', 'LSST_3x2pt_Y10_area1.321270e+04']

bary=['LPC_area7.623220e+03','LPC_area1.478630e+04', 'LPC_area9.931470e+03', 'LPC_area8.585430e+03', 'LPC_area1.768180e+04', 'LPC_area1.512690e+04', 'LPC_area9.747990e+03', 'LPC_area8.335080e+03', 'LPC_area9.533420e+03', 'LPC_area1.833130e+04', 'LPC_area1.286780e+04', 'LPC_area1.741890e+04', 'LPC_area1.978310e+04', 'LPC_area1.253880e+04', 'LPC_area1.526000e+04', 'LPC_area1.654070e+04', 'LPC_area1.963680e+04', 'LPC_area1.111270e+04', 'LPC_area1.038550e+04', 'LPC_area1.614020e+04', 'LPC_area1.892010e+04', 'LPC_area1.797620e+04', 'LPC_area1.135200e+04', 'LPC_area9.214770e+03', 'LPC_area1.691070e+04', 'LPC_area1.199560e+04', 'LPC_area1.619980e+04', 'LPC_area1.439510e+04', 'LPC_area8.133860e+03', 'LPC_area1.351050e+04', 'LPC_area1.912230e+04', 'LPC_area1.568450e+04', 'LPC_area1.201480e+04', 'LPC_area1.405970e+04', 'LPC_area1.091930e+04', 'LPC_area1.321270e+04']

source_z=['wl_redshift_model0_WLz01.880307e-01_WLalpha8.485694e-01.txt', 'wl_redshift_model1_WLz01.731166e-01_WLalpha7.662434e-01.txt', 'wl_redshift_model2_WLz01.846221e-01_WLalpha8.297540e-01.txt', 'wl_redshift_model3_WLz01.758423e-01_WLalpha7.812893e-01.txt', 'wl_redshift_model4_WLz01.721357e-01_WLalpha7.608290e-01.txt', 'wl_redshift_model5_WLz01.931476e-01_WLalpha8.768147e-01.txt', 'wl_redshift_model6_WLz01.919821e-01_WLalpha8.703814e-01.txt', 'wl_redshift_model7_WLz01.751912e-01_WLalpha7.776952e-01.txt', 'wl_redshift_model8_WLz01.741405e-01_WLalpha7.718958e-01.txt', 'wl_redshift_model9_WLz01.860048e-01_WLalpha8.373865e-01.txt', 'wl_redshift_model10_WLz01.888904e-01_WLalpha8.533151e-01.txt', 'wl_redshift_model11_WLz01.954564e-01_WLalpha8.895592e-01.txt', 'wl_redshift_model12_WLz01.945099e-01_WLalpha8.843347e-01.txt', 'wl_redshift_model13_WLz01.853877e-01_WLalpha8.339802e-01.txt', 'wl_redshift_model14_WLz01.775192e-01_WLalpha7.905458e-01.txt', 'wl_redshift_model15_WLz01.925611e-01_WLalpha8.735775e-01.txt', 'wl_redshift_model16_WLz01.708849e-01_WLalpha7.539247e-01.txt', 'wl_redshift_model17_WLz01.733832e-01_WLalpha7.677150e-01.txt', 'wl_redshift_model18_WLz01.820009e-01_WLalpha8.152851e-01.txt', 'wl_redshift_model19_WLz01.767381e-01_WLalpha7.862345e-01.txt', 'wl_redshift_model20_WLz01.866426e-01_WLalpha8.409073e-01.txt', 'wl_redshift_model21_WLz01.877261e-01_WLalpha8.468882e-01.txt', 'wl_redshift_model22_WLz01.826188e-01_WLalpha8.186960e-01.txt', 'wl_redshift_model23_WLz01.768086e-01_WLalpha7.866234e-01.txt', 'wl_redshift_model24_WLz01.802711e-01_WLalpha8.057362e-01.txt', 'wl_redshift_model25_WLz01.897506e-01_WLalpha8.580632e-01.txt', 'wl_redshift_model26_WLz01.831217e-01_WLalpha8.214720e-01.txt', 'wl_redshift_model27_WLz01.906926e-01_WLalpha8.632630e-01.txt', 'wl_redshift_model28_WLz01.714197e-01_WLalpha7.568766e-01.txt', 'wl_redshift_model29_WLz01.798667e-01_WLalpha8.035040e-01.txt', 'wl_redshift_model30_WLz01.842554e-01_WLalpha8.277299e-01.txt', 'wl_redshift_model31_WLz01.901797e-01_WLalpha8.604322e-01.txt', 'wl_redshift_model32_WLz01.937710e-01_WLalpha8.802561e-01.txt', 'wl_redshift_model33_WLz01.790701e-01_WLalpha7.991072e-01.txt', 'wl_redshift_model34_WLz01.811701e-01_WLalpha8.106987e-01.txt', 'wl_redshift_model35_WLz01.781525e-01_WLalpha7.940419e-01.txt'] 

lens_z=['LSS_redshift_model0_LSSz02.629496e-01_LSSalpha9.285983e-01.txt','LSS_redshift_model1_LSSz02.852923e-01_LSSalpha8.985943e-01.txt','LSS_redshift_model2_LSSz02.664822e-01_LSSalpha9.186036e-01.txt','LSS_redshift_model3_LSSz02.798758e-01_LSSalpha9.014201e-01.txt','LSS_redshift_model4_LSSz02.873873e-01_LSSalpha8.978683e-01.txt','LSS_redshift_model5_LSSz02.593970e-01_LSSalpha9.470921e-01.txt','LSS_redshift_model6_LSSz02.600213e-01_LSSalpha9.425114e-01.txt','LSS_redshift_model7_LSSz02.811155e-01_LSSalpha9.006370e-01.txt','LSS_redshift_model8_LSSz02.831875e-01_LSSalpha8.995165e-01.txt','LSS_redshift_model9_LSSz02.649368e-01_LSSalpha9.224338e-01.txt','LSS_redshift_model10_LSSz02.622058e-01_LSSalpha9.314128e-01.txt','LSS_redshift_model11_LSSz02.584820e-01_LSSalpha9.568082e-01.txt','LSS_redshift_model12_LSSz02.588053e-01_LSSalpha9.527220e-01.txt','LSS_redshift_model13_LSSz02.656075e-01_LSSalpha9.206866e-01.txt','LSS_redshift_model14_LSSz02.768397e-01_LSSalpha9.037492e-01.txt','LSS_redshift_model15_LSSz02.596975e-01_LSSalpha9.447600e-01.txt','LSS_redshift_model16_LSSz02.901709e-01_LSSalpha8.971658e-01.txt','LSS_redshift_model17_LSSz02.847362e-01_LSSalpha8.988183e-01.txt','LSS_redshift_model18_LSSz02.698330e-01_LSSalpha9.121821e-01.txt','LSS_redshift_model19_LSSz02.782257e-01_LSSalpha9.026084e-01.txt','LSS_redshift_model20_LSSz02.642756e-01_LSSalpha9.243038e-01.txt','LSS_redshift_model21_LSSz02.632273e-01_LSSalpha9.276296e-01.txt','LSS_redshift_model22_LSSz02.689934e-01_LSSalpha9.135968e-01.txt','LSS_redshift_model23_LSSz02.780987e-01_LSSalpha9.027073e-01.txt','LSS_redshift_model24_LSSz02.723464e-01_LSSalpha9.085463e-01.txt','LSS_redshift_model25_LSSz02.615210e-01_LSSalpha9.343470e-01.txt','LSS_redshift_model26_LSSz02.683327e-01_LSSalpha9.147934e-01.txt','LSS_redshift_model27_LSSz02.608392e-01_LSSalpha9.376962e-01.txt','LSS_redshift_model28_LSSz02.889655e-01_LSSalpha8.974355e-01.txt','LSS_redshift_model29_LSSz02.729686e-01_LSSalpha9.077654e-01.txt','LSS_redshift_model30_LSSz02.669178e-01_LSSalpha9.176391e-01.txt','LSS_redshift_model31_LSSz02.612016e-01_LSSalpha9.358553e-01.txt','LSS_redshift_model32_LSSz02.591077e-01_LSSalpha9.496317e-01.txt','LSS_redshift_model33_LSSz02.742326e-01_LSSalpha9.063038e-01.txt','LSS_redshift_model34_LSSz02.710103e-01_LSSalpha9.103760e-01.txt','LSS_redshift_model35_LSSz02.757517e-01_LSSalpha9.047459e-01.txt']

area_table=[7623.22,14786.3,9931.47,8585.43,17681.8,15126.9,9747.99,8335.08,9533.42,18331.3,12867.8,17418.9,19783.1,12538.8,15260.0,16540.7,19636.8,11112.7,10385.5,16140.2,18920.1,17976.2,11352.0,9214.77,16910.7,11995.6,16199.8,14395.1,8133.86,13510.5,19122.3,15684.5,12014.8,14059.7,10919.3,13212.7]

nsource_table=[13.991,33.3975,17.069,28.4875,35.3643,10.3802,11.1105,29.5904,31.4608,15.7463,13.3066,9.07218,9.58719,16.3234,25.8327,10.7415,38.0412,32.8821,19.8893,27.0369,15.1711,14.2418,19.1851,26.926,22.0012,12.6553,18.6304,11.9787,36.8728,22.5265,17.4381,12.3424,10.0095,23.5979,20.8771,24.8956]
  
nlens_table=[22.9726 ,61.5948 ,28.7809 ,51.4354 ,65.7226 ,16.3777 ,17.6899 ,53.6984 ,57.562 ,26.266 ,21.703 ,14.0587 ,14.9667 ,27.36 ,46.0364 ,17.0253 ,71.39 ,60.5184 ,34.2283 ,48.4767 ,25.1812 ,23.44 ,32.8578 ,48.2513 ,38.3766 ,20.5029 ,31.783 ,19.2647 ,68.9094 ,39.4168 ,29.4874 ,19.9292 ,15.7162 ,41.5487 ,36.1616 ,44.148]
  
shear_prior=[0.00891915,0.0104498 ,0.0145972 ,0.0191916 ,0.00450246 ,0.00567828 ,0.00294841 ,0.00530922 ,0.0118632 ,0.0151849 ,0.00410151 ,0.0170622 ,0.0197331 ,0.0106615 ,0.0124445 ,0.00994507 ,0.0136251 ,0.0143491 ,0.0164314 ,0.016962 ,0.0186608 ,0.00945903 ,0.0113246 ,0.0155225 ,0.00800846 ,0.00732104 ,0.00649453 ,0.00243976 ,0.0125932 ,0.0182587 ,0.00335859 ,0.00682287 ,0.0177269 ,0.0035219 ,0.00773304 ,0.0134886] 

delta_z_prior=[0.0032537,0.00135316,0.00168787,0.00215043,0.00406031,0.00222358,0.00334993,0.00255186,0.00266499,0.00159226,0.00183664,0.00384965,0.00427765,0.00314377,0.00456113,0.00347868,0.00487938,0.00418152,0.00469911,0.00367598,0.0028009,0.00234161,0.00194964,0.00200982,0.00122739,0.00310886,0.00275168,0.00492736,0.00437241,0.00113931,0.00104864,0.00292328,0.00452082,0.00394114,0.00150756,0.003613]

sigma_z_prior=[0.00331909,0.00529541,0.00478151,0.00437497,0.00443062,0.00486333,0.00467423,0.0036723,0.00426963,0.00515357,0.0054553,0.00310132,0.00305971,0.00406327,0.00594293,0.00348709,0.00562526,0.00396025,0.00540537,0.00500447,0.00318595,0.00460592,0.00412137,0.00336418,0.00524988,0.00390092,0.00498349,0.0056667,0.0036384,0.00455861,0.00554822,0.00381061,0.0057615,0.00357705,0.00590572,0.00422393]

sigma_z=[0.0849973, 0.0986032, 0.0875521, 0.0968222, 0.0225239, 0.0718278, 0.0733675, 0.0385274, 0.0425549, 0.0605867, 0.0178555, 0.0853407, 0.0124119, 0.0531027, 0.0304032, 0.0503145, 0.0132213, 0.0941765, 0.0416444, 0.0668198, 0.063227, 0.0291332, 0.0481633, 0.0595606, 0.0818742, 0.0472518, 0.0270185, 0.0767401, 0.0219945, 0.0902663, 0.0779705, 0.0337666, 0.0362358, 0.0692429, 0.0558841, 0.0150457]

survey_designation=['LSST']
tomo_binning_source=['source_std']
tomo_binning_lens=['LSST_gold']

model=7
file_source_z = os.path.join(dirname, "zdistris/",source_z[model])
file_lens_z = os.path.join(dirname, "zdistris/",lens_z[model])
data_file = os.path.join(dirname, "datav/",data[model]) 
cov_file = os.path.join(dirname, "inv/",inv[model])
#cov_file = os.path.join("/Users/timeifler/Dropbox/cosmolike_store/LSST_emu/inv/",inv[model])
#chain_file = os.path.join(dirname, "chains/like_LSST_3x2pt_sys_model_%d" %model)
chain_file = os.path.join(dirname, "/groups/timeifler/LSST_emu/chains/LSST_3x2pt_sys_model_%d" %model)
bary_file=os.path.join(dirname, "baryons/",bary[model])

initcosmo("halofit")
initbins(15,20.0,3000.0,3000.0,21.0,10,10)
initpriors(shear_prior[model],sigma_z[model],delta_z_prior[model],sigma_z_prior[model],sigma_z[model]*0.6,delta_z_prior[model],sigma_z_prior[model],3.0,1.2,3.8,2.0,16.0,5.0,0.8)
initsurvey(survey_designation[0],nsource_table[model],nlens_table[model],area_table[model])
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","source_std","LSST_gold")
initia("NLA_HF","GAMA")

# test also with
#initpriors("none","none","none","Planck")
#initpriors("none","none","none","random")
initprobes("3x2pt")
initdatainv(cov_file ,data_file, bary_file)

#sample_params=sample_LCDM_only()
#sample_params= sample_cosmology_only()
#sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear())
#sample_params = sample_cosmology_2pt_shear_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_photo_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
sample_params = sample_cosmology_2pt_nuisance_IA_bary_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 

sample_main(sample_params,sigma_z[model],sigma_z_prior[model],5000,840,1,chain_file, blind=False, pool=MPIPool())

