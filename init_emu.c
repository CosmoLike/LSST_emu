double invcov_read(int READ, int ci, int cj);
double data_read(int READ, int ci);
double bary_read(int READ, int PC, int cj);
void init_data_inv_bary(char *INV_FILE, char *DATA_FILE, char *BARY_FILE);
void init_priors(double Prior1, double Prior2, double Prior3, double Prior4);
void init_survey(char *surveyname);
void init_galaxies(char *SOURCE_ZFILE, char *LENS_ZFILE, char *lensphotoz, char *sourcephotoz, char *galsample);
void init_cosmo_runmode(char *runmode);
void init_binning_fourier(int Ncl, double lmin, double lmax, double lmax_shear, double Rmin_bias, int Ntomo_source, int Ntomo_lens);
void init_probes(char *probes);

void set_lens_galaxies_LSST();
void init_lens_sample(char *lensphotoz, char *galsample);
void init_source_sample(char *sourcephotoz);

void set_galaxies_source();
void set_galaxies_lens();

void set_clusters_LSST(); //set parameters for LSST/WFIRST forecasts

void init_clusterMobs();
void set_equal_tomo_bins();
void init_IA(char *model,char *lumfct);


int count_rows(char* filename,const char delimiter){
  FILE *file = fopen (filename, "r" );
  char line [1000];
  if (file != NULL) {
    fgets(line,sizeof line,file);
    fclose(file);
  }
  else{
    printf("count_rows: file %s not found.\nEXIT\n",filename);
    exit(1);
  }
  int count = 1;
  char *p;

  p = line;
  while (*p != '\0')
  {
    if (*p == delimiter){
        while (*p == delimiter){p++;}
        count++;
    }
      p++;
    }
   return count;
}


void init_data_inv_bary(char *INV_FILE, char *DATA_FILE, char *BARY_FILE)
{
  double init;
  printf("\n");
  printf("---------------------------------------\n");
  printf("Initializing data vector and covariance\n");
  printf("---------------------------------------\n");

  sprintf(like.INV_FILE,"%s",INV_FILE);
  printf("PATH TO INVCOV: %s\n",like.INV_FILE);
  sprintf(like.DATA_FILE,"%s",DATA_FILE);
  printf("PATH TO DATA: %s\n",like.DATA_FILE);
  sprintf(like.BARY_FILE,"%s",BARY_FILE);
  printf("PATH TO BARYONS: %s\n",like.BARY_FILE);
  init=data_read(0,1);
  init=bary_read(0,1,1);
  init=invcov_read(0,1,1);
}

double invcov_read(int READ, int ci, int cj)
{
  int i,j,intspace;
  static double **inv =0;

  if(READ==0 || inv == 0){
    inv   = create_double_matrix(0, like.Ndata-1, 0, like.Ndata-1);      
    FILE *F;
    F=fopen(like.INV_FILE,"r");
    for (i=0;i<like.Ndata; i++){
      for (j=0;j<like.Ndata; j++){
       fscanf(F,"%d %d %le\n",&intspace,&intspace,&inv[i][j]);  
     }
   }
   fclose(F);
   printf("FINISHED READING COVARIANCE\n");
 }    
 return inv[ci][cj];
}


double data_read(int READ, int ci)
{
  int i,intspace;
  static double *data = 0;
  
  if(READ==0 || data ==0){
    data  = create_double_vector(0, like.Ndata-1);      
    FILE *F;
    F=fopen(like.DATA_FILE,"r");
    for (i=0;i<like.Ndata; i++){  
      fscanf(F,"%d %le\n",&intspace,&data[i]);
    }
    fclose(F);
    printf("FINISHED READING DATA VECTOR\n");
  }    
  return data[ci];
}

double bary_read(int READ, int PC, int cj)
{
  int i,j,intspace, N_PC=6;
  static double **bary =0;

  if(READ==0 || bary == 0){
    bary   = create_double_matrix(0, N_PC-1, 0, like.Ndata-1);      
    FILE *F;
    F=fopen(like.BARY_FILE,"r");
    for (i=0;i<like.Ndata; i++){
      fscanf(F,"%d %le %le %le %le %le %le\n",&intspace,&bary[0][i],&bary[1][i],&bary[2][i],&bary[3][i],&bary[4][i],&bary[5][i]);  
    }
    fclose(F);
    printf("FINISHED READING BARYON MATRIX\n");
  }    
  return bary[PC][cj];
}



void init_cosmo_runmode(char *runmode)
{
  printf("\n");
  printf("-------------------------------------------\n");
  printf("Initializing Standard Runmode/Cosmology\n");
  printf("-------------------------------------------\n");
  set_cosmological_parameters_to_Planck_15_TT_TE_EE_lowP();
  //set_cosmological_parameters_to_Joe();
  sprintf(pdeltaparams.runmode,"%s",runmode);
  printf("pdeltaparams.runmode =%s\n",pdeltaparams.runmode);
}

void init_binning_fourier(int Ncl, double lmin, double lmax, double lmax_shear, double Rmin_bias, int Ntomo_source, int Ntomo_lens)
{
  printf("-------------------------------------------\n");
  printf("Initializing Binning\n");
  printf("-------------------------------------------\n");
  
  like.Rmin_bias=Rmin_bias;
  like.Ncl=Ncl;
  like.lmin= lmin; //std=20
  like.lmax= lmax; //15,000
  like.lmax_shear = lmax_shear; //5000
  tomo.shear_Nbin=Ntomo_source;
  tomo.clustering_Nbin=Ntomo_lens;
  //compute cluster ell bins acc to 2PCF l-bins
  double ell;
  int i,k=0;
  double logdl=(log(like.lmax)-log(like.lmin))/like.Ncl;
  for(i=0;i<like.Ncl;i++){
    ell=exp(log(like.lmin)+(i+0.5)*logdl);
    if (ell > like.lmax_shear){
      if (k==0) Cluster.l_min = ell;
      k=k+1;
    }
  } 
  Cluster.lbin = k;
  Cluster.l_max = lmax; //clusters go to highly nonlin as std
  like.lmax_kappacmb = 2999.;
  
  printf("number of ell bins Ncl: %d\n",like.Ncl);
  printf("minimum ell: %le\n",like.lmin);
  printf("maximum ell: %le\n",like.lmax);
}


void init_priors(double Prior1, double Prior2, double Prior3, double Prior4)
{
  int i;

  printf("\n");
  printf("---------------------------------------\n");
  printf("Initializing observational priors for marginalization\n");
  printf("---------------------------------------\n");
  
  for (i=0;i<tomo.shear_Nbin; i++){
    prior.shear_calibration_m[i][0] = 0.0;
    prior.shear_calibration_m[i][1] = Prior1;
  }
  like.shearcalib=1;
  
  for (i=0;i<tomo.shear_Nbin; i++){
    nuisance.bias_zphot_shear[i]=0.0;
    nuisance.sigma_zphot_shear[i]=Prior2; 
    // center of Gaussian priors
    prior.bias_zphot_shear[i][0]=nuisance.bias_zphot_shear[i];
    prior.sigma_zphot_shear[i][0]=nuisance.sigma_zphot_shear[i]; 
    // rms width of Gaussian priors
    prior.bias_zphot_shear[i][1] = Prior3;
    prior.sigma_zphot_shear[i][1]= Prior4;
  }
  like.wlphotoz=1;

  for (i=0;i<tomo.clustering_Nbin; i++){
    nuisance.bias_zphot_clustering[i]=0.0;
    nuisance.sigma_zphot_clustering[i]=Prior2*0.6; 
    // center of Gaussian priors
    prior.bias_zphot_clustering[i][0]=nuisance.bias_zphot_clustering[i];
    prior.sigma_zphot_clustering[i][0]=nuisance.sigma_zphot_clustering[i];
    // rms width of Gaussian priors
    prior.bias_zphot_clustering[i][1] = Prior3/2.;
    prior.sigma_zphot_clustering[i][1]= Prior4/2.;
  }
  like.clphotoz=1;

printf("\n");
printf("---------------------------------------\n");
printf("Shear Calibration Prior\n");
printf("Mean=%le, Sigma=%le\n",prior.shear_calibration_m[0][0],Prior1);

printf("\n");
printf("---------------------------------------\n");
printf("Photo-z priors Weak Lensing\n");
printf("Delta_z=%le, Sigma (Delta_z)=%le\n",prior.bias_zphot_shear[0][0],prior.bias_zphot_shear[0][1]);
printf("Sigma_z=%le, Sigma (Sigma_z)=%le\n",prior.sigma_zphot_shear[0][0],prior.sigma_zphot_shear[0][1]); 

printf("\n");
printf("---------------------------------------\n");
printf("Photo-z priors Clustering\n");
printf("Delta_z=%le, Sigma (Delta_z)=%le\n",prior.bias_zphot_clustering[0][0],prior.bias_zphot_clustering[0][1]);
printf("Sigma_z=%le, Sigma (Sigma_z)=%le\n",prior.sigma_zphot_clustering[0][0],prior.sigma_zphot_clustering[0][1]);

}




 

void init_survey(char *surveyname)
{
  printf("\n");
  printf("-------------------------------\n");
  printf("Initializing Survey Parameters\n");
  printf("-------------------------------\n");

  if(strcmp(surveyname,"LSST")==0) set_survey_parameters_to_LSST();
  if(strcmp(surveyname,"Euclid")==0) set_survey_parameters_to_Euclid();
  if(strcmp(surveyname,"WFIRST")==0) set_survey_parameters_to_WFIRST();
  printf("Survey set to %s\n",survey.name);
  printf("Survey area: %le deg^2\n",survey.area);
  printf("Source Galaxy Density: %le galaxies/arcmin^2\n",survey.n_gal); 
}


void init_galaxies(char *SOURCE_ZFILE, char *LENS_ZFILE, char *lensphotoz, char *sourcephotoz, char *galsample)
{
  printf("\n");
  printf("-----------------------------------\n");
  printf("Initializing galaxy samples\n");
  printf("-----------------------------------\n");
  
  sprintf(redshift.shear_REDSHIFT_FILE,"%s",SOURCE_ZFILE);
  printf("PATH TO SOURCE_ZFILE: %s\n",redshift.shear_REDSHIFT_FILE);
  
  init_source_sample(sourcephotoz);
  
  sprintf(redshift.clustering_REDSHIFT_FILE,"%s",LENS_ZFILE);
  printf("\n");
  printf("PATH TO LENS_ZFILE: %s\n",redshift.clustering_REDSHIFT_FILE);
  init_lens_sample(lensphotoz,galsample);
}

void init_clusters()
{
  printf("\n");
  printf("-----------------------------------\n");
  printf("Initializing clusters\n");
  printf("-----------------------------------\n");

  if (strcmp(survey.name,"LSST")==0 || strcmp(survey.name,"WFIRST")==0 || strcmp(survey.name,"HSC")==0) set_clusters_LSST();
  if (strcmp(survey.name,"Euclid")==0 || strcmp(survey.name,"DES")==0) {printf("Not in this init! Exit! \n"); exit(1);}

 set_clusterMobs_priors(); 
}


void init_probes(char *probes)
{
  printf("\n");
  printf("------------------------------\n");
  printf("Initializing Probes\n");
  printf("------------------------------\n"); 
  printf("tomo.cluster_Nbin=%d\n",tomo.cluster_Nbin);
  printf("tomo.cgl_Npowerspectra=%d\n",tomo.cgl_Npowerspectra);
  printf("Cluster.lbin=%d\n",Cluster.lbin);
  printf("Cluster.N200_Nbin=%d\n",Cluster.N200_Nbin);
  printf("like.Ncl=%d\n",like.Ncl);
  printf("tomo.shear_Npowerspectra=%d\n",tomo.shear_Npowerspectra);
  printf("tomo.ggl_Npowerspectra=%d\n",tomo.ggl_Npowerspectra);
  printf("tomo.clustering_Npowerspectra=%d\n",tomo.clustering_Npowerspectra);

  sprintf(like.probes,"%s",probes);
  if(strcmp(probes,"clusterN")==0){
    like.Ndata=tomo.cluster_Nbin*Cluster.N200_Nbin;
    like.clusterN=1;
    printf("Cluster Number Counts computation initialized\n");
  }
  if(strcmp(probes,"clusterN_clusterWL")==0){
    like.Ndata=tomo.cluster_Nbin*Cluster.N200_Nbin+tomo.cgl_Npowerspectra*Cluster.N200_Nbin*Cluster.lbin;
    like.clusterN=1;
    like.clusterWL=1;
    printf("Cluster Number Counts computation initialized\n");
    printf("Cluster weak lensing computation initialized\n");
  }
  if(strcmp(probes,"3x2pt_clusterN")==0){
    like.Ndata=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+tomo.cluster_Nbin*Cluster.N200_Nbin;
    like.shear_shear=1;
    like.shear_pos=1;
    like.pos_pos=1;
    like.clusterN=1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
    printf("Cluster Number Counts computation initialized\n");
  }

  if(strcmp(probes,"shear_shear")==0){
    like.Ndata=like.Ncl*tomo.shear_Npowerspectra;
    like.shear_shear=1;
    printf("Shear-Shear computation initialized\n");
  }
  if(strcmp(probes,"pos_pos")==0){
    like.Ndata= like.Ncl*tomo.clustering_Npowerspectra;
    like.pos_pos=1;
    printf("Position-Position computation initialized\n");
  }
  if(strcmp(probes,"ggl_cl")==0){
    like.Ndata=like.Ncl*(tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
    like.shear_pos=1;
    like.pos_pos=1;
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  }
  if(strcmp(probes,"3x2pt")==0){
    like.Ndata=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
    like.shear_shear=1;
    like.shear_pos=1;
    like.pos_pos=1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  } 
  if(strcmp(probes,"3x2pt_clusterN_clusterWL")==0){
    like.Ndata=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+tomo.cluster_Nbin*Cluster.N200_Nbin+tomo.cgl_Npowerspectra*Cluster.N200_Nbin*Cluster.lbin;
    like.shear_shear=1;
    like.shear_pos=1;
    like.pos_pos=1;
    like.clusterN=1;
    like.clusterWL=1;
    printf("%d\n",like.Ndata);
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
    printf("Cluster Number Counts computation initialized\n");
    printf("Cluster weak lensing computation initialized\n");
  }
   if (strcmp(probes,"LSSxCMB")==0) {
    like.Ndata = like.Ncl * (2*tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+tomo.shear_Npowerspectra);
    like.pos_pos = 1;
    like.gk = 1;
    like.shear_pos = 1;
    like.kk = 1;
    like.ks = 1;
    like.shear_shear = 1;
    printf("Initializing: gg, gk, gs, kk, ks, ss\n");
   }
   if (strcmp(probes,"gg_gk_gs")==0) {
    like.Ndata = like.Ncl * (2*tomo.clustering_Nbin+tomo.ggl_Npowerspectra);
    like.pos_pos = 1;
    like.gk = 1;
    like.shear_pos = 1;
    printf("Initializing: gg, gk, gs\n");
  }
  if (strcmp(probes,"kk_ks_ss")==0) {
    like.Ndata = like.Ncl * (1+tomo.shear_Nbin+tomo.shear_Npowerspectra);
    like.kk = 1;
    like.ks = 1;
    like.shear_shear = 1;
    printf("Initializing: kk, ks, ss\n");
  }
  printf("Total number of data points like.Ndata=%d\n",like.Ndata);
}



void init_lens_sample(char *lensphotoz, char *galsample)
{
  if(strcmp(lensphotoz,"none")==0) redshift.clustering_photoz=0;
  if(strcmp(lensphotoz,"voigt")==0) redshift.clustering_photoz=1;
  if(strcmp(lensphotoz,"voigt_out")==0) redshift.clustering_photoz=2;
  if(strcmp(lensphotoz,"gaussian")==0) redshift.clustering_photoz=3;
  if(strcmp(lensphotoz,"multihisto")==0) redshift.clustering_photoz=4;
  
  if ((redshift.clustering_photoz !=0) && (redshift.clustering_photoz !=1) && (redshift.clustering_photoz !=2) && (redshift.clustering_photoz !=3)) 
  {
    printf("init.c: init_lens_sample: redshift.clustering_photoz = %d not set properly!\nEXIT!\n",redshift.clustering_photoz);
    exit(1);
  }
  printf("Lens Sample Redshift Errors set to %s: redshift.clustering_photoz=%d\n",lensphotoz,redshift.clustering_photoz);
  
  if (strcmp(survey.name,"LSST")==0 || strcmp(survey.name,"WFIRST")==0 || strcmp(survey.name,"HSC")==0 || strcmp(survey.name,"Euclid")==0 ){ 
    if(strcmp(galsample,"source")==0){
      set_galaxies_lens();
    }
    if(strcmp(galsample,"LSST_Y10")==0){
      set_lens_galaxies_LSST();
    }
  }
  //call test_kmax once to initialize look-up tables at reference cosmology
  test_kmax(1000.,1);
}



void init_source_sample(char *sourcephotoz)
{
  if(strcmp(sourcephotoz,"none")==0) redshift.shear_photoz=0;
  if(strcmp(sourcephotoz,"voigt")==0) redshift.shear_photoz=1;
  if(strcmp(sourcephotoz,"voigt_out")==0) redshift.shear_photoz=2;
  if(strcmp(sourcephotoz,"gaussian")==0) redshift.shear_photoz=3;
  if(strcmp(sourcephotoz,"multihisto")==0) {
    printf("redshift.shear_photoz=4 not supported\n"); 
    exit(1);
  }
  if ((redshift.shear_photoz !=0) && (redshift.shear_photoz !=1) && (redshift.shear_photoz !=2) && (redshift.shear_photoz !=3)) 
  {
    printf("init.c: init_source_sample: redshift.shear_photoz = %d not set properly!\nEXIT!\n",redshift.shear_photoz);
    exit(1);
  }

  printf("Source Sample Redshift Errors set to %s: redshift.shear_photoz=%d\n",sourcephotoz,redshift.shear_photoz);
    set_galaxies_source();
}


void set_galaxies_source()
{
  int k,j;
  double frac, zi;
  
  tomo.shear_Npowerspectra=(int) (tomo.shear_Nbin*(tomo.shear_Nbin+1)/2);
  
  zdistr_histo_1(0.1, NULL);
  int zbins =2000;
  double da = (redshift.shear_zdistrpar_zmax-redshift.shear_zdistrpar_zmin)/(1.0*zbins);
  double *sum;
  sum=create_double_vector(0, zbins);
  
  sum[0] = 0.0;
  for (k = 0, zi = redshift.shear_zdistrpar_zmin; k<zbins; k++,zi+=da){
    sum[k+1] = sum[k]+zdistr_histo_1(zi, NULL);
  }
  
  tomo.shear_zmin[0] = redshift.shear_zdistrpar_zmin;
  tomo.shear_zmax[tomo.shear_Nbin-1] = redshift.shear_zdistrpar_zmax;
  printf("\n");
  printf("Source Sample - Tomographic Bin limits:\n");
  for(k=0;k<tomo.shear_Nbin-1;k++){
    frac=(k+1.)/(1.*tomo.shear_Nbin)*sum[zbins-1];
    j = 0;
    while (sum[j]< frac){
      j++;
    }
    tomo.shear_zmax[k] = redshift.shear_zdistrpar_zmin+j*da;
    tomo.shear_zmin[k+1] = redshift.shear_zdistrpar_zmin+j*da;
    printf("min=%le max=%le\n",tomo.shear_zmin[k],tomo.shear_zmax[k]);
  }
  printf("min=%le max=%le\n",tomo.shear_zmin[tomo.shear_Nbin-1],tomo.shear_zmax[tomo.shear_Nbin-1]);
  printf("redshift.shear_zdistrpar_zmin=%le max=%le\n",redshift.shear_zdistrpar_zmin,redshift.shear_zdistrpar_zmax);
  free_double_vector(sum,0,zbins);
}

void set_galaxies_lens()
{
  int k,j,n,i;
  double frac, zi;
  redshift.clustering_zdistrpar_zmin = 5.834332e-03;
  redshift.clustering_zdistrpar_zmax = 3.494167e+00;
  tomo.clustering_Npowerspectra=tomo.clustering_Nbin;
  tomo.clustering_zmin[0] = redshift.clustering_zdistrpar_zmin;
  tomo.clustering_zmax[tomo.clustering_Nbin-1] = redshift.clustering_zdistrpar_zmax;
  pf_histo(zi, NULL);
  int zbins =2000;
  double da = (redshift.clustering_zdistrpar_zmax-redshift.clustering_zdistrpar_zmin)/(1.0*zbins);
  double *sum;
  sum=create_double_vector(0, zbins);
  
  sum[0] = 0.0;
  for (k = 0, zi = redshift.clustering_zdistrpar_zmin; k<zbins; k++,zi+=da){
    sum[k+1] = sum[k]+pf_histo(zi, NULL);
  }
  
  printf("\n");
  printf("Lens Sample - Tomographic Bin limits:\n");
  for(k=0;k<tomo.clustering_Nbin-1;k++){
    frac=(k+1.)/(1.*tomo.clustering_Nbin)*sum[zbins-1];
    j = 0;
    while (sum[j]< frac){
      j++;
    }
    tomo.clustering_zmax[k] = redshift.clustering_zdistrpar_zmin+j*da;
    tomo.clustering_zmin[k+1] = redshift.clustering_zdistrpar_zmin+j*da;
    printf("min=%le max=%le\n",tomo.clustering_zmin[k],tomo.clustering_zmax[k]);
  }
  printf("min=%le max=%le\n",tomo.clustering_zmin[tomo.clustering_Nbin-1],tomo.clustering_zmax[tomo.clustering_Nbin-1]);
  printf("redshift.clustering_zdistrpar_zmin=%le max=%le\n",redshift.clustering_zdistrpar_zmin,redshift.clustering_zdistrpar_zmax);
  free_double_vector(sum,0,zbins);
    gbias.b1_function = &b1_per_bin;
  for (i =0; i < tomo.clustering_Nbin ; i++){
    // gbias.b[i] = 0.95/(growfac(1./(1.+(tomo.clustering_zmax[i]+tomo.clustering_zmin[i]/2.)))/growfac(1.));
    gbias.b[i] = 1.3+0.1*i;
    printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
  }
  n=0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
      printf("GGL combinations zl=%d zs=%d accept=%d\n",i,j,test_zoverlap(i,j));
    }
  }
  tomo.ggl_Npowerspectra = n;
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
}

void set_lens_galaxies_LSST()
{
  int i,j,n,k;
  double frac, zi;
  redshift.clustering_zdistrpar_zmin = 0.2;
  redshift.clustering_zdistrpar_zmax = 1.2;
  tomo.clustering_Npowerspectra=tomo.clustering_Nbin;
  tomo.clustering_zmin[0] = redshift.clustering_zdistrpar_zmin;
  tomo.clustering_zmax[tomo.clustering_Nbin-1] = redshift.clustering_zdistrpar_zmax;

  int zbins =2000;
  double da = (redshift.clustering_zdistrpar_zmax-redshift.clustering_zdistrpar_zmin)/(1.0*zbins);
  double *sum;
  sum=create_double_vector(0, zbins);
  
  sum[0] = 0.0;
  for (k = 0, zi = redshift.clustering_zdistrpar_zmin; k<zbins; k++,zi+=da){
    sum[k+1] = sum[k]+pf_histo(zi, NULL);
  }
  printf("\n");
  printf("Source Sample - Tomographic Bin limits:\n");
  for(k=0;k<tomo.clustering_Nbin-1;k++){
    frac=(k+1.)/(1.*tomo.clustering_Nbin)*sum[zbins-1];
    j = 0;
    while (sum[j]< frac){
      j++;
    }
    tomo.clustering_zmax[k] = redshift.clustering_zdistrpar_zmin+j*da;
    tomo.clustering_zmin[k+1] = redshift.clustering_zdistrpar_zmin+j*da;
    printf("min=%le max=%le\n",tomo.clustering_zmin[k],tomo.clustering_zmax[k]);
  }
  printf("min=%le max=%le\n",tomo.clustering_zmin[tomo.clustering_Nbin-1],tomo.clustering_zmax[tomo.clustering_Nbin-1]);
  printf("redshift.clustering_zdistrpar_zmin=%le max=%le\n",redshift.clustering_zdistrpar_zmin,redshift.clustering_zdistrpar_zmax);
  free_double_vector(sum,0,zbins);
    gbias.b1_function = &b1_per_bin;
  for (i =0; i < tomo.clustering_Nbin ; i++){
    gbias.b[i] = 0.95/(growfac(1./(1.+(tomo.clustering_zmax[i]+tomo.clustering_zmin[i]/2.)))/growfac(1.));
    // gbias.b[i] = 1.3+0.1*i;
    printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
  }
  n=0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
      printf("GGL combinations zl=%d zs=%d accept=%d\n",i,j,test_zoverlap(i,j));
    }
  }
  tomo.ggl_Npowerspectra = n;
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
}  


void set_clusters_LSST(){
  int i,j;
  //N200->M relationship from Murata et al. (2018)
  nuisance.cluster_Mobs_lgN0 = 3.207; //fiducial: 3.207, flat prior [0.5, 5.0]
  nuisance.cluster_Mobs_alpha = 0.993; //fiducial: 0.993, flat prior [0.0, 2.0]
  nuisance.cluster_Mobs_beta = 0.0; //fiducial: 0.0, flat prior [-1.5, 1.5]
  nuisance.cluster_Mobs_sigma0 = 0.456; //fiducial: 0.456, flat prior [0.0, 1.5]
  nuisance.cluster_Mobs_sigma_qm = 0.0; //fiducial: -0.169, flat prior [-1.5, 1.5]
  nuisance.cluster_Mobs_sigma_qz = 0.0; //fiducial: 0.0, flat prior [-1.5, 1.5]
  //Compliteness parameters are not marinilized, but just fixed to 1.
  nuisance.cluster_completeness[0] = 1.0;
  nuisance.cluster_completeness[1] = 1.0;
  nuisance.cluster_completeness[2] = 1.0;
  nuisance.cluster_completeness[3] = 1.0;

  //no miscentering so far
  nuisance.cluster_centering_f0 = 1.0;
  nuisance.cluster_centering_alpha = 0;
  nuisance.cluster_centering_sigma = 0;
  nuisance.cluster_centering_M_pivot = 1.e+14;
  printf("%e %e %e %e %e %e\n",nuisance.cluster_Mobs_lgN0, nuisance.cluster_Mobs_alpha, nuisance.cluster_Mobs_beta, nuisance.cluster_Mobs_sigma0,  nuisance.cluster_Mobs_sigma_qm, nuisance.cluster_Mobs_sigma_qz);
  tomo.cluster_Nbin = 4; // number of cluster redshift bins
  tomo.cluster_zmin[0] = 0.4;
  tomo.cluster_zmax[0] = 0.6;
  tomo.cluster_zmin[1] = 0.6;
  tomo.cluster_zmax[1] = 0.8;
  tomo.cluster_zmin[2] = 0.8;
  tomo.cluster_zmax[2] = 1.0;
  tomo.cluster_zmin[3] = 1.0;
  tomo.cluster_zmax[3] = 1.2;
  tomo.cgl_Npowerspectra = 0;// number of cluster-lensing tomography combinations
  for (i = 0; i < tomo.cluster_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      tomo.cgl_Npowerspectra += test_zoverlap_c(i,j);
    }
  }
  
  Cluster.N200_min = 25.; //formerly 20
  Cluster.N200_max = 220.;
  Cluster.N200_Nbin = 4;
  strcpy(Cluster.model,"Murata_etal_2018");
  //upper bin boundaries - note that bin boundaries need to be integers!
  int Nlist[4] = {40,80,120,Cluster.N200_max}; //formerly 40,80,120,...
  Cluster.N_min[0] = Cluster.N200_min;
  Cluster.N_max[0] = Nlist[0];
  for (i = 1; i < Cluster.N200_Nbin; i++){
    Cluster.N_min[i] = Nlist[i-1];
    Cluster.N_max[i] = Nlist[i];
  }
 for (i = 0; i < Cluster.N200_Nbin; i++){
    printf ("Richness bin %d: %e - %e, N(z = 0.3) = %e, N(z = 0.7) = %e\n", i,Cluster.N_min[i],Cluster.N_max[i],N_N200(0,i),N_N200(2,i));
  }
  printf("Clusters set to LSST Y10\n");
  printf("Clusters cgl_Npowerspectra=%d\n",tomo.cgl_Npowerspectra);
}


void init_IA(char *model,char *lumfct)
{  
  if(strcmp(lumfct,"GAMA")==0) set_LF_GAMA();
  else if(strcmp(lumfct,"DEEP2")==0) set_LF_DEEP2();
  else {
    printf("init.c:init_IA: %s lumfct not defined\n",lumfct);
    printf("USING GAMA LF INSTEAD\n");
    set_LF_GAMA();
  }
  printf("SET LUMINOSITY FUNCTION=%s\n",lumfct);
  
  nuisance.oneplusz0_ia=1.3; 
  //z0=0.3 is arbitrary pivot redshift J11 p18
  nuisance.c1rhocrit_ia=0.0134; 
  // J11 p.8
  
  if(strcmp(model,"none")==0)  like.IA=0;
  else if(strcmp(model,"NLA_HF")==0)  like.IA=1;
  else if(strcmp(model,"lin")==0)  like.IA=2;
  else{
    printf("init.c:init_IA: %s IA model not defined\n",model);
    exit(1);
  }
  printf("SET IA MODEL=%s\n",model);
  set_ia_priors();
  log_like_f_red();
}




void init_HOD_rm(){
  set_HOD_redmagic_priors();
  like.Rmin_bias = 0.1;//use halo+HOD model down to 100 kpc/h
  redm.parameterization = 0; //Zehavi et al. 2011 HOD parameterization
  redm.cg = 1.0;
  redm.fc = 0.2;
  redm.hod[0] = 12.1;
  redm.hod[1] = 0.4;
  redm.hod[2] = 13.65;
  redm.hod[3] = 12.2;
  redm.hod[4] = 1.0;
}

