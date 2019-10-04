#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include <fftw3.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "../cosmolike_core/theory/basics.c"
#include "../cosmolike_core/theory/structs.c"
#include "../cosmolike_core/theory/parameters.c"
#include "../cosmolike_core/emu17/P_cb/emu.c"
#include "../cosmolike_core/theory/recompute.c"
#include "../cosmolike_core/theory/cosmo3D.c"
#include "../cosmolike_core/theory/redshift_spline.c"
#include "../cosmolike_core/theory/halo.c"
#include "../cosmolike_core/theory/HOD.c"
#include "../cosmolike_core/theory/pt.c"
#include "../cosmolike_core/theory/cosmo2D_fourier.c"
#include "../cosmolike_core/theory/IA.c"
#include "../cosmolike_core/theory/cluster.c"
#include "../cosmolike_core/theory/BAO.c"
#include "../cosmolike_core/theory/external_prior.c"
#include "../cosmolike_core/theory/covariances_3D.c"
#include "../cosmolike_core/theory/covariances_fourier.c"
#include "../cosmolike_core/theory/covariances_cluster.c"
#include "init_emu.c"

void run_cov_N_N (char *OUTFILE, char *PATH, int nzc1, int nzc2,int start);
void run_cov_cgl_N (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster,int N1, int nzc2, int start);
void run_cov_cgl_cgl (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int start);
void run_cov_cgl_cgl_all (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster);
void run_cov_shear_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start);
void run_cov_shear_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start);
void run_cov_ggl_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start);
void run_cov_ggl_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start);
void run_cov_cl_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start);
void run_cov_cl_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start);

void run_cov_ggl_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_clustering_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_clustering_ggl(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_clustering(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_ggl(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_shear_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);


void run_cov_N_N (char *OUTFILE, char *PATH, int nzc1, int nzc2,int start)
{
  int nN1, nN2,i,j;
  double cov;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  for (nN1 = 0; nN1 < Cluster.N200_Nbin; nN1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      i += Cluster.N200_Nbin*nzc1+nN1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;

      cov =cov_N_N(nzc1,nN1, nzc2, nN2);
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,0.0,0.0, nzc1, nN1, nzc2, nN2,cov,0.0);
    }
  }
  fclose(F1);
}

void run_cov_cgl_N (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster,int N1, int nzc2, int start)
{
  int nN1, nN2, nl1, nzc1, nzs1,i,j;
  double cov;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzc1 = ZC(N1);
  nzs1 = ZSC(N1);
  for (nN1 = 0; nN1 < Cluster.N200_Nbin; nN1 ++){
    for( nl1 = 0; nl1 < Cluster.lbin; nl1 ++){
     for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
       i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
       i += (N1*Cluster.N200_Nbin+nN1)*Cluster.lbin +nl1;
       j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
       j += Cluster.N200_Nbin*nzc2+nN2;

       cov =cov_cgl_N(ell_Cluster[nl1],nzc1,nN1, nzs1, nzc2, nN2);
       fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell_Cluster[nl1], 0., nzc1, nzs1, nzc2, nN2,cov,0.);
     }
   }
 }
 fclose(F1);
}

void run_cov_cgl_cgl (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int start)
{
  int nN1, nN2, nl1, nzc1, nzs1, nl2, nzc2, nzs2,i,j;
  double c_g, c_ng;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzc1 = ZC(N1);
  nzs1 = ZSC(N1);
  nzc2 = ZC(N2);
  nzs2 = ZSC(N2);
  for (nN1 = 0; nN1 < Cluster.N200_Nbin; nN1 ++){
    for( nl1 = 0; nl1 < Cluster.lbin; nl1 ++){
      for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
        for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
          i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
          i += (N1*Cluster.N200_Nbin+nN1)*Cluster.lbin +nl1;
          j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
          j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;

          c_g = 0;
          c_ng = cov_NG_cgl_cgl(ell_Cluster[nl1],ell_Cluster[nl2],nzc1,nN1, nzs1, nzc2, nN2,nzs2);
          if (nl2 == nl1){c_g =cov_G_cgl_cgl(ell_Cluster[nl1],dell_Cluster[nl1],nzc1,nN1, nzs1, nzc2, nN2,nzs2);}
          fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell_Cluster[nl1],ell_Cluster[nl2], nzc1, nzs1, nzc2, nzs2,c_g, c_ng);
        }
      }
    }
  }
  fclose(F1);
}

void run_cov_cgl_cgl_all (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster)
{
  int nN1, nN2, nl1, nzc1, nzs1, nl2, nzc2, nzs2,i,j, N1,N2;
  double c_g, c_ng;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s",PATH,OUTFILE);
  F1 =fopen(filename,"w");
  for (N1 = 0; N1 < tomo.cgl_Npowerspectra; N1 ++){
    for (N2 = 0; N2 < tomo.cgl_Npowerspectra; N2 ++){
      nzc1 = ZC(N1);
      nzs1 = ZSC(N1);
      nzc2 = ZC(N2);
      nzs2 = ZSC(N2);
      for (nN1 = 0; nN1 < Cluster.N200_Nbin; nN1 ++){
        for( nl1 = 0; nl1 < Cluster.lbin; nl1 ++){
          for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
            for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
              i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
              i += (N1*Cluster.N200_Nbin+nN1)*Cluster.lbin +nl1;
              j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
              j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;         
              c_g = 0;
              c_ng = cov_NG_cgl_cgl(ell_Cluster[nl1],ell_Cluster[nl2],nzc1,nN1, nzs1, nzc2, nN2,nzs2);
              if (nl2 == nl1){c_g =cov_G_cgl_cgl(ell_Cluster[nl1],dell_Cluster[nl1],nzc1,nN1, nzs1, nzc2, nN2,nzs2);}
              fprintf(F1,"%d %d %e %e %d %d %d  %d %d %d  %e %e\n",i,j,ell_Cluster[nl1],ell_Cluster[nl2], nzc1, nN1,nzs1, nzc2, nN2,nzs2,c_g, c_ng);
            }
          }
        }
      }
    }
  }
  fclose(F1);
}

void run_cov_shear_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start)
{
  int nz1,nz2, nN2, nl1, nzc1, i,j;
  double cov;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nz1 = Z1(N1);
  nz2 = Z2(N1);
  for( nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      cov = 0.;
      i = like.Ncl*N1+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;

      if (ell[nl1] < like.lmax_shear){cov =cov_shear_N(ell[nl1],nz1,nz2, nzc2, nN2);}
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell[nl1], 0., nz1, nz2, nzc2, nN2,cov,0.);
    }
  }
  fclose(F1);
}

void run_cov_shear_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start)
{
  int nN1, nN2, nzs1, nzs2,nl2, nzc2, nzs3,i,j;
  double c_g, c_ng;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzs1 = Z1(N1);
  nzs2 = Z2(N1);
  nzc2 = ZC(N2);
  nzs3 = ZSC(N2);
  for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
        i = like.Ncl*N1+nl1;
        j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
        j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;          
        c_g = 0.;
        c_ng = 0.;
        if (ell[nl1] < like.lmax_shear){
          c_ng = cov_NG_shear_cgl(ell[nl1],ell_Cluster[nl2],nzs1, nzs2, nzc2, nN2,nzs3);
          if (fabs(ell[nl1]/ell_Cluster[nl2] -1.) < 0.001){ 
            c_g =cov_G_shear_cgl(ell[nl1],dell_Cluster[nl2],nzs1,nzs2, nzc2, nN2,nzs3);
          }
        }
        fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], nzs1, nzs2, nzc2, nzs3,c_g, c_ng);
      }
    }
  }
  fclose(F1);
}

void run_cov_ggl_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start)
{
  int zl,zs, nN2, nl1, nzc1, i,j;
  double cov,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  zl = ZL(N1);
  zs = ZS(N1);
  for( nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      i = like.Ncl*(tomo.shear_Npowerspectra+N1)+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;
      cov = 0.;
      weight = test_kmax(ell[nl1],zl);
      if (weight){
        cov =cov_ggl_N(ell[nl1],zl,zs, nzc2, nN2);
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell[nl1], 0., zl, zs, nzc2, nN2,cov,0.);
    }
  }
  fclose(F1);
}

void run_cov_ggl_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start)
{
  int nN2, zl, zs, nzs1, nl2, nzc2, nzs3,i,j;
  double c_g, c_ng,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  zl = ZL(N1);
  zs = ZS(N1);
  nzc2 = ZC(N2);
  nzs3 = ZSC(N2);
  for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
        i = like.Ncl*(tomo.shear_Npowerspectra+N1)+nl1;
        j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
        j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;

        c_g = 0; c_ng = 0.;
        weight = test_kmax(ell[nl1],zl);
        if (weight){
          c_ng = cov_NG_ggl_cgl(ell[nl1],ell_Cluster[nl2],zl,zs, nzc2, nN2,nzs3);
          if (fabs(ell[nl1]/ell_Cluster[nl2] -1.) < 0.1){c_g =cov_G_ggl_cgl(ell[nl1],dell_Cluster[nl2],zl,zs, nzc2, nN2,nzs3);}
        }
        fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], zl, zs, nzc2, nzs3,c_g, c_ng);
      }
    }
  }
  fclose(F1);
}

void run_cov_cl_N (char *OUTFILE, char *PATH, double *ell, double *dell,int N1, int nzc2, int start)
{
  int zl,zs, nN2, nl1, nzc1, i,j;
  double cov,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  for( nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+N1)+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;
      cov = 0.;
      weight = test_kmax(ell[nl1],N1);
      if (weight){
        cov =cov_cl_N(ell[nl1],N1,N1,nzc2,nN2);
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell[nl1], 0., N1, N1, nzc2, nN2,cov,0.);
    }
  }
  fclose(F1);
}

void run_cov_cl_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start)
{
  int nN2,nzc2, nzs3,i,j,nl2;
  double c_g, c_ng,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzc2 = ZC(N2);
  nzs3 = ZSC(N2);
  for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
        i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+N1)+nl1;
        j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
        j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;

        c_g = 0; c_ng = 0.;
        weight = test_kmax(ell[nl1],N1);
        if (weight){
          c_ng = cov_NG_cl_cgl(ell[nl1],ell_Cluster[nl2],N1,N1, nzc2, nN2,nzs3);
          if (fabs(ell[nl1]/ell_Cluster[nl2] -1.) < 0.1){
            c_g =cov_G_cl_cgl(ell[nl1],dell_Cluster[nl2],N1,N1, nzc2, nN2,nzs3);
          }
        }
        fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], N1,N1, nzc2, nzs3,c_g, c_ng);
        //printf("%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], N1,N1, nzc2, nzs3,c_g, c_ng);
      }
    }
  }
  fclose(F1);
}

void run_cov_ggl_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int zl,zs,z3,z4,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w"); 
  zl = ZL(n1); zs = ZS(n1);
  printf("\nN_ggl = %d (%d, %d)\n", n1,zl,zs);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 <like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],zl);
      if (weight && ell[nl2] < like.lmax_shear){
        if (test_zoverlap(zl,z3)*test_zoverlap(zl,z4)){
          c_ng = cov_NG_gl_shear_tomo(ell[nl1],ell[nl2],zl,zs,z3,z4);
        }
        if (nl1 == nl2){
          c_g =  cov_G_gl_shear_tomo(ell[nl1],dell[nl1],zl,zs,z3,z4);
        }
      }
      fprintf(F1, "%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+n1)+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],zl,zs,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_clustering_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl = %d \n", n1);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 <like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],z1);
      if (weight && ell[nl2] < like.lmax_shear){
        if (test_zoverlap(z1,z3)*test_zoverlap(z1,z4)){
          c_ng = cov_NG_cl_shear_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);
        }
        if (nl1 == nl2){
          c_g =  cov_G_cl_shear_tomo(ell[nl1],dell[nl1],z1,z2,z3,z4);
        }
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_clustering_ggl(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,zl,zs,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl_1 = %d \n", n1);
  zl = ZL(n2); zs = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl,zs);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      if (z1 == zl){
        weight = test_kmax(ell[nl1],z1)*test_kmax(ell[nl2],zl);
        if (weight){
          c_ng = cov_NG_cl_gl_tomo(ell[nl1],ell[nl2],z1,z2,zl,zs);
          if (nl1 == nl2){
            c_g =  cov_G_cl_gl_tomo(ell[nl1],dell[nl1],z1,z2,zl,zs);
          }
        }
      }
      fprintf(F1, "%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,like.Ncl*(tomo.shear_Npowerspectra+n2)+nl2, ell[nl1],ell[nl2],z1,z2,zl,zs,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_clustering(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl_1 = %d \n", n1);
  z3 = n2; z4 = n2;
  printf("N_cl_2 = %d\n", n2);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      if (z1 == z3){
        weight = test_kmax(ell[nl1],z1)*test_kmax(ell[nl2],z3);
        if (weight) {
          c_ng = cov_NG_cl_cl_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);
        }
        if (nl1 == nl2){
          c_g =  cov_G_cl_cl_tomo(ell[nl1],dell[nl1],z1,z2,z3,z4);
        }
      }
      fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra + n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_ggl(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int zl1,zl2,zs1,zs2,nl1,nl2, weight;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");

  zl1 = ZL(n1); zs1 = ZS(n1);
  printf("\nN_tomo_1 = %d (%d, %d)\n", n1,zl1,zs1);
  zl2 = ZL(n2); zs2 = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl2,zs2);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],zl1)*test_kmax(ell[nl2],zl2);
      if (weight && zl1 == zl2) {
        c_ng = cov_NG_gl_gl_tomo(ell[nl1],ell[nl2],zl1,zs1,zl2,zs2);
      }
      if (nl1 == nl2){
        c_g =  cov_G_gl_gl_tomo(ell[nl1],dell[nl1],zl1,zs1,zl2,zs2);
      }
      if (weight ==0 && n2 != n1){
        c_g = 0;
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+n1)+nl1,like.Ncl*(tomo.shear_Npowerspectra+n2)+nl2, ell[nl1],ell[nl2],zl1,zs1,zl2,zs2,c_g,c_ng);   
    }
  }
  fclose(F1);
}


void run_cov_shear_shear(char *OUTFILE, char *PATH, double *ell, double *dell,int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  z1 = Z1(n1); z2 = Z2(n1);
  printf("N_shear = %d\n", n1);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n",n2,z3,z4);
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      if (ell[nl1] < like.lmax_shear && ell[nl2] < like.lmax_shear){
        c_ng = cov_NG_shear_shear_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);
      }
      if (nl1 == nl2){
        c_g =  cov_G_shear_shear_tomo(ell[nl1],dell[nl1],z1,z2,z3,z4);
        if (ell[nl1] > like.lmax_shear && n1!=n2){c_g = 0.;} 
      }         
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",like.Ncl*n1+nl1,like.Ncl*(n2)+nl2,ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
      //printf("%d %d %e %e %d %d %d %d %e %e\n", like.Ncl*n1+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}


int main(int argc, char** argv)
{
  
  int i,l,m,n,o,s,p,nl1,t,k;
  char OUTFILE[400],filename[400],arg1[400],arg2[400];
  
  int N_scenarios=36;
  double area_table[36]={7623.22,14786.3,9931.47,8585.43,17681.8,15126.9,9747.99,8335.08,9533.42,18331.3,12867.8,17418.9,19783.1,12538.8,15260.0,16540.7,19636.8,11112.7,10385.5,16140.2,18920.1,17976.2,11352.0,9214.77,16910.7,11995.6,16199.8,14395.1,8133.86,13510.5,19122.3,15684.5,12014.8,14059.7,10919.3,13212.7};
  double nsource_table[36]={13.991,33.3975,17.069,28.4875,35.3643,10.3802,11.1105,29.5904,31.4608,15.7463,13.3066,9.07218,9.58719,16.3234,25.8327,10.7415,38.0412,32.8821,19.8893,27.0369,15.1711,14.2418,19.1851,26.926,22.0012,12.6553,18.6304,11.9787,36.8728,22.5265,17.4381,12.3424,10.0095,23.5979,20.8771,24.8956};
  
  double nlens_table[36]={22.9726 ,61.5948 ,28.7809 ,51.4354 ,65.7226 ,16.3777 ,17.6899 ,53.6984 ,57.562 ,26.266 ,21.703 ,14.0587 ,14.9667 ,27.36 ,46.0364 ,17.0253 ,71.39 ,60.5184 ,34.2283 ,48.4767 ,25.1812 ,23.44 ,32.8578 ,48.2513 ,38.3766 ,20.5029 ,31.783 ,19.2647 ,68.9094 ,39.4168 ,29.4874 ,19.9292 ,15.7162 ,41.5487 ,36.1616 ,44.148};
  char survey_designation[1][200]={"LSST_Y10"};
  
  char source_zfile[36][400]={"wl_redshift_model0_WLz01.880307e-01_WLalpha8.485694e-01.txt", "wl_redshift_model1_WLz01.731166e-01_WLalpha7.662434e-01.txt", "wl_redshift_model2_WLz01.846221e-01_WLalpha8.297540e-01.txt", "wl_redshift_model3_WLz01.758423e-01_WLalpha7.812893e-01.txt", "wl_redshift_model4_WLz01.721357e-01_WLalpha7.608290e-01.txt", "wl_redshift_model5_WLz01.931476e-01_WLalpha8.768147e-01.txt", "wl_redshift_model6_WLz01.919821e-01_WLalpha8.703814e-01.txt", "wl_redshift_model7_WLz01.751912e-01_WLalpha7.776952e-01.txt", "wl_redshift_model8_WLz01.741405e-01_WLalpha7.718958e-01.txt", "wl_redshift_model9_WLz01.860048e-01_WLalpha8.373865e-01.txt", "wl_redshift_model10_WLz01.888904e-01_WLalpha8.533151e-01.txt", "wl_redshift_model11_WLz01.954564e-01_WLalpha8.895592e-01.txt", "wl_redshift_model12_WLz01.945099e-01_WLalpha8.843347e-01.txt", "wl_redshift_model13_WLz01.853877e-01_WLalpha8.339802e-01.txt", "wl_redshift_model14_WLz01.775192e-01_WLalpha7.905458e-01.txt", "wl_redshift_model15_WLz01.925611e-01_WLalpha8.735775e-01.txt", "wl_redshift_model16_WLz01.708849e-01_WLalpha7.539247e-01.txt", "wl_redshift_model17_WLz01.733832e-01_WLalpha7.677150e-01.txt", "wl_redshift_model18_WLz01.820009e-01_WLalpha8.152851e-01.txt", "wl_redshift_model19_WLz01.767381e-01_WLalpha7.862345e-01.txt", "wl_redshift_model20_WLz01.866426e-01_WLalpha8.409073e-01.txt", "wl_redshift_model21_WLz01.877261e-01_WLalpha8.468882e-01.txt", "wl_redshift_model22_WLz01.826188e-01_WLalpha8.186960e-01.txt", "wl_redshift_model23_WLz01.768086e-01_WLalpha7.866234e-01.txt", "wl_redshift_model24_WLz01.802711e-01_WLalpha8.057362e-01.txt", "wl_redshift_model25_WLz01.897506e-01_WLalpha8.580632e-01.txt", "wl_redshift_model26_WLz01.831217e-01_WLalpha8.214720e-01.txt", "wl_redshift_model27_WLz01.906926e-01_WLalpha8.632630e-01.txt", "wl_redshift_model28_WLz01.714197e-01_WLalpha7.568766e-01.txt", "wl_redshift_model29_WLz01.798667e-01_WLalpha8.035040e-01.txt", "wl_redshift_model30_WLz01.842554e-01_WLalpha8.277299e-01.txt", "wl_redshift_model31_WLz01.901797e-01_WLalpha8.604322e-01.txt", "wl_redshift_model32_WLz01.937710e-01_WLalpha8.802561e-01.txt", "wl_redshift_model33_WLz01.790701e-01_WLalpha7.991072e-01.txt", "wl_redshift_model34_WLz01.811701e-01_WLalpha8.106987e-01.txt", "wl_redshift_model35_WLz01.781525e-01_WLalpha7.940419e-01.txt"};

  char lens_zfile[36][400]={"LSS_redshift_model0_LSSz02.629496e-01_LSSalpha9.285983e-01.txt","LSS_redshift_model1_LSSz02.852923e-01_LSSalpha8.985943e-01.txt","LSS_redshift_model2_LSSz02.664822e-01_LSSalpha9.186036e-01.txt","LSS_redshift_model3_LSSz02.798758e-01_LSSalpha9.014201e-01.txt","LSS_redshift_model4_LSSz02.873873e-01_LSSalpha8.978683e-01.txt","LSS_redshift_model5_LSSz02.593970e-01_LSSalpha9.470921e-01.txt","LSS_redshift_model6_LSSz02.600213e-01_LSSalpha9.425114e-01.txt","LSS_redshift_model7_LSSz02.811155e-01_LSSalpha9.006370e-01.txt","LSS_redshift_model8_LSSz02.831875e-01_LSSalpha8.995165e-01.txt","LSS_redshift_model9_LSSz02.649368e-01_LSSalpha9.224338e-01.txt","LSS_redshift_model10_LSSz02.622058e-01_LSSalpha9.314128e-01.txt","LSS_redshift_model11_LSSz02.584820e-01_LSSalpha9.568082e-01.txt","LSS_redshift_model12_LSSz02.588053e-01_LSSalpha9.527220e-01.txt","LSS_redshift_model13_LSSz02.656075e-01_LSSalpha9.206866e-01.txt","LSS_redshift_model14_LSSz02.768397e-01_LSSalpha9.037492e-01.txt","LSS_redshift_model15_LSSz02.596975e-01_LSSalpha9.447600e-01.txt","LSS_redshift_model16_LSSz02.901709e-01_LSSalpha8.971658e-01.txt","LSS_redshift_model17_LSSz02.847362e-01_LSSalpha8.988183e-01.txt","LSS_redshift_model18_LSSz02.698330e-01_LSSalpha9.121821e-01.txt","LSS_redshift_model19_LSSz02.782257e-01_LSSalpha9.026084e-01.txt","LSS_redshift_model20_LSSz02.642756e-01_LSSalpha9.243038e-01.txt","LSS_redshift_model21_LSSz02.632273e-01_LSSalpha9.276296e-01.txt","LSS_redshift_model22_LSSz02.689934e-01_LSSalpha9.135968e-01.txt","LSS_redshift_model23_LSSz02.780987e-01_LSSalpha9.027073e-01.txt","LSS_redshift_model24_LSSz02.723464e-01_LSSalpha9.085463e-01.txt","LSS_redshift_model25_LSSz02.615210e-01_LSSalpha9.343470e-01.txt","LSS_redshift_model26_LSSz02.683327e-01_LSSalpha9.147934e-01.txt","LSS_redshift_model27_LSSz02.608392e-01_LSSalpha9.376962e-01.txt","LSS_redshift_model28_LSSz02.889655e-01_LSSalpha8.974355e-01.txt","LSS_redshift_model29_LSSz02.729686e-01_LSSalpha9.077654e-01.txt","LSS_redshift_model30_LSSz02.669178e-01_LSSalpha9.176391e-01.txt","LSS_redshift_model31_LSSz02.612016e-01_LSSalpha9.358553e-01.txt","LSS_redshift_model32_LSSz02.591077e-01_LSSalpha9.496317e-01.txt","LSS_redshift_model33_LSSz02.742326e-01_LSSalpha9.063038e-01.txt","LSS_redshift_model34_LSSz02.710103e-01_LSSalpha9.103760e-01.txt","LSS_redshift_model35_LSSz02.757517e-01_LSSalpha9.047459e-01.txt"};


    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // use this remapping only if files fail !!!!!!!!! 
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // int fail[5]={1134, 259, 497, 623, 718};
    // hit=fail[hit-1];
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
  int hit=atoi(argv[1]);
  Ntable.N_a=20;
  k=1;
   
  for(t=8;t<9;t++){
   
    //RUN MODE setup
    init_cosmo_runmode("emu");
    init_binning_fourier(15,20.0,3000.0,3000.0,21.0,10,10);
    init_survey("LSST");
    sprintf(arg1,"zdistris/%s",source_zfile[t]);
    sprintf(arg2,"zdistris/%s",lens_zfile[t]); 
    init_galaxies(arg1,arg2,"none","none","LSST_Y10");
    init_clusters();
    init_IA("none", "GAMA"); 
    init_probes("3x2pt");
   

    //set l-bins for shear, ggl, clustering, clusterWL
    double logdl=(log(like.lmax)-log(like.lmin))/like.Ncl;
    double *ell, *dell, *ell_Cluster, *dell_Cluster;
    ell=create_double_vector(0,like.Ncl-1);
    dell=create_double_vector(0,like.Ncl-1);
    ell_Cluster=create_double_vector(0,Cluster.lbin-1);
    dell_Cluster=create_double_vector(0,Cluster.lbin-1);
    int j=0;
    for(i=0;i<like.Ncl;i++){
      ell[i]=exp(log(like.lmin)+(i+0.5)*logdl);
      dell[i]=exp(log(like.lmin)+(i+1)*logdl)-exp(log(like.lmin)+(i*logdl));
      if(ell[i]<like.lmax_shear) printf("%le\n",ell[i]);
      if(ell[i]>like.lmax_shear){
        ell_Cluster[j]=ell[i];
        dell_Cluster[j]=dell[i];
        printf("%le %le\n",ell[i],ell_Cluster[j]);
        j++;
      }
    } 


    printf("----------------------------------\n");  
    survey.area=area_table[t];
    survey.n_gal=nsource_table[t];
    survey.n_lens=nlens_table[t]; 
    
    sprintf(survey.name,"%s_area%le_ng%le_nl%le",survey_designation[0],survey.area,survey.n_gal,survey.n_lens);
    printf("area: %le n_source: %le n_lens: %le\n",survey.area,survey.n_gal,survey.n_lens);

//    sprintf(covparams.outdir,"/home/u17/timeifler/covparallel/"); 
    sprintf(covparams.outdir,"/aurora_nobackup/cosmos/teifler/covparallel/");

    printf("----------------------------------\n");  
    sprintf(OUTFILE,"%s_ssss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Npowerspectra; l++){
      for (m=l;m<tomo.shear_Npowerspectra; m++){
        if(k==hit){ 
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {
            run_cov_shear_shear(OUTFILE,covparams.outdir,ell,dell,l,m,k);
          }
        }
        k=k+1;
        //printf("%d\n",k);
      }
    }

    sprintf(OUTFILE,"%s_lsls_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.ggl_Npowerspectra; l++){
      for (m=l;m<tomo.ggl_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {
            run_cov_ggl(OUTFILE,covparams.outdir,ell,dell,l,m,k);
          }
        }
       //printf("%d\n",k);
        k=k+1;
      }
    }
    sprintf(OUTFILE,"%s_llll_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Npowerspectra; l++){ //auto bins only for now!
      for (m=l;m<tomo.clustering_Npowerspectra; m++){
        if(k==hit){ 
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {
            run_cov_clustering(OUTFILE,covparams.outdir,ell,dell,l,m,k);
          }
        }
        k=k+1;
        //printf("%d %d %d\n",l,m,k);
      }
    }
    sprintf(OUTFILE,"%s_llss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Npowerspectra; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {
            run_cov_clustering_shear(OUTFILE,covparams.outdir,ell,dell,l,m,k);
          }
        }
        k=k+1;
        //printf("%d\n",k);
      }
    }
    sprintf(OUTFILE,"%s_llls_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Npowerspectra; l++){
      for (m=0;m<tomo.ggl_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {
            run_cov_clustering_ggl(OUTFILE,covparams.outdir,ell,dell,l,m,k);
          }
        }
        k=k+1;
        //printf("%d\n",k);
      }
    }
    sprintf(OUTFILE,"%s_lsss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.ggl_Npowerspectra; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {
            run_cov_ggl_shear(OUTFILE,covparams.outdir,ell,dell,l,m,k);
          }
        }
        k=k+1;
        //printf("%d\n",k);
      }
    }

    //****************************** 
    //******cluster covariance****** 
    //******************************
    
   //  sprintf(OUTFILE,"%s_nn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   //  for (l=0;l<tomo.cluster_Nbin; l++){
   //    for (m=0;m<tomo.cluster_Nbin; m++){
   //      if(k==hit){
   //        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //        if (fopen(filename, "r") != NULL){exit(1);}
   //        else {
   //          run_cov_N_N (OUTFILE,covparams.outdir,l,m,k);
   //        }
   //      }
   //      k=k+1;
   //      //printf("%d\n",k);
   //    }
   //  }
   //  sprintf(OUTFILE,"%s_cscs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   //  for (l=0;l<tomo.cgl_Npowerspectra; l++){
   //    for (m=0;m<tomo.cgl_Npowerspectra; m++){
   //      if(k==hit){
   //        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //        if (fopen(filename, "r") != NULL){exit(1);}
   //        else {
   //          run_cov_cgl_cgl (OUTFILE,covparams.outdir,ell_Cluster,dell_Cluster,l,m,k);
   //        }
   //      }
   //      k=k+1;
   //    } 
   //  }
   //  sprintf(OUTFILE,"%s_csn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   //  for (l=0;l<tomo.cgl_Npowerspectra; l++){
   //    for (m=0;m<tomo.cluster_Nbin; m++){
   //      if(k==hit){
   //        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //        if (fopen(filename, "r") != NULL){exit(1);}
   //        else {
   //          run_cov_cgl_N (OUTFILE,covparams.outdir,ell_Cluster,dell_Cluster,l,m,k);
   //        }
   //      }
   //      k=k+1;
   //    }
   //  }
   //  //shear X cluster
   //  sprintf(OUTFILE,"%s_ssn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   //  for (l=0;l<tomo.shear_Npowerspectra; l++){
   //    for (m=0;m<tomo.cluster_Nbin; m++){
   //      if(k==hit){
   //        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //        if (fopen(filename, "r") != NULL){exit(1);}
   //        else {
   //         run_cov_shear_N (OUTFILE,covparams.outdir,ell,dell,l,m,k);
   //       }
   //     }
   //     k=k+1;
   //   }
   // } 
   // sprintf(OUTFILE,"%s_sscs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   // for (l=0;l<tomo.shear_Npowerspectra; l++){
   //    for (m=0;m<tomo.cgl_Npowerspectra; m++){
   //      //for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
   //        if(k==hit){
   //          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //          if (fopen(filename, "r") != NULL){exit(1);}
   //          else {
   //            run_cov_shear_cgl (OUTFILE,covparams.outdir,ell,dell,ell_Cluster,dell_Cluster,l,m,nl1,k);
   //          }
   //        }
   //        k=k+1;
   //      //}
   //    }
   //  }
   //  // ggl X cluster
   //  sprintf(OUTFILE,"%s_lsn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   //  for (l=0;l<tomo.ggl_Npowerspectra; l++){
   //    for (m=0;m<tomo.cluster_Nbin; m++){
   //      if(k==hit){
   //        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //        if (fopen(filename, "r") != NULL){exit(1);}
   //        else {
   //         run_cov_ggl_N (OUTFILE,covparams.outdir,ell,dell,l,m,k);
   //        }
   //      }
   //      k=k+1;
   //    }
   //  }
   //  sprintf(OUTFILE,"%s_lscs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   //  for (l=0;l<tomo.ggl_Npowerspectra; l++){
   //    for (m=0;m<tomo.cgl_Npowerspectra; m++){
   //      //for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
   //        if(k==hit){
   //          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //          if (fopen(filename, "r") != NULL){exit(1);}
   //          else {
   //            run_cov_ggl_cgl (OUTFILE,covparams.outdir,ell,dell,ell_Cluster,dell_Cluster,l,m,nl1,k);
   //          }
   //        }       
   //        k=k+1;
   //      //}
   //    }
   //  }
   //  // clustering X cluster
   //  sprintf(OUTFILE,"%s_lln_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   //  for (l=0;l<tomo.clustering_Npowerspectra; l++){
   //    for (m=0;m<tomo.cluster_Nbin; m++){
   //      if(k==hit){
   //        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //        if (fopen(filename, "r") != NULL){exit(1);}
   //        else {
   //         run_cov_cl_N (OUTFILE,covparams.outdir,ell,dell,l,m,k);          
   //        }
   //      } 
   //      //printf("%d\n",k);
   //      k=k+1;
   //    }
   //  }
   //  sprintf(OUTFILE,"%s_llcs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   //  for (l=0;l<tomo.clustering_Npowerspectra; l++){
   //    for (m=0;m<tomo.cgl_Npowerspectra; m++){
   //      //for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
   //        if(k==hit){
   //          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //          if (fopen(filename, "r") != NULL){exit(1);}
   //          else {
   //            run_cov_cl_cgl (OUTFILE,covparams.outdir,ell,dell,ell_Cluster,dell_Cluster,l,m,nl1,k);
   //          }
   //        }
   //        k=k+1;
   //      //}
   //    }
   //  }
  }
  printf("number of cov blocks for parallelization: %d\n",k-1); 
  printf("-----------------\n");
  printf("PROGRAM EXECUTED\n");
  printf("-----------------\n");
  return 0;   
}

