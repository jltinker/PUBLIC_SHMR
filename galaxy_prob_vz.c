#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"


/* These are the parameters for creating the pdf table.
 */
int nv1,nr1,nphi1;
void create_pdf_table();
double ***ppp,*vvv,*rbin,*phibin,***yy2,FIRSTFLAG=1;
double BINSIZE;

double galaxy_prob_vz(double vel, double rad, double theta)
{
  static int flag=0,prev_cosmo=0;
  static double *y,*ya,*yb,*y2b,*xb,rv[21],rmin,rmax,dlogr;

  int i,j,k,irad,iphi;
  double vprob,m1,m2;

  if(!flag || RESET_PVZ || RESET_COSMOLOGY!=prev_cosmo)
    {
      rmin=1.7*pow(3*HOD.M_low/(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
      rmax=XI_MAX_RADIUS;
      dlogr=log(rmax/rmin)/(nr1-1);

      if((RESET_PVZ || RESET_COSMOLOGY!=prev_cosmo) && OUTPUT)
	{
	  printf("RESET: creating new table for:\n");
	  printf(" > OMEGA_M = %f\n",OMEGA_M);
	  printf(" > OMEGA_Z = %f\n",OMEGA_Z);
	  printf(" > SIGMA_8 = %f\n",SIGMA_8);
	  printf(" > VBIAS   = %f\n",VBIAS);
	  printf(" > VBIAS_C = %f\n",VBIAS_C);
	  fflush(stdout);
	}
      RESET_PVZ=0;
      prev_cosmo=RESET_COSMOLOGY;

      create_pdf_table();
      
      if(!flag)
	{
	  y=dvector(1,nphi1);
	  ya=dvector(1,nphi1);
	  yb=dvector(1,4);
	  y2b=dvector(1,4);
	  xb=dvector(1,4);
	}

      flag=1;

      for(i=1;i<=nr1;++i)
	for(j=1;j<=nphi1;++j)
	  spline(vvv,ppp[i][j],nv1,1.0E+30,1.0E+30,yy2[i][j]);

    }
  if(fabs(vel)>MAXVEL)return(0);

  if(theta<0)theta*=-1;

  for(i=1;i<=nr1;++i)
    if(rad<rbin[i])break;
  if(i<3)i=3;
  if(i>nr1-1)i=nr1-1;
  irad=i;
  while(rbin[irad-2]==0)irad++;
  if(irad>nr1-1) 
    {
      printf("ERROR CANNOT ITERPOLATE FOR r=%f\n",rad);
      exit(0);
    }

  for(j=1,i=irad-2;i<=irad+1;++i,++j)
    {
      for(k=1;k<=nphi1;++k)
	{
	  splint(vvv,ppp[i][k],yy2[i][k],nv1,vel,&y[k]);
	}
      spline(phibin,y,nphi1,1.0E+30,1.0E+30,ya);
      splint(phibin,y,ya,nphi1,theta,&yb[j]);
      xb[j]=rbin[i];
    }
  spline(xb,yb,4,1.0E30,1.0E30,y2b);
  splint(xb,yb,y2b,4,rad,&vprob);

  if(vprob<0)vprob=0;
  return(vprob);



  irad=i;
  if(irad<2)irad=2;
  if(irad>nr1)irad=nr1;
  iphi = (int)theta/(PI/18.0)+1;
  if(iphi==nphi1)iphi=nphi1-1;


  for(j=1,i=irad-1;i<=irad;++i)
    {
      for(k=iphi;k<=iphi+1;++k)
	{
	  splint(vvv,ppp[i][k],yy2[i][k],nv1,vel,&y[j++]);
	}
    }
  m1 = (y[1]-y[2])/(phibin[iphi]-phibin[iphi+1])*(theta-phibin[iphi])+y[1];
  m2 = (y[3]-y[4])/(phibin[iphi]-phibin[iphi+1])*(theta-phibin[iphi])+y[3];
  vprob = (m2-m1)/(rbin[irad]-rbin[irad-1])*(rad-rbin[irad-1])+m1;
  if(vprob<0)vprob=0;
  return(vprob);

}

/* For consistency, the MAXVEL and BINSIZE variables will be scaled by
 * OMEGA_M^0.6.
 *
 */
void create_pdf_table()
{
#define NZMASS 19

  int TEST=0;

  FILE *fp;
  double dlogr,vbias1,vbias2;
  static int CNT=0;

  double sgal[NZMASS][NZMASS],sgal_cs[NZMASS][NZMASS],sgal_cc[NZMASS][NZMASS],
    *cf,fdat[10],s2gal[4],wgal[3],**mtemp,
    ***ptemp,*vtemp,*temp,*p1temp,*p2temp,*precv,*rvir,*ngal,**ngal2,**ngal2x,*mbin;
  double binsize,ptot,vprob,p0,p1,p2,p3,s1,s2,s3,s4,v1,v,
    mass,mlo,mhi,sgal1[NZMASS][NZMASS],hmass[NZMASS],
    sgal2[NZMASS][NZMASS],w1[NZMASS][NZMASS],w2[NZMASS][NZMASS],mass1,mass2,wsum,fac,
    w3[NZMASS][NZMASS],t0,tsum=0,t1,tdiff,frac1,frac2,frac3,frac4,x1,x2,deltar,
    haloweight,immax,mmax,rmax,rmin,exfac,xx1=0,weightsum;
  int i,j,k,n,nzmass,imass=0,irad,im1,im2,istep=1,istart=1,jstart=1,jstep=1,
    idat[10],nitems;
  float f1;
  
  CNT++;

#ifdef PARALLEL
  istep = 1;
  istart = 1;
  jstart = ThisTask + 1;
  jstep = NTask;
#endif

  BINSIZE=20*pow(OMEGA_Z/0.3,GAMMA);

  rmin=R_MIN_2HALO;
  rmin=1.1*pow(3*HOD.M_low/(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  rmax=XI_MAX_RADIUS;
  nr1=30;
  nv1=301;
  nphi1=9;
  MAXVEL=BINSIZE*floor(nv1/2);


  /* This assumes that even if the cosmology is changed, that
   * the dimensions of these arrays won't change.
   */
  if(FIRSTFLAG)
    {
      vvv=dvector(1,nv1);
      rbin=dvector(1,nr1);
      phibin=dvector(1,nphi1);
      yy2=d3tensor(1,nr1,1,nphi1,1,nv1);
      ppp=d3tensor(1,nr1,1,nphi1,1,nv1);
      FIRSTFLAG=0;
    }

  if(TEST)
    {
      HOD.pdfs=0;
      HOD.M_min = HOD.M_low = 4.4e11;
      rmin=1.1*pow(3*HOD.M_low/(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
    }

  dlogr=log(rmax/rmin)/(nr1-1);
  for(i=1;i<=nr1;++i)
    rbin[i]=rmin*exp((i-1)*dlogr);
  for(i=1;i<=nphi1;++i)
    phibin[i]=(i-0.5)/nphi1*PI/2;
  for(i=1;i<=nv1;++i)
    vvv[i]=-MAXVEL+(i-1)*BINSIZE;

  nzmass=NUM_POW2MASS_BINS-1;

  ptemp=d3tensor(0,nzmass-1,0,nzmass-1,1,nv1);
  vtemp=dvector(1,nv1);
  temp=dvector(1,nv1);  
  p1temp=dvector(1,nv1);
  p2temp=dvector(1,nv1);
  precv=dvector(1,nv1);
  rvir=dvector(0,nzmass);
  mbin=dvector(0,nzmass);


  for(i=1;i<=nr1;++i)
    for(j=1;j<=nphi1;++j)
      for(k=1;k<=nv1;++k)
	ppp[i][j][k]=0;

  binsize=BINSIZE;  
  
  fac=sqrt(4.499E-48/2.0)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;
  fac=fac*fac;

  /* w1 -> number of i_sat + j_cen pairs
   * w2 -> number of j_sat + i_cen pairs
   * w3 -> number of i_sat + j_sat pairs
   *
   * Okay, I know that all you have to do is switch the indices on
   * w1 to get w2, but just to keep things really clear, I'll have three
   * different matrices.
   */

  for(i=0;i<nzmass;++i)
    for(j=0;j<nzmass;++j)
      {
	mass1=1.5*pow(2.0,i)*HOD.M_low;
	vbias1 = (VBIAS_SLOPE*(log(mass1)-34.5)+VBIAS);
	if(vbias1<0.1)vbias1=0.1;

	
	//vbias1 = 0.4/3*log(mass1)/LOGE_10 - 1;
	//if(vbias1<0.6)vbias1=0.6;

	/* TEMP-> mass threshold
	 */
	/*
	vbias1=1;
	if(mass1<VBIAS_MASS_THRESHOLD)
	  vbias1=VBIAS;
	*/
	vbias1 *= vbias1;

	mass2=1.5*pow(2.0,j)*HOD.M_low;
	vbias2 = (VBIAS_SLOPE*(log(mass2)-34.5)+VBIAS);
	if(vbias2<0.1)vbias2=0.1;

	//vbias2 = 0.4/3*log(mass2)/LOGE_10 - 1;
	//if(vbias2<0.6)vbias2=0.6;

	/* TEMP-> mass threshold
	 */
	/*
	vbias2=1;
	if(mass2<VBIAS_MASS_THRESHOLD)
	  vbias2=VBIAS;
	*/
	vbias2 *= vbias2;
	mbin[i]=mass1;

	s1=vbias1*fac*pow(mass1,0.66666666667);
	s2=vbias2*fac*pow(mass2,0.66666666667);	
	sgal[i][j]=s1+s2;

	s1=fac*pow(mass1,0.66666666667)*VBIAS_C*VBIAS_C;
	sgal_cs[i][j]=s1+s2;
	s2=fac*pow(mass2,0.66666666667)*VBIAS_C*VBIAS_C;
	sgal_cc[i][j]=s1+s2;

	wsum=N_avg(mass1)*dndM_interp(mass1)*N_avg(mass2)*dndM_interp(mass2);
	w1[i][j]=N_sat(mass1)*dndM_interp(mass1)*N_cen(mass2)*dndM_interp(mass2);
	w2[i][j]=N_sat(mass2)*dndM_interp(mass2)*N_cen(mass1)*dndM_interp(mass1);
	w3[i][j]=N_sat(mass1)*dndM_interp(mass1)*N_sat(mass2)*dndM_interp(mass2);

	// this is for testing sat-sat 2halo term
	//w1[i][j] = w2[i][j] = 0;
	//wsum = w3[i][j];
	//w1[i][j] = w2[i][j] = w3[i][j] = 0;
	//wsum = w3[i][j];
	//w1[i][j] = w2[i][j] = 0;

	w1[i][j]/=wsum;
	w2[i][j]/=wsum;
	w3[i][j]/=wsum;
      }

  /* Calculate the average number of galaxy
   * pairs in each combination of halo mass bins.
   */
  ngal=dvector(0,nzmass-1);
  ngal2=dmatrix(0,nzmass-1,0,nzmass-1);
  ngal2x=dmatrix(0,nzmass-1,0,nzmass-1);

  for(i=0;i<nzmass;++i)
    { 
      mlo=pow(2.0,i)*HOD.M_low;
      mhi=mlo*2;
      ngal[i]=qromo(func_galaxy_density,log(mlo),log(mhi),midpnt);      
      ngal2[i][i]=ngal[i]*ngal[i]*0.5;
    }
  for(i=0;i<nzmass;++i)
    for(j=i+1;j<nzmass;++j)
      ngal2[i][j]=ngal[i]*ngal[j];

  /* Calculate the virial radius for each halo bin.
   */
  for(i=0;i<nzmass;++i) 
    rvir[i]=pow(3*pow(2.0,i)*HOD.M_low*ROOT2/
		(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
  for(i=0;i<nzmass;++i) 
    hmass[i] = 1.5*pow(2.0,i)*HOD.M_low;

  for(k=1;k<=nv1;++k)
    vtemp[k]=vvv[k];
  
  /* TEMP TEMP TEMP
   */
  for(i=0;i<nzmass;++i)
    for(j=i;j<nzmass;++j)
      //ngal2x[i][j] = N_cen(hmass[i])*N_sat(hmass[j]) + N_cen(hmass[j])*N_sat(hmass[i]);
      //ngal2x[i][j] = N_sat(hmass[i])*N_sat(hmass[j]);
      ngal2x[i][j] = N_cen(hmass[i])*N_cen(hmass[j]);
  for(i=0;i<nzmass;++i)
    ngal2x[i][i]/=2;

  /* For each r/phi combination, first calculate the pdf for all possible 
   * halo-halo pairs. Once the h-h pdf is calculates, convolve with the
   * galaxy Gaussian pdf.
   */
  j=jstart;
  if(NTask==1 && OUTPUT)
    {
      printf("TASK %d starting: %d %d\n",ThisTask,istart,jstart);
      fflush(stdout);
    }
  t0=second();
  for(i=istart;i<=nr1;)
    {
      while(j>nphi1) { i++; j=j-nphi1; }
      if(i>nr1)break;
      jstart=j;
      for(j=jstart;j<=nphi1;j+=jstep)
	{
	  t1=second();
	  tdiff=timediff(t0,t1);
	  tsum+=tdiff;
	  t0=second();
	  if(NTask==1 && OUTPUT)
	    {
	      fprintf(stdout,"P(v_z) r[%d]= %5.2f phi[%d]=%4.1f [dt= %.2f T=%7.2f] (task= %d)\n",
		      i,rbin[i],j,phibin[j]*180/PI,tdiff,tsum,ThisTask);
	      fflush(stdout);
	    }
	  for(im1=0;im1<nzmass;++im1)
	    {
	      for(im2=im1;im2<nzmass;++im2)
		{
		  /* NEW VERSION, with varying VBIAS_C
		   */
		  s1=sgal_cs[im2][im1];
		  s2=sgal_cs[im1][im2];
		  s3=sgal[im1][im2];
		  s4=sgal_cc[im1][im2];

		  s2gal[0]=s1;
		  s2gal[1]=s2;
		  s2gal[2]=s3;
		  s2gal[3]=s4;
		  wgal[0]=w1[im1][im2];
		  wgal[1]=w2[im1][im2];
		  wgal[2]=w3[im1][im2];

		  vdelta_v4(mbin[im1],mbin[im2],rbin[i],phibin[j],BINSIZE,-MAXVEL,nv1,cf,
			    temp,p1temp,ptemp[im1][im2],wgal,s2gal);
		  
		}
	    }
	  

	  /* Output some stuff for testing purposes
	   */
	  if(TEST && (j==1 || j==9))
	    //if(j==1 || j==9)
	    {
	      for(im1=0;im1<nzmass;++im1)
		{
		  p0=p1=p2=0;
		  for(k=1;k<=nv1;++k)
		    {
		      v = -MAXVEL + (k-1)*BINSIZE;
		      p0 += ptemp[im1][im1][k]*BINSIZE;
		      p1 += ptemp[im1][im1][k]*v*BINSIZE;
		      p2 += ptemp[im1][im1][k]*v*v*BINSIZE;
		      //printf("xTEST%d MASS%d %f %e\n",j,im1+1,ptemp[im1][im1][k]);
		    }
		  p1 = p1/p0;
		  p2 = sqrt(p2/p0 - p1*p1);
		  //if(!isnan(p0+p1+p2))
		  //printf("TEST%d MASS%d %e %f %f %f\n",j,im1+1,hmass[im1],rbin[i],-p1,p2);
		}
	    }

	  /* Now calculate the galaxy-number-wieghted average
	   * of the v_z pdf by summing over all halo-halo pairs.
	   * (Using sum of virial radii for halo exclusion.)
	   */
	  ptot=0;
	  weightsum=0;
	  for(k=1;k<=nv1;++k)
	    for(ppp[i][j][k]=0,im1=0;im1<nzmass;++im1)
	      for(im2=im1;im2<nzmass;++im2)
		{
		  if(rbin[i]>rvir[im2]*0.75) 
		    {
		      ppp[i][j][k]+=ptemp[im1][im2][k]*ngal2[im1][im2]*
			(bias_interp(hmass[im1],-1)*bias_interp(hmass[im2],-1));
		      if(k==1)weightsum+=ngal2[im1][im2]*
			(bias_interp(hmass[im1],-1)*bias_interp(hmass[im2],-1));
		    }
		  continue;

		  if(rbin[i]/(rvir[im1]+rvir[im2])>0.5)
		    {
		      exfac=1;
		      if(rbin[i]/(rvir[im1]+rvir[im2])<1.5)
			exfac=ellipsoidal_exclusion_probability(rvir[im2]/rvir[im1],
								rbin[i]/(rvir[im1]+rvir[im2]));    
		      /*
		      if(rbin[i]<(rvir[im1]+rvir[im2]))
			exfac = 0;
		      if(rvir[im2]/rvir[im1]>pow(100.0,.333) && rbin[i]>(rvir[im1]+rvir[im2])*0.5)
			exfac = 1;
		      */
		      ppp[i][j][k]+=ptemp[im1][im2][k]*ngal2[im1][im2]*exfac*1*
			sqrt(bias_interp(hmass[im1],-1)*bias_interp(hmass[im2],-1));
		      if(k==1)weightsum+=ngal2[im1][im2]*exfac*1*
			sqrt(bias_interp(hmass[im1],-1)*bias_interp(hmass[im2],-1));
		    }
		}
	  //printf("WEIGHTSUM %d %d %e %.0f\n",
	  //	 i,j,weightsum,weightsum*BOX_SIZE*BOX_SIZE*BOX_SIZE);
	  if(weightsum>0)
	    {
	      for(k=1;k<=nv1;++k)
		ppp[i][j][k]/=weightsum;
	      //  for(k=1;k<=nv1;++k)
	      //printf("yTEST %d %d %f %e\n",i,j,vvv[k],ppp[i][j][k]);
	      for(ptot=0,k=1;k<=nv1;++k)
		ptot+=ppp[i][j][k]*binsize;
	    }
	  else
	    for(k=1;k<=nv1;++k)
	      ppp[i][j][k]=0;
	    
	}
    }
  
#ifdef PARALLEL
  for(i=1;i<=nr1;++i)
    for(j=1;j<=nphi1;++j)
      {
	for(k=1;k<=nv1;++k)
	  temp[k]=ppp[i][j][k];
	MPI_Allreduce(&temp[1],&precv[1],nv1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	for(k=1;k<=nv1;++k)
	  ppp[i][j][k]=precv[k];
      }
#endif

  free_d3tensor(ptemp,0,nzmass-1,0,nzmass-1,1,nv1);
  free_dvector(vtemp,1,nv1);
  free_dvector(temp,1,nv1);  
  free_dvector(p1temp,1,nv1);
  free_dvector(p2temp,1,nv1);
  free_dvector(precv,1,nv1);
  free_dvector(rvir,0,nzmass);
  free_dvector(mbin,0,nzmass);
  free_dvector(ngal,0,nzmass-1);
  free_dmatrix(ngal2,0,nzmass-1,0,nzmass-1);

}


