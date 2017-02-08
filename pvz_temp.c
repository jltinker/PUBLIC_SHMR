#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "header.h"

double delta_pdf(double delta, double r);

void vdelta_v4(double m0, double m1, double rad, double theta, double binsize, 
	       double v0, int nv, double *a, double *pv, double *pt, double *pz,
	       double wgal[3], double sgal[4])
{
  static double pnorm=-1,flag=0;
  static double *rho,*mu,*pdelta;
  static int prev_cosmo=0;

  int i,j,k,nrho=75,halo_pair_only=0,use_file_coeff=0;
  double dlogrho,sintheta,costheta,tan2theta,rho_min,rho_max,w,b0,fac,
    collapse_factor,rho200,ptot,v,sigv,rho0,x1,x2,x3,rvir0,rho200t,sigvt,alpha,alphat,
    rvir1,rvirmax,rvir2,vfac,sr1,sr2,sr3,sr4,st1,st2,st3,st4,sz1,sz2,sz3,sz4,
    gnorm1,gnorm2,gnorm3,gnorm4,vv,pmax,rcrit,sig20,sig10,sig5;

  int direction,iend;
  double sigz1[76],sigz2[76],sigz3[76],sigz4[76],
    g1[76],g2[76],g3[76],g4[76],minp,xfit[3],yfit[3],afit,bfit;

  double t0,t1,tt1,tt2;
  int ii=0;

  double vfac0 = 0.97; /* factor om rho200[t,r]. was 0.97 */

  for(i=1;i<=nv;++i)
    /* pt[i]=pv[i]=*/ pz[i]=0;

  rvir0=(pow(3*m0/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),1.0/3.0) + 
	 pow(3*m1/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),1.0/3.0));
  rvir1=pow(3*m1/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
  rvir2=pow(3*m0/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
  b0=(bias_interp(m0,-1)+bias_interp(m1,-1));

  /* new spherical overdensity model
   */
  if(rad<rvir1*0.5)return;

  /* 0.5 comes from ellipsoidal halo exclusion.
   */
  //  if(rad<rvir0*0.5)return;
  /*if(rad<rvir0)rad=rvir0;*/
  rho_min=0.01;

  alpha=alphat=pow(rad/35,0.1);
  rcrit=4*rvir1;

  rcrit=0;
  if(rad<=rcrit)
    {
      x1=rcrit;
      x2=pow(x1/35,0.1);
      x3=x1*pow(x2,-1.0/(+0.3));
      alpha=alphat=pow(rad/x3,+0.3);
    }
  rcrit=4*rvir1;

  rho200 = pow(rad/5.0/sqrt(rvir1),-4.0) + pow(rad/11.5/sqrt(rvir0),-1.3) + 0.5;
  rho200t = pow(rad/7.2/sqrt(rvir1),-2.5) + pow(rad/12.6/sqrt(rvir0),-0.8) + 0.48;
  rho0 = 1.41*b0 + pow(rad/9.4/rvir1,-2.2);

  /* TEMP TEMP TEMP
   */
  /* rho0 = 1.38*b0 + pow(rad/9.4/rvir1,-2.2); */
  vfac0 = 1.01;

  /* If using the bias parameters that fit Warren's sims, use:
   */  
  rho0 = 1.43*b0 + pow(rad/9.4/rvir1,-2.2);
  vfac0 = 0.95;


  /* If using the bias parameters for the linear P(k), use:
   */
  if(LINEAR_PSP)
    rho0 = 1.28*b0 + pow(rad/9.4/rvir1,-2.2);
    

  /* The NEW SPHERICAL OVERDENSITY FIT
   */
  rho200 = pow(rad/5.0/sqrt(rvir1),-3.0) + pow(rad/11.5/sqrt(rvir0),-1.45) + 0.5;
  rho200t = pow(rad/7.2/sqrt(rvir1),-2.5) + pow(rad/12.6/sqrt(rvir0),-1.2) + 0.48; 
  rho0 = 1.43*b0 + pow(rad/9.4/rvir1,-2.3);
  vfac0 = 1.0;

  vfac0 = 1.2; //testing to see how much i need.
  rho200t = pow(rad/5.6/sqrt(rvir1),-2.5) + pow(rad/16.6/sqrt(rvir0),-1.2) + 1.2;
  rho200 = pow(rad/3.0/sqrt(rvir1),-3.0) + pow(rad/15.0/sqrt(rvir0),-1.45) + 1.0;
  vfac0 = 1.05;
  rho0 = 1.20*b0 + pow(rad/9.4/rvir1,-2.3);

  sintheta=sin(theta);
  costheta=cos(theta);
  tan2theta=sintheta*sintheta/(costheta*costheta);

  if(!flag)
    {
      flag=1;
      mu=dvector(1,nrho);
      rho=dvector(1,nrho);
      pdelta=dvector(1,nrho);
    }

  /* Determine the proper limits of the integral over rho
   */
  dlogrho=log(1.0E+3/0.01)/(50-1);
  pmax=0;
  i=0;

  for(i=1;i<=50;++i)
    {
      rho[i]=rho_min*exp((i-1)*dlogrho);
      pdelta[i]=exp(-rho0/rho[i])*delta_pdf(rho[i]-1,rad)*rho[i];
      if(pdelta[i]>pmax){ pmax=pdelta[i]; j=i; }
    }

  pmax*=1.0E-4;
  for(i=1;i<j;++i)
    if(pdelta[i]>=pmax)break;
  rho_min=rho[i];
  for(i=j+1;i<50;++i)
    if(pdelta[i]<=pmax)break;
  rho_max=rho[i];

  /* Tabulate the spherical collapse model, density, and p(delta).
   */
  dlogrho=log(rho_max/rho_min)/(nrho-1);
  ptot=0;

  for(i=1;i<=nrho;++i)
    {
      rho[i]=rho_min*exp((i-1)*dlogrho);
      pdelta[i]=exp(-rho0/rho[i])*delta_pdf(rho[i]-1,rad)*rho[i];
      ptot+=pdelta[i]*dlogrho;
      if(pdelta[i]>pdelta[j])j=i;

      x1=-(100*rad)*pow(OMEGA_Z,GAMMA)*(rho[i]-1)/3*HUBBLEZ;
      
      /* If underdense region, use linear theory.
       */
      if(rho[i]<=1)
	{
	  mu[i]=x1;
	  continue;
	}

      x2=(rad)*spherical_collapse_model(rho[i]-1)*exp(-(4.5/rad/rho[i])*(4.5/rad/rho[i]));
      if(x2>0)x2=0;

      /* If small separation, use spherical collapse model only.
       */
      if(rad<=4)
	{
	  mu[i]=x2;
	  continue;
	}

      /* If intermediate separation, use linear combination of v_sph & v_lin.
       */
      if(rad<=20)
	{
	  w=-0.62*log(rad)+1.86;
	  mu[i]=w*x2 + (1-w)*x1;
	  continue;
	}
      
      /* In linear regime. Use linear theory.
       */
      mu[i]=x1;
    }

  if(rad<=rcrit)
    for(i=1;i<=nrho;++i)
      mu[i]=mu[j];
  
  for(i=1;i<=nrho;++i)
    pdelta[i]/=ptot;

  vfac=pow(OMEGA_Z/0.3,GAMMA)*(SIGMA_8/0.8)*vfac0;

  /* If Halos only, then set all the galaxy velocity dispersions
   * and weights to be zero.
   */
  if(!HOD.pdfs)
    {
      sgal[0]=sgal[1]=sgal[2]=0;
      wgal[0]=wgal[1]=wgal[2]=0;
    }


  fac=0.4*sqrt(4.499E-48/2.0)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;
  fac = fac*fac;

  for(i=1;i<=nrho;++i)
    {
      sigv=200*pow(rho[i]/rho200,alpha)*vfac;
      sigvt=200*pow(rho[i]/rho200t,alphat)*vfac;
      //      sigv = sqrt(sigv*sigv + fac*pow(m1,0.6666) + fac*pow(m0,0.666));
      //sigvt = sqrt(sigvt*sigvt + fac*pow(m1,0.6666) + fac*pow(m0,0.666));


      sr1=sigv*sigv + sgal[0];
      sr2=sigv*sigv + sgal[1];
      sr3=sigv*sigv + sgal[2];
      sr4=sigv*sigv + sgal[3];

      st1=sigvt*sigvt + sgal[0];
      st2=sigvt*sigvt + sgal[1];
      st3=sigvt*sigvt + sgal[2];
      st4=sigvt*sigvt + sgal[3];

      sz1=2*(st1+tan2theta*sr1)*costheta*costheta;
      sz2=2*(st2+tan2theta*sr2)*costheta*costheta;
      sz3=2*(st3+tan2theta*sr3)*costheta*costheta;
      sz4=2*(st4+tan2theta*sr4)*costheta*costheta;

      sigz1[i]=sz1;
      sigz2[i]=sz2;
      sigz3[i]=sz3;
      sigz4[i]=sz4;

      gnorm1=wgal[0]/(RT2PI*sqrt(st1*sr1)*sqrt(1.0/(sr1) + tan2theta/(st1)));
      gnorm2=wgal[1]/(RT2PI*sqrt(st2*sr2)*sqrt(1.0/(sr2) + tan2theta/(st2)));
      gnorm3=wgal[2]/(RT2PI*sqrt(st3*sr3)*sqrt(1.0/(sr3) + tan2theta/(st3)));
      gnorm4=(1-wgal[0]-wgal[1]-wgal[2])/
	(RT2PI*sqrt(st4*sr4)*sqrt(1.0/(sr4) + tan2theta/(st4)));

      g1[i]=gnorm1;
      g2[i]=gnorm2;
      g3[i]=gnorm3;
      g4[i]=gnorm4;
    }

  /* Find the mode of the distribution. Start at vz=0 and work either way.
   */
  for(k=1,j=nv/2;j<=nv;++j,++k)
    {
      for(i=1;i<=nrho;++i)
	{
	  v=v0+(j-1)*binsize;
	  vv=(v-sintheta*mu[i])*(v-sintheta*mu[i]);
	  pz[j]+=pdelta[i]*dlogrho/costheta*
	    (g1[i]*exp(-vv/(sigz1[i])) + g2[i]*exp(-vv/(sigz2[i])) 
	     + g3[i]*exp(-vv/(sigz3[i])) + g4[i]*exp(-vv/(sigz4[i])));
	}
      if(k==2){
	if(pz[j-1]>pz[j])direction=-1;
	else direction=1;
	break;
      }
    }

  direction=1;
  if(direction==1)iend=nv;
  else iend=1;

  for(k=1,j=nv/2-direction;j!=iend;j+=direction,++k)
    {
      if(pz[j]>0)goto SKIP1;
      for(i=1;i<=nrho;++i)
	{
	  v=v0+(j-1)*binsize;
	  vv=(v-sintheta*mu[i])*(v-sintheta*mu[i]);
	  pz[j]+=pdelta[i]*dlogrho/costheta*
	    (g1[i]*exp(-vv/(sigz1[i])) + g2[i]*exp(-vv/(sigz2[i])) 
	     + g3[i]*exp(-vv/(sigz3[i])) + g4[i]*exp(-vv/(sigz4[i])));
	}
    SKIP1:
      if(k<3)continue;
      if(direction>0) {
	/*printf("%d %e %e %e\n",j,pz[j],pz[j-1],pz[j-2]);*/
	if(pz[j-1]>pz[j] && pz[j-1]>pz[j-2]) {
	  minp = pz[j-1]*0.001; break; }
      } else {
	/*printf("%d %e %e %e\n",j,pz[j],pz[j+1],pz[j+2]);*/
	if(pz[j+1]>pz[j] && pz[j+1]>pz[j+2]) {
	  minp = pz[j+1]*0.001; break; }
      }
    }

  for(j=nv/2;j<=nv;++j)
    {
      if(pz[j]>0)continue;
      for(i=1;i<=nrho;++i)
	{
	  v=v0+(j-1)*binsize;
	  vv=(v-sintheta*mu[i])*(v-sintheta*mu[i]);
	  pz[j]+=pdelta[i]*dlogrho/costheta*
	    (g1[i]*exp(-vv/(sigz1[i])) + g2[i]*exp(-vv/(sigz2[i])) 
	     + g3[i]*exp(-vv/(sigz3[i])) + g4[i]*exp(-vv/(sigz4[i])));
	}
      if(pz[j]<=minp && nv-j>=3)break;
    }

  
  /* Now extrapolate the rest of the PVD from the last three points.
   */
  for(k=0,i=j-2;i<=j;++i,++k)
    {
      yfit[k]=log(pz[i]);
      xfit[k]=v0+(i-1)*binsize;
    }
  least_squares(xfit,yfit,3,&afit,&bfit);


  for(i=j+1;i<=nv;++i) 
    pz[i] = exp(afit + bfit*(v0+(i-1)*binsize));

  /* Now go from the mode to i=0
   */
  for(j=nv/2-1;j>=1;--j)
    {
      if(pz[j]>0)continue;
      for(i=1;i<=nrho;++i)
	{
	  v=v0+(j-1)*binsize;
	  vv=(v-sintheta*mu[i])*(v-sintheta*mu[i]);
	  pz[j]+=pdelta[i]*dlogrho/costheta*
	    (g1[i]*exp(-vv/(sigz1[i])) + g2[i]*exp(-vv/(sigz2[i])) 
	     + g3[i]*exp(-vv/(sigz3[i])) + g4[i]*exp(-vv/(sigz4[i])));
	}
      if(pz[j]<=minp && j>=3)break;
    }

  
  /* Now extrapolate the rest of the PVD from the last three points.
   */
  for(k=2,i=j+2;i>=j;--i,--k)
    {
      yfit[k]=log(pz[i]);
      xfit[k]=v0+(i-1)*binsize;
    }
  least_squares(xfit,yfit,3,&afit,&bfit);

  for(i=j-1;i>0;--i)
    pz[i] = exp(afit + bfit*(v0+(i-1)*binsize));


  /* Do this just to get rid of any numerical roundoff errors and such.
   * (It whould already be one by construction if I didn't have the
   * linear extrapolation at the end of each side, but the correction
   * shouldn't be much.)
   */
  for(ptot=0,i=1;i<=nv;++i)
    ptot+=pz[i]*binsize;

  for(i=1;i<=nv;++i)
    pz[i]/=ptot;


  if(isnan(ptot))
    {
      printf("NAN ptot %f %f %f %f %f %f\n",rho200,rho200t,alpha,rho0,sr4,st4);
      fflush(stdout);
      for(i=1;i<=nv;++i)
	pt[i]=pz[i]=pv[i]=0;
      return;
    }


  return;

}

/* This is the standard log-normal 1-pt distribution of dark
 * matter densities (non-linear) at top-hat smoothing scale r.
 */
double delta_pdf(double delta, double r)
{
  static int model_prev=0;
  static double pnorm=-1,rprev=0,sig0,sig1sqr;
  double pnorm1;

  if(pnorm<0 || RESET_COSMOLOGY!=model_prev)
    {
      pnorm1=SIGMA_8/sigmac(8.0);
      /*printf("NORM %e %e %e\n",pnorm1,sigmac(8.0),nonlinear_sigmac(8.0));*/
      pnorm=pnorm1*sigmac(80.0)/nonlinear_sigmac(80.0);
      rprev=0;
      /*printf("NORM %e %e %e %e %f\n",pnorm,sigmac(8.0),sigmac(80.0),nonlinear_sigmac(80.0),SIGMA_8);*/
    }
  if(r!=rprev)
    {
      sig0=pnorm*nonlinear_sigmac(r);
      sig1sqr=log(1+sig0*sig0);
      /*printf("NORM %f %f %e %e %e\n",r,sig0,pnorm,nonlinear_sigmac(8.0),sigmac(8.0));*/
    }  
  rprev=r;
  model_prev=RESET_COSMOLOGY;
  return(1.0/(RT2PI*sqrt(sig1sqr))*exp(-pow(log((1+delta))+sig1sqr*0.5,2.0)/
				       (2*sig1sqr))/(1+delta));
  
}
