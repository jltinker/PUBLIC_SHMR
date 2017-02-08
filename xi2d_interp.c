#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#ifdef PARALLEL
#include <mpi.h>
#endif
#include "header.h"

void calc_xi2d(double **xi2d, int n, double rmax, double rmin);

void splin2(double x1a[], double x2a[], double **ya, double **y2a, 
	    int m, int n,double x1, double x2, double *y);
void splie2(double x1a[], double x2a[], double **ya, int m, int n, double **y2a);


double xi2d_interp(double rs1, double rp1, double rs2, double rp2)
{
  static double **xi2d,rmax,rmin,*rsig,*rpi,**yy2d;
  static int flag=1,n,prev;
  int i,j,ng=20;
  double dx,dy,dlogr,xsum,rs,rp,x;

  if(flag || prev!=RESET_COSMOLOGY)
    {
      two_halo(1,1);

      prev=RESET_COSMOLOGY;
      n=35;
      if(flag)
	{
	  xi2d=dmatrix(1,n,1,n);
	  yy2d=dmatrix(1,n,1,n);
	  rsig=dvector(1,n);
	  rpi=dvector(1,n);
	}
      flag=0;
      rmax=log(40.5);
      rmin=log(0.05);
      dlogr = (rmax-rmin)/(n-1);
      for(i=1;i<=n;++i)
	rsig[i]=rpi[i]=exp((i-1)*dlogr+rmin);
      calc_xi2d(xi2d,n,rmax,rmin);
      splie2(rsig,rpi,xi2d,n,n,yy2d);
    }
  if(rs2<0 && rp2<0)
    {
      splin2(rsig,rpi,xi2d,yy2d,n,n,(rs1),(rp1),&x);
      return(x);
    }
  

  dx=(rs2-rs1)/ng;
  dy=(rp2-rp1)/ng;
  xsum=0;
  for(i=1;i<=ng;++i)
    for(j=1;j<=ng;++j)
      {
	rs=rs1+(i-0.5)*dx;
	rp=rp1+(j-0.5)*dy;
	splin2(rsig,rpi,xi2d,yy2d,n,n,(rs),(rp),&x);
	xsum+=2*dy*rs*dx*(x);
      }
  xsum/=((rs2*rs2-rs1*rs1)*(rp2-rp1));

  return(xsum);

}

/* The angular grid here is set such that
 * it's only used for the small-r SDSS multipoles, i<=3, or r<~ 1 Mpc/h
 */
double xi2d_interp_polar(double rs1, double rs2, double phi1, double phi2)
{
  static double **xi2d,rmax,rmin,*rsig,*rpi,**yy2d;
  static int flag=1,n,prev;
  int i,j,ng=20;
  double dx,dy,dlogr,xsum,rs,rp,x,vsum,dr,dphi,phi,r;

  if(flag || prev!=RESET_COSMOLOGY)
    {
      two_halo(1,1);

      prev=RESET_COSMOLOGY;
      n=15;
      if(flag)
	{
	  xi2d=dmatrix(1,n,1,n);
	  yy2d=dmatrix(1,n,1,n);
	  rsig=dvector(1,n);
	  rpi=dvector(1,n);
	}
      flag=0;
      rmax=log(2.5);
      rmin=log(0.05);
      dlogr = (rmax-rmin)/(n-1);
      for(i=1;i<=n;++i)
	rsig[i]=rpi[i]=exp((i-1)*dlogr+rmin);
      calc_xi2d(xi2d,n,rmax,rmin);
      splie2(rsig,rpi,xi2d,n,n,yy2d);
    }

  rs1 = log(rs1);
  rs2 = log(rs2);
  dr=(rs2-rs1)/ng;
  dphi=(phi2-phi1)/ng;
  xsum=vsum=0;
  for(i=1;i<=ng;++i)
    for(j=1;j<=ng;++j)
      {
	r=exp(rs1+(i-0.5)*dr);
	phi=phi1+(j-0.5)*dphi;
	rp = r*sin(phi);
	rs = r*cos(phi);
	splin2(rsig,rpi,xi2d,yy2d,n,n,(rs),(rp),&x);
	xsum+=sin(phi)*r*r*r*(x);
	vsum+=r*r*r*sin(phi);
      }
  xsum/=vsum;

  return(xsum);

}

void calc_xi2d(double **xi2d, int n, double rmax, double rmin)
{
  static int ii=0;
  double dlogr,rs,rp,*send,*recv,x1,x2;
  int i1,j1,i,j,istart=1,istep=1;

#ifdef PARALLEL
  istep = NTask;
  istart = ThisTask + 1;
#endif

  printf("CPUz %d here\n",ThisTask);fflush(stdout);

  dlogr = (rmax-rmin)/(n-1);

  for(i=1;i<=n;++i)
    for(j=1;j<=n;++j)
      xi2d[i][j]=0;

  for(i1=istart;i1<=n*n;i1+=istep)
    {
      i = (i1-1)%n+1;
      j = (i1-1)/n+1;
      rs = exp((i-1)*dlogr+rmin);
      rp = exp((j-1)*dlogr+rmin);
      if(KAISER)
	{
	  xi2d[i][j] = kaiser_distortion(rs,rp);
	  continue;
	}
      x1 = one_halo(rs,rp);
      x2 = two_halo(rs,rp);
      xi2d[i][j] = (x1+x2);
      if(OUTPUT)
	printf("ZZ%d %03d %d %d %f %f %e %e %e\n",ii,i1,i,j,rs,rp,(xi2d[i][j]),x1,x2);
    }
#ifdef PARALLEL
  send=dvector(1,n);
  recv=dvector(1,n);
  for(i=1;i<=n;++i)
    {
      for(j=1;j<=n;++j)
	send[j]=xi2d[i][j];
      printf("CPUa %d %d\n",ThisTask,i);
      fflush(stdout);
      MPI_Allreduce(&send[1],&recv[1],n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      printf("CPUb %d %d\n",ThisTask,i);
      fflush(stdout);
      for(j=1;j<=n;++j)
	xi2d[i][j]=recv[j];
      }
#endif
}
