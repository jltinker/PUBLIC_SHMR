#include <stdio.h>
#include <math.h>
#include "header.h"

/* 2005-1-9. New Romberg integrator.  Supply limits (a,b) and two flags,
   normally (1,1) meaning the integral is closed; the function can be 
   evaluated right on the limits (a,b).  If the flags are set to zero then
   the function cannot be evaluated right on that limit and this routine
   performs an open (or semi-open) interval integration.

   More details: The code can be divided into two parts, the trapezoidal
   integrator and "Richardson's deferred approach to limit." See NR Fig 4.2.1: 
   the row labelled "N=2" is where we start with k=0.  Then the next row 
   (labelled N=3) is k=1, and so on, so that the number of function 
   evaluations per k step is n = 2^k.  dx is the spacing between the points.
   fxna,fxnb,fxn0, and fxn1 are basically temporary variables which are 
   shuffled around depending on the closed flags.  

   2005-1-10 adding a new criterion for convergence, because occasionally
   I get convergent wrong answers.  Require the last two answer[k]'s to 
   agree to a certain fraction.  This doesn't have be a very small number.
   Hm, actually it does.  We are dealing with an function that can be very
   tricky if you sample at even intervals.  In order to improve, ideally
   you would couple this integrator with a monte carlo abscissa integrator,
   and make sure the results match.
*/

#define kmaxx 21
#define nfit 8

double integrate(double (*fxn)(double), double a, double b,int closed_a,int closed_b) {
  int i,j,k,n;
  double fxna,fxnb,dx,halfdx,fxn0,fxn1,sum,totsum,tentative,cdterm;
  double answer[kmaxx],h[kmaxx],c[nfit],d[nfit];
  double tolerance = 1.e-7, tolerance2=1.e-4, answer_discrep;

  n = 1;
  totsum = 0.;
  h[0] = 1.;
  if (closed_a) fxna = (*fxn)(a);
  if (closed_b) fxnb = (*fxn)(b);
  for (k=0;k<kmaxx;k++) {
    if (k) h[k] = .25 * h[k-1];
    dx = (b-a) / (double)n;
    halfdx = dx/2.;
    fxn0 = (*fxn)(a + halfdx);
    fxn1 = (*fxn)(b - halfdx);
    if (n==1) sum = fxn0; else sum = fxn0 + fxn1;
    for (i=1;i<n-1;i++) {
      sum += (*fxn)(a + halfdx + i*dx);
    }
    totsum += sum;
    if (!closed_a) fxna = fxn0;
    if (!closed_b) fxnb = fxn1;
    answer[k] = (dx/2.) * (totsum + .5*(fxna+fxnb)); // the tentative answer
    tentative = answer[k];
    if(HOD.free[0]==-100)
      printf("%d\t%g\t%g\t%f\n",k,answer[k],a,b);
    if (k >= nfit-1) { // This is the Neville part (NR 3.1)
      answer_discrep = fabs(answer[k]-answer[k-1]);
      if (answer_discrep <= tolerance2*fabs(tentative)) {// passes my new test
      for (i=0;i<nfit;i++) {
	c[i] = answer[k-nfit+1 + i];
	d[i] = answer[k-nfit+1 + i];
      }
      for (i=1;i<nfit;i++) {
	for (j=0;j<nfit-i;j++) {
	  cdterm = (c[j+1]-d[j]) / (h[j]-h[j+i]);
	  d[j] = h[j+i]*cdterm;
	  c[j] = h[j]  *cdterm;
	}
	tentative += d[j-1];
	if (fabs(d[j]) < tolerance*fabs(tentative)) return tentative;
      } }
    }
    n *= 2;
  }

  fprintf(stderr,"too many romberg steps\n");
  return 0.;
}

#undef kmaxx
#undef nfit
