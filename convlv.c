#define NRANSI
#include "nrutil.h"
#include "stdlib.h"
#include "stdio.h"

void convlv(double data[], unsigned long n, double respns[], unsigned long m,
	int isign, double ans[])
{
	void realft(double data[], unsigned long n, int isign);
	void twofft(double data1[], double data2[], double fft1[], double fft2[],
		unsigned long n);
	unsigned long i,no2;
	double dum,mag2,*fft;

	fft=dvector(1,n<<1);
	for (i=1;i<=(m-1)/2;i++)
		respns[n+1-i]=respns[m+1-i];
	for (i=(m+3)/2;i<=n-(m-1)/2;i++)
		respns[i]=0.0;
	/*
	for(i=1;i<=-n;++i)
	  fprintf(stdout,"RESP %d %e\n",i,respns[i]);
	fflush(stdout);
	*/
	twofft(data,respns,fft,ans,n);
	/*
	for(i=2;i<=n+2;i+=2)
	  fprintf(stdout,"FFT %d %e %e %e\n",i,ans[i-1],ans[i],DSQR(ans[i-1])+DSQR(ans[i]));
	*/
	for(i=2;i<=n+2;i+=2)
	  {
	    //ans[i] = 0;
	    //if(i>=60)ans[i-1] = 0;
	  }
	fflush(stdout);
	
	no2=n>>1;
	for (i=2;i<=n+2;i+=2) {
		if (isign == 1) {
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
			ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;
		} else if (isign == -1) {
		  if ((mag2=DSQR(ans[i-1])+DSQR(ans[i])) == 0.0) {
		    //printf("FFT %d %e %e %e\n",i,ans[i-1],ans[i],DSQR(ans[i-1])+DSQR(ans[i]));
		    nrerror("Deconvolving at response zero in convlv"); }
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;
			ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
		} else nrerror("No meaning for isign in convlv");
	}
	ans[2]=ans[n+1];
	//realft(fft,n,-1);
	realft(ans,n,-1);
	free_dvector(fft,1,n<<1);
}
#undef NRANSI
