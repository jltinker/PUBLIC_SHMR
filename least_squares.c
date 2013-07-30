#include "header.h"

/* Just a little ditty to do a least squares fit to the arrays given
 * (assuming no errors).
 */
void least_squares(double *x, double *y, int n, double *a, double *b)
{
  int i,j;
  double delta,sx=0,sy=0,sxx=0,sxy=0,syy=0;
  
  for(i=0;i<n;++i)
    {
      sx+=x[i];
      sy+=y[i];
      sxx+=x[i]*x[i];
      syy+=y[i]*y[i];
      sxy+=x[i]*y[i];
    }
  delta=n*sxx-sx*sx;
  *a=(sxx*sy-sx*sxy)/delta;
  *b=(n*sxy-sx*sy)/delta;
}
