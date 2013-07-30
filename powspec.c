#include "header.h"

double sig8;         // value of sigma_8
double rho;          // current density of universe (solar mass/Mpc^3)
double G0;           // present day density cst
double Gl;           // Cosmological Constant densty
double H0;           // Hubble cst km/s.MPc
double pow_n;        // Scalar spectral index
double bfrac;        // Baryon fraction of total density
double Gamma;        // P(k) parameter - normally unused
double pow_norm;     // power spectrum normalisation
double NB_density;   // Nbody density being modelled
int pow_spec_type=2; // choice of power spectrum

// set up rho and power spectrum normalisation
void nm_power_eh()
{
  double r, mass, sigmasq;

  rho = 3.*G0*(3.236e-20*H0)*(3.236e-20*H0)/(8.*pi*G_cst*M_sun/(Mpc*Mpc*Mpc));
  fprintf(stderr,"h=%g, rho=%g\n",H0/100.,rho);
 
  r = 8.0;
  mass = 4./3.*pi*r*r*r*rho; 
  pow_norm=1.;
  sigmasq=qmidinf(Sigth_eh,1.0,1.e40,r)/(2.*pi*pi);
  pow_norm=(sig8*sig8)/sigmasq;
  fprintf(stderr,"r=%g, mass=%g, sig8=%g, nor=%g\n",r,mass,sig8, pow_norm);
}

// Sigmasq no filter : return D^2(k)
double Sig_eh(double lnk) { 
  double k;
  k = pow(10.,lnk);
  return k*k*k*Pk(k)/(2.*pi*pi); 
} 

// Sigmasq for a top hat filter - integrate 0..infinity
double Sigth_eh( double kplus1, double R )
{
  double k, x, W;

  k = kplus1-1.0;
  x = k*R;
  W = 3.*(sin(x) - x*cos(x))/(x*x*x);

  return W*W*k*k*Pk(k);
}

// Power spectrum
float Pk(float k)  
{  
  float pn=0.0;
  float q, tf, tfsq=0.0;  
  float tk_eh(float);
  if(k==0.) return 0.;

  switch(pow_spec_type) {

  // BBKS
  case 1:
    q = k/Gamma;
    tf = log(1.+2.34*q)/(2.34*q);
    tfsq = tf*tf/sqrt(1.0+q*(3.89+q*(259.21+q*(162.77+q*2027.17))));
    pn = pow_norm*pow(k,pow_n)*tfsq;
    break;

  // Eisenstein & Hu
  case 2:
    tf = tk_eh(k);
    tfsq = tf*tf;
    pn = pow_norm*pow(k,pow_n)*tfsq;
    break;

  // Power law k^(pow_n)
  case 3:
    tfsq = 1.0;
    pn = pow_norm*pow(k,pow_n)*tfsq;
    break;

  // DEFW P(k)
  case 4:
    q = k/Gamma;
    tf = pow(1.+pow(6.4*q+pow(3.0*q,1.5)+pow(1.7*q,2.),1.13),(-1./1.13));
    tfsq = tf*tf;
    pn = pow_norm*pow(k,pow_n)*tfsq;
    break;

  // baryonic fit
    //  case 5:
    //    pn = sin(fit_phi+k/fit_kc);
    //    pn = pow_norm*(1.+fit_a*k+fit_b*k*k)*pn*pn;
    //    break;

  // Eisenstein & Hu + non-linear approximation
  case 6:
    tf = tk_eh(k);
    tfsq = tf*tf;
    pn = pow_norm*pow(k,pow_n)*tfsq; 
    if(k>0.2) pn = exp( 2.0*(1.+log(k)-log(0.2))*log(pn) );
    break;

  // Eisenstein & Hu + shot noise
  case 7:
    tf = tk_eh(k);
    tfsq = tf*tf;
    pn = pow_norm*pow(k,pow_n)*tfsq;
    if(k<=(2.*2.*pi*pow(NB_density,1./3.))) pn += 1.0/NB_density;
    break;

  }

  return pn;
}

// Transfer function of Eisenstein & Hu 1998 
// (Equation numbers refer to this paper)

float tk_eh(float k)
{
  float rk,e,thet,thetsq,thetpf,b1,b2,zd,ze,rd,re,rke,s,rks,q,y,g;
  float ab,a1,a2,ac,bc,f,c1,c2,tc,bb,bn,ss,tb,tk_eh;
  float h,hsq,om_mhsq,om_b,om_m;

  // set up cosmology
  h    = H0/100.0;
  om_m = G0;
  om_b = bfrac*om_m;

  // convert k to Mpc^-1 rather than hMpc^-1
  rk=k*h;
  hsq=h*h;
  om_mhsq=om_m*hsq;

  // constants
  e=exp(1.);      
  thet=2.728/2.7;
  thetsq=thet*thet;
  thetpf=thetsq*thetsq;

  // Equation 4 - redshift of drag epoch
  b1=0.313*pow(om_mhsq,-0.419)*(1.+0.607*pow(om_mhsq,0.674));
  b2=0.238*pow(om_mhsq,0.223);
  zd=1291.*(1.+b1*pow(om_b*hsq,b2))*pow(om_mhsq,0.251)
    /(1.+0.659*pow(om_mhsq,0.828));

  // Equation 2 - redshift of matter-radiation equality
  ze=2.50e4*om_mhsq/thetpf;

  // value of R=(ratio of baryon-photon momentum density) at drag epoch
  rd=31500.*om_b*hsq/(thetpf*zd);

  // value of R=(ratio of baryon-photon momentum density) at epoch of matter-radiation equality
  re=31500.*om_b*hsq/(thetpf*ze);

  // Equation 3 - scale of ptcle horizon at matter-radiation equality
  rke=7.46e-2*om_mhsq/(thetsq);

  // Equation 6 - sound horizon at drag epoch
  s=(2./3./rke)*sqrt(6./re)*log((sqrt(1.+rd)+sqrt(rd+re))/(1.+sqrt(re)));

  // Equation 7 - silk damping scale
  rks=1.6*pow(om_b*hsq,0.52)*pow(om_mhsq,0.73)*(1.+pow(10.4*om_mhsq,-0.95));

  // Equation 10  - define q
  q=rk/13.41/rke;
      
  // Equations 11 - CDM transfer function fits
  a1=pow(46.9*om_mhsq,0.670)*(1.+pow(32.1*om_mhsq,-0.532));
  a2=pow(12.0*om_mhsq,0.424)*(1.+pow(45.0*om_mhsq,-0.582));
  ac=pow(a1,(-om_b/om_m))*pow(a2,pow(-(om_b/om_m),3.));

  // Equations 12 - CDM transfer function fits
  b1=0.944/(1.+pow(458.*om_mhsq,-0.708));
  b2=pow(0.395*om_mhsq,-0.0266);
  bc=1./(1.+b1*(pow(1.-om_b/om_m,b2)-1.));

  // Equation 18
  f=1./(1.+pow(rk*s/5.4,4.));

  // Equation 20
  c1=14.2 + 386./(1.+69.9*pow(q,1.08));
  c2=14.2/ac + 386./(1.+69.9*pow(q,1.08));

  // Equation 17 - CDM transfer function
  tc=f*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c1*q*q) +
    (1.-f)*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c2*q*q);

  // Equation 15
  y=(1.+ze)/(1.+zd);
  g=y*(-6.*sqrt(1.+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)));

  // Equation 14
  ab=g*2.07*rke*s/pow(1.+rd,0.75);

  // Equation 23
  bn=8.41*pow(om_mhsq,0.435);

  // Equation 22
  ss=s/pow(1.+pow(bn/rk/s,3.),1./3.);

  // Equation 24
  bb=0.5+(om_b/om_m) + (3.-2.*om_b/om_m)*sqrt(pow(17.2*om_mhsq,2.)+1.);

  // Equations 19 & 21
  tb=log(e+1.8*q)/(log(e+1.8*q)+c1*q*q)/(1+pow(rk*s/5.2,2.));
  tb=(tb+ab*exp(-pow(rk/rks,1.4))/(1.+pow(bb/rk/s,3.)))*sin(rk*ss)/rk/ss;
    
  // Equation 8
  tk_eh=(om_b/om_m)*tb+(1.-om_b/om_m)*tc;
  
  return tk_eh;
}

