double qp_input2a()
{
  int i,j,k,n,nfiles=1000,itail;
  char fname[1000], aa[1000], froot[100], fname1[1000], command[1000];
  float x1,x2,x3,xx[100],xibar,xb[100];
  FILE *fp;

  sprintf(froot,"ximulti");

  if(ISKIP_RSD)
    {
      sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_XI/%s_%04d.dat",froot,1);
      fp = openfile(fname);
      rsdq.np = filesize(fp);
      fclose(fp);
      sprintf(froot,"ximulti_xx");
      itail = rsdq.np - ISKIP_RSD;
      for(i=1;i<=nfiles;++i)
	{
	  sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_XI/ximulti_%04d.dat",i);
	  sprintf(fname1,"/Users/tinker/NSF_2013_WORK/QPM_XI/%s_%04d.dat",froot,i);
	  sprintf(command,"tail -%d %s > %s",itail,fname,fname1);
	  fprintf(stdout,"input_rsd> [%s]\n",command);
	  system(command);
	} 
    }

  // get the number of data points;
  sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_XI/%s_%04d.dat",froot,1);
  fp = openfile(fname);
  rsdq.np = filesize(fp);
  if(rsdq.np>ICUT_RSD)rsdq.np = ICUT_RSD;
  rsdq.covar = dmatrix(1,rsdq.np,1,rsdq.np);
  rsdq.r = dvector(1,rsdq.np);
  rsdq.x = dvector(1,rsdq.np);
  rsdq.e = dvector(1,rsdq.np);



  for(i=1;i<=rsdq.np;++i)
    for(j=1;j<=rsdq.np;++j)
      rsdq.covar[i][j] = 0;

  n = nfiles;
  xibar = 0;
  for(j=1;j<=rsdq.np;++j)
    {
      fscanf(fp,"%f %f %f",&x1,&x2,&x3);
      fgets(aa,1000,fp);
      xibar += x1*x1*x2;
      rsdq.x[j] = x3/(x2-3/(j*j*j)*xbiar);
      rsdq.r[j] = x1; 
    }
  fclose(fp);


  for(i=2;i<=n;++i)
    {
      xibar = 0;
      sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_XI/%s_%04d.dat",froot,i);
      fp = openfile(fname);
      for(j=1;j<=rsdq.np;++j)
	{
	  fscanf(fp,"%f %f %f",&x1,&x2,&x3);
	  fgets(aa,1000,fp);
	  xibar += x1*x1*x2;
	  rsdq.x[j] = x3/(x2-3/(j*j*j)*xbiar);
	  rsdq.r[j] = x1; 
	}
      fclose(fp);
      //muh(i);
    }

  for(k=1;k<=n;++k)
    {
      sprintf(fname,"/Users/tinker/NSF_2013_WORK/QPM_XI/%s_%04d.dat",froot,k);
      fp = openfile(fname);
      for(i=1;i<=rsdq.np;++i)
	{
	  fscanf(fp,"%f %f %f",&x1,&x2,&x3);
	  xibar += x1*x1*x2;
	  xx[j] = x3/(x2-3/(j*j*j)*xbiar);
	  fgets(aa,1000,fp);
	}
      fclose(fp);
	    
      for(i=1;i<=rsdq.np;++i)
	for(j=1;j<=rsdq.np;++j)
	    rsdq.covar[i][j] += (xx[i]-rsdq.x[i])*(xx[j]-rsdq.x[j])/n;
    }
  for(i=1;i<=rsdq.np;++i)
    rsdq.e[i] = sqrt(rsdq.covar[i][i]);

  for(i=1;i<=-rsdq.np;++i)
    printf("INPUT %e %e %e\n",rsdq.r[i],rsdq.x[i],rsdq.e[i]);

  // lets do some PCA
  rsdq.eval = dvector(1,rsdq.np);
  rsdq.evect = dmatrix(1,rsdq.np,1,rsdq.np);
  rsdq.npca = 50;
      
}


