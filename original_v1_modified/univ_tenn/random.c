#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "random.h"

void freedist( distribution *d ) /* frees dist */
{
	free( d->a );
	free( d->p );
	free( d );
}

void sv_rnd(Ulong *v) /* save state; also saves for normal, thus v[59] */
{
  int i;

  for (i = 0; i < 55; i++) v[i] = rtab[i];
  v[i]   = rndx;
  v[i+1] = nflg; *((double *)(v+i+2)) = nv;
}

void rst_rnd(Ulong *v) /* restore state; also restores for normal, thus v[59] */
{
  int i;

  for (i = 0; i < 55; i++) rtab[i] = v[i];
  rndx = v[i];
  nflg =  v[i+1]; nv = *((double *)(v+i+2));
}

int nrndm()
{
  int i;

  for (i =  0; i < 24; i++) rtab[i] -= rtab[i+31];
  for (i = 24; i < 55; i++) rtab[i] -= rtab[i-24];
  return 0;
}

void initrand(Ulong j) /* Initialize the 32 bit random number generator with seed j */
{
  int h,i;  Ulong k;

  for (rtab[54] = j |= (k = i = 1); i < 55; i++)
    h = (21*i)%55, rtab[--h] = k, k = j - k, j = rtab[h];

  while (i--) nrndm();
  rndx = 0;
}

distribution *allocdist(long n) /* allocate a distribution over {0..n-1} */
{
  distribution *d;

  d = (distribution *) malloc(sizeof(distribution));
  d->n = n;
  d->a = (Ulong  *)malloc(d->n * sizeof(Ulong ));
  d->p = (double *)malloc(d->n * sizeof(double));
  return d;
}

#define getsmall { while (p[j] >= q) if ((++j) == stop) goto end; t = j++; }
#define getlarge   while (p[k] <  q) if ((++k) == stop) goto cleanup;

distribution *initdist(distribution *d, double s) /* Initialize the distribution d */
{
  /*
    Note: d->p must have d->n elements which sum to s on entry to initdist.
    The elements of d->p and d->a are overwritten by the initialization process.
  */

  Ulong j,k,t,stop,*a;  double q,*p;

  stop = d->n, q = s/stop, j = k = 0;

  d->m1 = stop/TWO_LONGSZ;
  d->m2 = s/(stop * TWO_LONGSZ);

  a = d->a;
  p = d->p;

  getsmall; getlarge;

 loop:
    
  a[t]  = k;
  p[k] += p[t] - q;

  if (p[k] >= q) { 
    if (j == stop) goto end;
    getsmall;
    goto loop;
  }
  t = k++;
  if (k == stop) goto cleanup;
  if (j < k) getsmall;
  getlarge;
  goto loop;

 cleanup:

  a[t] = t;
  while (j < stop) { a[j] = j; j++; }

 end: 
  return d;
}

Ulong drand(distribution *d) /* Return element from {0..d->n-1} according to d */
{
  int j;  Ulong r1 = rndm;

  if ((rndm+(r1&0xffff)/65536.0)*(d->m2) < d->p[j=r1*(d->m1)]) return j;
  return d->a[j];
}

double normal(double m, double s) // normaly distributed; mean m standard deviation s
{                                 // nflg nd nv are globals (so state can easily be saved/restored)
  double u,r;

  if (!nflg){
    do
      u = 2*U01 - 1., nv = 2*U01 - 1.;
    while ((r = u*u + nv*nv) >= 1 || (r == 0));
    r = sqrt( -2.*log(r)/r );
    u *= r;
    nv *= r;
    nflg = 1;
    return s*u + m;
  } else {
    nflg = 0;
    return s*nv + m;
  }
}

double gammln(double xx)
{
  double cof[6]={76.18009173,-86.50532033,24.01409822,-1.231739516,.00120858003,-.00000536382};
  double stp=2.50662827465,half=.5,one=1.,fpf=5.5,x,tmp,ser;
  int i;

  x=xx-one;
  tmp=x+fpf;
  tmp=(x+half)*log(tmp)-tmp;
  ser=one;
  for (i = 0; i < 6; i++){
    x=x+one;
    ser=ser+cof[i]/x;
  }
  return tmp+log(stp+ser); 
}

int poisson(double xm)
{
  double g,em,t,sq,alxm,y; int tmp;

  if (xm < 12.) {
    g=exp(-xm);
    em=-1.;
    t=1.;
  top: 
    em=em+1.;     
    t=t*U01;
    if (t > g) goto top;
  } else {
    sq=sqrt(2.*xm);
    alxm=log(xm);
    g=xm*alxm-gammln(xm+1.);
  o:  
    y=tan(M_PI*U01);
    em=sq*y+xm;
    if (em < 0.) goto o;
    tmp=em;
    em=tmp;
    t=.9*(1.+y*y)*exp(em*alxm-gammln(em+1.)-g);
    if (U01 > t) goto o;    
  }  
  tmp=em;
  return tmp;
}

#if 0
int main (int argc, char **argv)
{
  int i,count=0; double mean=atof(argv[2]),desired=atof(argv[3]);

  initrand(7);

  for (i = 0; i < atoi(argv[1]); i++) if (poisson(atof(argv[2])) == atol(argv[3])) count++;
  printf("%f %f\n",pow(mean*exp((-mean/desired)+1.)/desired,desired)/sqrt(2.*PI*desired),(double)count/(double)atoi(argv[1]));
  
  //printf("%d\n",poisson(atof(argv[1])));
  return 0;
}
#endif
