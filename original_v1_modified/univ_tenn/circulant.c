#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <fftw3.h>
#include "complex.h"
#include "random.c"

#define SEED 37

#define SYNTHETIC 0

#if SYNTHETIC
#include "param"
#else
#define GX_DIM    512
#define GY_DIM    640     /*  GX_DIM*GY_DIM = # rows of H ; WVDIM*XDIM*YDIM = # columns of H */
#define WVDIM (int)    1
#define XDIM  (int)    132
#define YDIM  (int)    132
#define DX    (int)     1
#define DY    (int)     1
#define GX    (int)  512
#endif

Cx        *Tc0,*Tc1;
fftw_plan PlanF,PlanI;

#define malloc(x)   fftw_malloc(x)
#define calloc(n,s) memset(malloc(n*s),0,n*s)

#define I_reverse_cpy(n,in,out) fftw_execute_dft(PlanI,in,out)
#define F_reverse_cpy(n,in,out) fftw_execute_dft(PlanF,in,out)

#define ifft(n,x) I_##x
#define fft(n,x)  F_##x

void error(int x, char *s, ...)
{
  va_list args;
  va_start(args,s);
  vfprintf(stderr,s,args);
  if (x) exit(x);
}

char *syntheticLine()
{
  static char b[64];  static int h,r, c = 0, v = 0, x = XDIM, y = YDIM;
  
  if (!v){
    if ((v = GX_DIM*GY_DIM-YDIM*GX) < 2) error(1,"syntheticLine\n");
    printf("v = %d\n",v);
  n:
    h = 1+rnd(v);
    r = 0;
    sprintf(b,"col %d\n",c++);
    return b;
  }
  if (r < h){
    sprintf(b,"%d %lf\n",r++,U01); // r th value in current column
    return b;
  }
  if (r++ == h){
    sprintf(b,"\n");               // end of current column
    return b;
  }
  r = 0;                           // reset row
  if (--x){
  t:
    sprintf(b,"col %d\n",c++);
    return b;
  }
  x = XDIM;                        // reset small block
  if (--y){
    goto t;
  }
  y = YDIM;                        // reset large block
  if (c == XDIM*YDIM*WVDIM)
    return NULL;
  goto n;
}

#if SYNTHETIC
void computeC(int n, char *c) // n = GX_DIM*GY_DIM
{
  char b[64]; Cx *f,*g; int h,i,j,k; FILE *file;

  if ((i = GX_DIM*GY_DIM-YDIM*GX) < 2) error(1,"syntheticLine\n");
  printf("i = %d\n",i);
  k = 1+rnd(i);

  f = calloc(n,sizeof(Cx));
  g = calloc(n,sizeof(Cx));

  for (h = 0; h < WVDIM; h++){
    for (j = 0; j < k; j++){
      f[j][0] = U01;
    }
    fft(n,reverse_cpy(n,f,g));          // transform column --> diagonal
    sprintf(b,"A_%d",h);
    file = fopen(b,"w");
    write(fileno(file),g,n*sizeof(Cx)); // transformed (diagonal) representation of A_...
    fclose(file);
    for (i = 0; i < n; i++){
      f[i][0] = 0.;
    }
  }
  free(f); free(g);
}

#else
/*

This function expects to see the columns of the system matrix H in 
the following format:

col x
y z
.
.

where x is a column number, y is a row number and z is a system matrix value.  
Only nonzero system matrix values are given.  The first column for each wavelength 
is loaded.  All other columns, if given, are ignored.

The inverse transform (as defined by the paper) of each column times a scaling factor 
is stored on disk.  In other words, 

sqrt(n) F* C_k e_0

The reason why calling the fftw forward transform results in the above computation 
is because of the difference between the way the paper defines the forward transform 
and the way the fftw library defines the forward transform.

Paper:

F_i,j = 1/sqrt(n) e^(2 pi sqrt(-1) ij/n)

fftw:

F_i,j = e^-(2 pi sqrt(-1) ij/n)

*/
void computeC(int n, char *c) // n = GX_DIM*GY_DIM
{
  char b[64]; double t; Cx *f,*g; int h,i,j; FILE *file0, *file1;

#if SYNTHETIC
  char *q;
#define fgets(...) ((q = syntheticLine())? strcpy(b,q): NULL)
#else
  if (!(file0 = fopen(c,"r"))) error(1,"computeC : can't open %s\n",c);
#endif

  f = calloc(n,sizeof(Cx));
  g = calloc(n,sizeof(Cx));
  h = -1;

  while (fgets(b,64,file0)){                        // read columns
    if (1 == sscanf(b,"col %d",&i)){
      //      printf("col %d\n",i);
      if (h < 0){
        if (i%(XDIM*YDIM)) error(0,"first col = %d, XDIM*YDIM = %d\n",i,XDIM*YDIM);
        h = i/(XDIM*YDIM);
      }
      while (fgets(b,64,file0)){                    // column entries
        if (b[0] == '\n') break;

        if (i%(XDIM*YDIM)) continue;                // not first column of new group
        if (2 == sscanf(b,"%d %lf",&j,&t))
          f[j][0] = t;                              // save column
        else
          error(1,"computeC i = %d, j = %d\n",i,j);
      }
      if (!(i%(XDIM*YDIM))){
        h++;                                        // next slice
        //      printf("%d\n",h);
      }
    } else
      error(1,"computeC i = %d\n",i);
    if (!(i%(XDIM*YDIM))){
      fft(n,reverse_cpy(n,f,g));                    // transform column --> diagonal
      sprintf(b,"A_%d",h-1);
      file1 = fopen(b,"w");
      write(fileno(file1),g,n*sizeof(Cx));          // transformed (diagonal) representation of A_...
      fclose(file1);
      for (i = 0; i < n; i++){
        f[i][0]  = f[i][1] = 0.;
      }
    }
  }
#if !SYNTHETIC
  fclose(file0);
#endif  
  free(f); free(g);
}

#endif

int main(int argc, char **argv)
{
  char b[124]; int n;  FILE *file;
 
  n = GX_DIM*GY_DIM;                 // # rows in H

  Tc0 = malloc(n*  sizeof(Cx    ));
  Tc1 = malloc(n*  sizeof(Cx    ));
  sprintf(b,"wisdom_%d",n);
  if ((file = fopen(b,"r"))){
    fftw_import_wisdom_from_file(file);
    fclose(file);
  }
  PlanF = fftw_plan_dft_1d(n,Tc0,Tc1,FFTW_FORWARD,FFTW_MEASURE);
  PlanI = fftw_plan_dft_1d(n,Tc0,Tc1,FFTW_BACKWARD,FFTW_MEASURE);
  if ((file = fopen(b,"w"))){
    fftw_export_wisdom_to_file(file);
    fclose(file);
  }
  initrand(SEED);
  computeC(GX_DIM*GY_DIM,"save-columns-Hmatrix.h!");
  return 0;
}
