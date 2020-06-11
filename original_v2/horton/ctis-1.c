#define  _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include <semaphore.h>
#include <sys/types.h>
#include <signal.h> 
#include <sys/stat.h>
#include <sys/time.h>
#include <dirent.h>              
#include <fftw3.h>
#include "complex.h"
#include "random.c"



// these #defines are to be removed later

#ifndef GENERATE_G
  #define GENERATE_G   0                                         // generate synthetic g (focal plane image) with noise
#endif

#ifndef SYNTHETIC_F
  #define SYNTHETIC_F  "synthetic_f.pgm"                         // 24 copies of the same image
#endif              



// these #defines stay

#ifndef THREADS
  #define THREADS 1                                              // number of threads
#endif

#ifndef FTHREADS
  #define FTHREADS 1                                             // number of threads for fftw
#endif

#ifndef FOCAL_PLANE
  #define FOCAL_PLANE  "focal_plane"                             // g
#endif

#ifndef COMPUTE_PT
  #define COMPUTE_PT  1                                          // compute P^T and store to disk
#endif

#ifndef REMOVE_NOISE
  #define REMOVE_NOISE 0                                         // remove noise 
#endif

#ifndef REPORT
  #define REPORT       1                                         // report error at each cg step
#endif

#ifndef PICTURE_OUT                                                                
  #define PICTURE_OUT  "reconstructed_image"                     // reconstructed image
#endif

#ifndef PICTURE
  #define PICTURE      1                                         // output reconstructed image to pgm file
#endif

#ifndef SEED
  #define SEED         7                                        // random number generator seed
#endif

#ifndef EPSILON
  #define EPSILON      (1e-8)                                    // cg error threshold
#endif

#ifndef MU
  #define MU           .01                                       // regularization parameter
#endif

#ifndef LIM
  #define LIM           4                                        // cg iteration threshold
#endif

#ifndef RESTART
  #define RESTART       4                                        // # of vose-horton steps
#endif

#include "param"                                                // focal plane size; field stop size; system matrix 'gap' size

Cx             **A,**C,*Tc0,*Tc1,*Tc2,*Tc3[THREADS];            // diagonalized representation of script_H, P^T; complex scratch space
double          *Vec,*Td4,*B,**M,*Alpha;        
volatile double *In,*Out,E;
volatile int     J,Z;
sem_t            S[(THREADS)-1];                                // thread semaphores (start)
sem_t            F[(THREADS)-1];                                // thread semaphores (finished)
fftw_plan        PlanF,PlanI;                                   // fftw forward, reverse fft setup
volatile int     Lim,Quit,R;

pthread_mutex_t  Mutex;

double norm (double *v, int n)
{
  double nrm=0.0; int i;

  for (i = 0; i < n; i++){
    nrm += v[i]*v[i];
  }
  return sqrt(nrm);
}

#define malloc(x)   fftw_malloc(x)              
#define calloc(n,s) memset(malloc(n*s),0,n*s)     

void quit(int s)
{
  Quit = Lim = R = 0;
}

void error(int x, char *s, ...)
{
  va_list args;
  va_start(args,s);
  vfprintf(stderr,s,args);
  if (x) exit(x);
}

/*

This function does a number of things:

1) malloc space

2) read the diagonal representation of the system matrix off the disk

   Dbar = sqrt(n) F* C e_0

3) read fftw wisdom from off the disk

   fftw is capable of remembering information from previous runs to make 
   the current run faster

*/

void init(int n, int m, int w)
{
  char b[64]; int i,k; FILE *file;

  if (n   < m) error(1,"n (%d) < m (%d)\n",n,m);
  if (n*w < R) error(1,"n (%d) < R (%d)\n",n,R);


  if (fftw_init_threads() == 0){
    printf("error with fftw_init_threads\n");
    exit(1);
  }
  fftw_plan_with_nthreads(FTHREADS);

  pthread_mutex_init(&Mutex,NULL);

  for (i = 0; i < (THREADS)-1; i++){
    sem_init(S+i,0,0);
    sem_init(F+i,0,0);
  }
  initrand(SEED);
  A      = malloc(w*sizeof(Cx     *));
  C      = malloc(w*sizeof(Cx     *));
  Vec    = malloc(w*n*sizeof(double));
  M      = malloc(R*sizeof(double *));
  B      = malloc(R*sizeof(double  ));
  Alpha  = malloc(R*sizeof(double  ));
  Tc0    = malloc(n*sizeof(Cx      ));
  Tc1    = malloc(n*sizeof(Cx      ));
  Tc2    = malloc(n*sizeof(Cx      ));
  Tc3[0] = malloc(n*sizeof(Cx      ));
  Td4    = malloc(n*sizeof(double  ));
  for(i = 0; i < R; i++) M[i] = malloc(R*sizeof(double));
  for(i = 0; i < w; i++){
    C[i] = malloc(n*sizeof(Cx));                                 // "diagonal" representation of P^T
    A[i] = calloc(n,sizeof(Cx));                                 // "diagonal" representation of script_H
    sprintf(b,"A_%d",i+5);

printf("A_%d\n",i+5);

    if ((file = fopen(b,"r"))){
      k=read(fileno(file),A[i],n*sizeof(Cx));     
      fclose(file);
    } else {
      printf("could not open %s\n",b);
    }
  }
  sprintf(b,"wisdom_%d",n);
  if ((file = fopen(b,"r"))){
    fftw_import_wisdom_from_file(file);
    fclose(file);
  } else {
    printf("could not open %s\n",b);
  }

printf("making fftw_plan\n");

#if 0
  PlanF = fftw_plan_dft_1d(n,Tc0,Tc1,FFTW_FORWARD,FFTW_EXHAUSTIVE);
  PlanI = fftw_plan_dft_1d(n,Tc0,Tc1,FFTW_BACKWARD,FFTW_EXHAUSTIVE);
#else

#if 0
  PlanF = fftw_plan_dft_1d(n,Tc0,Tc1,FFTW_FORWARD,FFTW_MEASURE);
  PlanI = fftw_plan_dft_1d(n,Tc0,Tc1,FFTW_BACKWARD,FFTW_MEASURE);
#else
  PlanF = fftw_plan_dft_1d(n,Tc0,Tc1,FFTW_FORWARD,FFTW_ESTIMATE);
  PlanI = fftw_plan_dft_1d(n,Tc0,Tc1,FFTW_BACKWARD,FFTW_ESTIMATE);
#endif

#endif

printf("leaving init\n");

}

void swap_out(char *s, void *a, int n)
{
  FILE *f = fopen(s,"w"); int k;

  if (!f) error(1,"could not open %s for writing\n",s);
  k=write(fileno(f),a,n);
}

void swap_in(char *s, void *a, int n)
{
  FILE *f = fopen(s,"r"); int k;

  if (!f) error(1,"could not open %s for reading\n",s);
  k=read(fileno(f),a,n);
}

// fftw_execute_dft

#define I_reverse_cpy(n,in,out) fftw_execute_dft(PlanI,in,out)
#define F_reverse_cpy(n,in,out) fftw_execute_dft(PlanF,in,out)

#define ifft(n,x) I_##x
#define fft(n,x)  F_##x

/*

This function multiplies the m * 1 input vector in  by the n * m 
system matrix H to produce the n * 1 output vector out.

It does so by:

1) spreading (adding zeros in the right places) in such that it is w*n * 1
2) multiplying each n*1 piece times a different circulant matrix 
3) adding results together into an n*1 vector

The above is another way of saying that instead of performing the multiplication: Hf, 
the multiplication H_script (I_w cross E) f is performed.  The equivilance of these two 
multiplications is worked out in the paper.

===

When multiplying a circulant matrix times a vector, the following is leveraged:

C = F* D F

conjugate both sides =>

Cbar = F Dbar F*

Since C is real =>

C = F Dbar F* 

D = sqrt(n) diag (F C e_0)

conjugate both sides =>

Dbar = sqrt(n) diag (F* C e_0)

Now =>

C = F sqrt(n) diag (F* C e_0) F* 

The F above is the forward transform as defined by the paper.  However, 

the forward transform as defined by the fftw library is different.

---

Paper:
    
F_i,j = 1/sqrt(n) e^(2 pi sqrt(-1) ij/n)

fftw:

F_i,j = e^-(2 pi sqrt(-1) ij/n)

---

The rightmost F* above is computed in the code by calling the fftw forward transfrom.  

A[j] in the code is sqrt(n) diag (F* C e_0).

The leftmost F above is computed in the code by calling the fftw inverse transform.

Because the forward and inverse fftw transforms also differ from the paper by 1/sqrt(n), 
the code multiplies by 1/n.

*/
double *multiply_H(double *in, double *out, int n, int m, int w)
{                                                           
  int i,j,k;  double s = 1./(double)n;  Cx t;           
                                                         
  for (k = j = 0; j < w; j++){
    for (i = 0; i < n; i++)
      if ((i/GX >= YDIM)||(XDIM <= (i%GX)))
        Tc1[i][0] = Tc1[i][1] = 0.;
      else {
        Tc1[i][0] = in[k++]*s; Tc1[i][1] = 0.;                   // * 1/n because fft and ifft each missing 1/sqrt(n)
      }
    fft(n,reverse_cpy(n,Tc1,Tc0));                               // F* as defined by paper (except for 1/sqrt(n))
    if (!j)
      for (i = 0; i < n; i++) c_mul(A[j][i],Tc0[i],Tc2[i]);      // multiply A[j]: sqrt(n) F* C e_0
    else 
      for (i = 0; i < n; i++)
        c_add(Tc2[i],c_mul(A[j][i],Tc0[i],t),Tc2[i]);            // multiply A[j] and accumulate
  }
  ifft(n,reverse_cpy(n,Tc2,Tc1));                                // F as defined by paper (except for 1/sqrt(n))
  for (i = 0; i < n; i++) out[i] = Tc1[i][0];
  if (k != m) error(1,"multiply_H: m = %d, k = %d\n",m,k);
  return out;
}

#define OUT(j,i) out[n*j+i]

/*

This function computes script_H^T x

---

C = F* D F

conjugate transpose =>

C^T = F* Dbar F <-- C is real, D is symmetric

conjugate =>

C^T = F D F* <-- C is real

The above is precisely what is being done in the code.

See circulant.c to understand why A[j] is Dbar for the jth circulant matrix of H_script.

See multiply_H to understand why the calls to fftw below are equivalent to the above computation.

*/

void embed(double *in, double *out, int n, int w)                // embed and myltiply by A^T
{
  int i,j;  double s = 1./(double)n;

  for (i = 0; i < n; i++){
    Tc1[i][0] = in[i]*s; Tc1[i][1] = 0.;                         // * 1/n because fftw missing 1/sqrt(n) for F, F*
  }                                                      
  fft(n,reverse_cpy(n,Tc1,Tc0));                                 // F* as defined in the paper (except for 1/sqrt(n))
  for (j = 0; j < w; j++){
    for (i = 0; i < n; i++) c_cjm(A[j][i],Tc0[i],Tc1[i]);        // A[j] = Dbar_j = sqrt(n) F* C_j e_0; conjugate to get D
    ifft(n,reverse_cpy(n,Tc1,Tc2));                              // F as defined in the paper (except for 1/sqrt(n))
    for (i = 0; i < n; i++) OUT(j,i) = Tc2[i][0];         
  }
}

#define IN(j,i) in[n*j+i]

/*

This function performs the computation script_H^T script_H y_infinity.

See inv_C for an explanation of why it is that conjugation yields the transpose of the 
the diagonal representation of a circulant matrix.

See multiply_H for explanation of why it is that the computations below correspond to 
the multiplication of circulant matrices by vectors.

*/

void multiply_A(double *in, double *out, int n, int w) // sans MU * I ("H^T H" but in "higher" dimension)
{
  int i,j,k;  double s = 1./(double)n;  Cx t;

  In = in; Out = out; 

  pthread_mutex_lock(&Mutex);
  J = 0;
  pthread_mutex_unlock(&Mutex);

  for (i = 0; i < ((THREADS)-1); i++){
    sem_post(S+i);
  }

  for (j = 0; j < w; j++){
    if (j%(THREADS) == 0){
      for (i = 0; i < n; i++){
        Tc1[i][0] = IN(j,i)*s; Tc1[i][1] = 0.;                   // complexify input
        OUT(j,i) = 0.;
      }
      fft(n,reverse_cpy(n,Tc1,Tc0));                             // transform input
      if (!j)
        for (i = 0; i < n; i++) c_mul(A[j][i],Tc0[i],Tc3[0][i]); // multiply A[j]
      else 
        for (i = 0; i < n; i++)
          c_add(Tc3[0][i],c_mul(A[j][i],Tc0[i],t),Tc3[0][i]);    // multiply A[j] and accumulate
    }
  }
  for (i = 0; i < ((THREADS)-1); i++){
    sem_wait(F+i);
  }

  for (i = 0; i < n; i++){                                       // global accumulate
    for (k = 1; k < THREADS; k++){
      c_add(Tc3[0][i],Tc3[k][i],Tc3[0][i]);
    }
  }
  for (i = 0; i < ((THREADS)-1); i++){
    sem_post(S+i);
  }

  //sem_post(S+0); sem_post(S+1); sem_post(S+2);
  for (j = 0; j < w; j++){
    if (j%(THREADS) == 0){
      for (i = 0; i < n; i++) c_cjm(A[j][i],Tc3[0][i],Tc2[i]);   // multiply conjugate A[j] (equivalent to A^T)
      ifft(n,reverse_cpy(n,Tc2,Tc1));
      for (i = 0; i < n; i++) OUT(j,i) += Tc1[i][0];
    }
  }
  for (i = 0; i < ((THREADS)-1); i++){
    sem_wait(F+i);
  }

}

void multiply_A_thread(int k, int n, int w, Cx *tc1, Cx *tc2, Cx *tc3, double *in, double *out)
{
  int h,i,j; double s = 1./(double)n;  Cx t;

  for (h = j = 0; j < WVDIM; j++){
    if (j%(THREADS) == k+1){
      for (i = 0; i < n; i++){
        tc1[i][0] = IN(j,i)*s; tc1[i][1] = 0.;                   // complexify input
        OUT(j,i) = 0.;
      }
      fft(n,reverse_cpy(n,tc1,tc2));                             // transform input
      if (!h++)
        for (i = 0; i < n; i++) c_mul(A[j][i],tc2[i],tc3[i]);    // multiply A[j]
      else 
        for (i = 0; i < n; i++)
          c_add(tc3[i],c_mul(A[j][i],tc2[i],t),tc3[i]);          // multiply A[j] and accumulate
    }
  }
  sem_post(F+k);                                                 // signal finished
  sem_wait(S+k);                                                 // wait to start
  for (j = 0; j < w; j++){
    if (j%(THREADS) == k+1){
      for (i = 0; i < n; i++) c_cjm(A[j][i],Tc3[0][i],tc2[i]);   // multiply conjugate A[j]
      ifft(n,reverse_cpy(n,tc2,tc1));
      for (i = 0; i < n; i++) OUT(j,i) += tc1[i][0];
    }
  }
  sem_post(F+k);                                                 // signal finished
}

#define ZERO(i) (((i%n)/GX >= YDIM)||(XDIM <= ((i%n)%GX)))

/*

This function performs the multiplication ((1/mu)I - P Z' P^T) h where P is a column of w n*n circulant 
matrice, P^T is a row of w n*n circulant matrices and h is an wn*1 column vector.  

In this algorithm, whenever a row of w n*n circulant matrices is multiplied by a wn*1 column 
vector, a circulant n*n matrix is multiplied by an n*1 piece of a wn*1 column vector w times.  
The resulting w n*1 column vectors are added together. The result is a single n*1 column vector.

In this algorithm, whenever a column of w n*n circulant matrices is multiplied by an n*1 column 
vector, a circulant n*n matrix is multiplied by an n*1 column vector w times.  The result is 
a wn*1 column vector.

Multiplicaton by Z' amounts to ensuring that the resulting column vector has zeros where it should.

See inv_C for an explanation of why it is that C contains the diagonal representation of P^T.

See inv_C for an explanation of why it is that the conjugate of C contains the 
diagonal representation of P.

See multiply_H for an explanation of why it is that the operations below correspond 
to the multiplication of circulant matrices by column vectors.

*/        

void multiply_PTP(double *in, double *out, int n, int w) 
{
  int i,j,k;  double s = 1./(double)n;

  pthread_mutex_lock(&Mutex);
  J = 3;
  pthread_mutex_unlock(&Mutex);
  for(s *= s, i = 0; i < n; i++){
    Tc0[i][0] = s*in[i]; Tc0[i][1] = 0.;                           // complexify and scale
    out[i] = in[i]/MU;
  }
  fft(n,reverse_cpy(n,Tc0,Tc3[0]));                                // transform input
  for (i = 0; i < ((THREADS)-1); i++){
    sem_post(S+i);
  }

  for (j = 0; j < w; j++){
    if (j%(THREADS) == 0){
      for(i = 0; i < n; i++) c_cjm(C[j][i],Tc3[0][i],Tc0[i]);      // multiply conjugate
      ifft(n,reverse_cpy(n,Tc0,Tc1));                              // inverse transform
      for(i = 0; i < n; i++){                                      // multiply Z'
        Tc1[i][1] = 0.;
        if (!ZERO(i)) Tc1[i][0] = 0.;
      }
      fft(n,reverse_cpy(n,Tc1,Tc0));                               // transform input
      for(i = 0; i < n; i++) c_mul(Tc0[i],C[j][i],Tc1[i]);         // multiply matrix
      ifft(n,reverse_cpy(n,Tc1,Tc0));                              // inverse transform
      for(i = 0; i < n; i++) out[i] -= Tc0[i][0];                  // accumulate answer
    }
  }
  for (i = 0; i < ((THREADS)-1); i++){
    sem_wait(F+i);
  }

  for (i = 0; i < n; i++){                                         // global accumulate
    for (k = 1; k < THREADS; k++){
      out[i] -= Tc3[k][i][0];
    }
  }
}

void multiply_PTP_thread(int k, int n, int w, Cx *tc0, Cx *tc1, Cx *tc2, Cx *tc3)
{
  int h,i,j;

  for (h = j = 0; j < w; j++){
    if (j%(THREADS) == k+1){
      for(i = 0; i < n; i++) c_cjm(C[j][i],tc0[i],tc1[i]);       // multiply conjugate
      ifft(n,reverse_cpy(n,tc1,tc2));                            // inverse transform
      for(i = 0; i < n; i++){                                    // multiply Z'
        tc2[i][1] = 0.;
        if (!ZERO(i)) tc2[i][0] = 0.;
      }
      fft(n,reverse_cpy(n,tc2,tc1));                             // transform input
      for(i = 0; i < n; i++) c_mul(tc1[i],C[j][i],tc2[i]);       // multiply matrix
      ifft(n,reverse_cpy(n,tc2,tc1));                            // inverse transform
      if (h++){
        for(i = 0; i < n; i++) tc3[i][0] += tc1[i][0];           // accumulate answer
      } else {
        for(i = 0; i < n; i++) tc3[i][0]  = tc1[i][0];           // copy answer
      }
    }
  }
  sem_post(F+k);
}

/*

This function performs the multiplication P h where P is a column of w n*n circulant 
matrices and h is an n*1 column vector.  In this algorithm, whenever a column of w n*n 
circulant matrices is multiplied by an n*1 column vector, a circulant n*n matrix is 
multiplied by an n*1 column vector w times.  The result is a wn*1 column vector.

See inv_C for an explanation of why it is that the conjugate of C contains the 
diagonal representation of P.

See multiply_H for an explanation of why it is that the operations below correspond 
to the multiplication of circulant matrices by column vectors.

*/

void multiply_P(double *in, double *out, int n, int w) 
{
  int i,j;  double s = 1./(double)n;

  Out = out; 
  pthread_mutex_lock(&Mutex);
  J = 4;
  pthread_mutex_unlock(&Mutex);
  for (i = 0; i < n; i++){
    Tc0[i][0] = s*in[i]; Tc0[i][1] = 0.;                         // complexify and scale
  }
  fft(n,reverse_cpy(n,Tc0,Tc3[0]));                              // transform input
  for (i = 0; i < ((THREADS)-1); i++){
    sem_post(S+i);
  }

  for (j = 0; j < w; j++){                                       // process blocks
    if (j%(THREADS) == 0){
      for (i = 0; i < n; i++) c_cjm(C[j][i],Tc3[0][i],Tc0[i]);   // multiply conjugate
      ifft(n,reverse_cpy(n,Tc0,Tc2));                            // inverse transform
      for (i = 0; i < n; i++){                                   // place answer
        OUT(j,i) = Tc2[i][0];
      }
    }
  }
  for (i = 0; i < ((THREADS)-1); i++){
    sem_wait(F+i);
  }

}

void multiply_P_thread(int k, int n, int w, Cx *tc0, Cx *tc1, Cx *tc2)
{
  int i,j;

  for (j = 0; j < w; j++){                                    // process blocks
    if (j%(THREADS) == k+1){
      for (i = 0; i < n; i++) c_cjm(C[j][i],tc0[i],tc1[i]);   // multiply conjugate
      ifft(n,reverse_cpy(n,tc1,tc2));                         // inverse transform
      for (i = 0; i < n; i++){                                // place answer
        Out[n*j+i] = tc2[i][0];
      }
    }
  }
  sem_post(F+k);
}

/*

This function performs the multiplication P^T h where P^T is a row of w n*n circulant 
matrices and h is an wn*1 column vector.  In this algorithm, whenever a row of w n*n 
circulant matrices is multiplied by a wn*1 column vector, a circulant n*n matrix is 
multiplied by an n*1 piece of a wn*1 column vector w times.  The resulting w n*1 column 
vectors are added together. The result is a single n*1 column vector.

See inv_C for an explanation of why it is that C contains the diagonal representation of P^T.

See multiply_H for an explanation of why it is that the operations below correspond 
to the multiplication of circulant matrices by column vectors.

*/

void multiply_PT(double *in, double *out, int n, int w) 
{
  int i,j,k;  double s = 1./(double)n;  Cx t;

  In = in; 
  pthread_mutex_lock(&Mutex);
  J = 5;
  pthread_mutex_unlock(&Mutex);
  for (i = 0; i < ((THREADS)-1); i++){
    sem_post(S+i);
  }

  for (j = 0; j < w; j++){
    if (j%(THREADS) == 0){
      for (i = 0; i < n; i++){
        Tc1[i][0] = IN(j,i)*s; Tc1[i][1] = 0.;                   // (complexify input)
      }
      fft(n,reverse_cpy(n,Tc1,Tc0));                             // (transform input)
      if (!j)
        for (i = 0; i < n; i++) c_mul(C[j][i],Tc0[i],Tc2[i]);    // (multiply C[j])
      else 
        for (i = 0; i < n; i++)
          c_add(Tc2[i],c_mul(C[j][i],Tc0[i],t),Tc2[i]);          // (multiply C[j] and accumulate)
    }
  }
  ifft(n,reverse_cpy(n,Tc2,Tc3[0]));
  for (i = 0; i < ((THREADS)-1); i++){
    sem_wait(F+i);
  }

  for (i = 0; i < n; i++){
    out[i]  = Tc3[0][i][0];
    for (k = 1; k < THREADS; k++){
      out[i] += Tc3[k][i][0];
    }
  }
}

void multiply_PT_thread(int k, int n, int w, Cx *tc1, Cx *tc2, Cx *tc3)
{
  int h,i,j;  double s = 1./(double)n;  Cx t;

  for (h = j = 0; j < w; j++){
    if (j%(THREADS) == k+1){
      for (i = 0; i < n; i++){
        tc1[i][0] = In[n*j+i]*s; tc1[i][1] = 0.;                 // (complexify input)
      }
      fft(n,reverse_cpy(n,tc1,tc3));                             // (transform input)
      if (!h++)
        for (i = 0; i < n; i++) c_mul(C[j][i],tc3[i],tc2[i]);    // (multiply C[j])
      else 
        for (i = 0; i < n; i++)
          c_add(tc2[i],c_mul(C[j][i],tc3[i],t),tc2[i]);          // (multiply C[j] and accumulate)
    }
  }
  ifft(n,reverse_cpy(n,tc2,tc3));
  sem_post(F+k);
}

/*

The ith row of out is the diagonal of

F* Cbar_i F (I + 1/mu sum ((F C_j F*)(F C^T_j F*)))^{-1/2}

---

The ith row of in is 

sqrt(n) F* C_i e_0

---

C = F* D F

conjugate =>

Cbar = F Dbar F* =>

Dbar = F* Cbar F

conjugate transpose =>

D = F* C^T F

---

D = sqrt(n) diag (F C e_0)

conjugate =>

Dbar = sqrt(n) diag (F* C e_0)

This explains F* Cbar_i F above

It also explains (F C_j F*)(D C^T_j F*) because 

(I + 1/mu sum ((F* Cbar_j F)(F* C^T F))) is real so, conjugate changes nothing =>

(I + 1/mu sum ((F C_j F*)(F C^T F*)))

---

So the ith row of out is the diagonal of

F* Cbar_i F (F(I + 1/mu sum (C_j C^T_j))F*)^{-1/2}

= (F F*) F* Cbar_i F (F F*)(F(I + 1/mu sum (C_j C^T_j))F*)^{-1/2} (F F*)

= F C^T_i {F* (F(I + 1/mu sum (C_j C^T_j))F*)^{-1/2} F} F* <-- from below

---

C = F* D F

conjugate transpose =>

C^T = F* Dbar F <-- C is real, D is symmetric
 
= F* F* Cbar F F

---

Check:

(F* (F(I + 1/mu sum (C_j C^T_j))F*)^{-1/2} F)^2 (I + 1/mu sum(C_j C^T_j)) =? I

F* (..)^{-1/2} F F* (..)^{-1/2} F

= F* (..)^{-1} F

= F* (F (I + 1/mu sum (C_j C^T_j))F*)^{-1} F

= F* F (I + 1/mu sum (C_j C^T_j))^{-1} F* F  <-- (A_0 .. A_n)^{-1} = A_n^{-1} .. A_0^{-1}

= (I + 1/mu sum (C_j C^T_j))^{-1}

---

Therefore the ith row of out is the diagonal of 

F C^T_i (I + 1/mu sum (C_j C^T_j))^{-1/2} F*

___

Later, when the above is used as the diagonal representation of P, it is first conjugated =>

F* C^T_i (I + 1/mu sum (C_j C^T_j))^{-1/2} F

The above is precisely the diagonal representation of P as needed by this method.

In other words it is Dbar in the expression C = F Dbar F* where C is the circulant matrix P.

When the above is used as the diagonal representation of P^T, it is without conjugation.  This is because 

C^T = F* Dbar F 

then conjugate both sides =>
 
C^T = F D F* <-- C^T is real

*/
double inv_C(Cx **in, double mu, Cx **out, int n, int w)
{
  int i,j;

  for (j = 0; j < w; j++){
    if (!j)
      for (i = 0; i < n; i++){
        Td4[i] = 1. + c_msq(in[j][i])/mu;                    
        c_set(in[j][i],out[j][i]);                           
      }
    else 
      for (i = 0; i < n; i++){
        Td4[i] += c_msq(in[j][i])/mu;
        c_set(in[j][i],out[j][i]);
      }
  }
  for (i = 0; i < n; i++) Td4[i] = 1./(sqrt(Td4[i])*mu);     
  for (j = 0; j < w; j++)  
    for (i = 0; i < n; i++){
      out[j][i][0] *= Td4[i];
      out[j][i][1] *= Td4[i];
    }
  return 1./mu;
}

#if 0
void show_P(int n, int w)
{
  int i,j,k;  double s = 1./n;

  for (k = 0; k < w; k++){
    ifft(n,reverse_cpy(n,C[k],Tc0));
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
        Pm[k][j][i] = s*Tc0[(n+i-j)%n][0];
  }
  for (k = 0; k < w; k++){
    for (i = 0; i < n; i++){
      for (j = 0; j < n; j++){
        //      printf("%f ",Pm[k][i][j]);
      }
      putchar('\n');
    }
    putchar('\n');    
  }
}

void mul_P(double *in, double *out, int n, int w)
{
  int i,j;  double s;

  for(i = 0; i < n; i++){
    Tc0[i][0] = in[i];
  }
  for(i = 0; i < w*n; i++){
    for(s = j = 0; j < n; j++)
      s += Pm[i/n][i%n][j]*Tc0[j][0];
    out[i] = s;
  }
}
#endif

/*

This function computes Mh where M = (1/mu)I - PP^T, P is a column of w n*n circulant 
matrice, P^T is a row of w n*n circulant matrices, and h is a wn*1 column vector.  

In this algorithm, whenever a row of w n*n circulant matrices is multiplied by a wn*1 column 
vector, a circulant n*n matrix is multiplied by an n*1 piece of a wn*1 column vector w times.  
The resulting w n*1 column vectors are added together. The result is a single n*1 column vector.

In this algorithm, whenever a column of w n*n circulant matrices is multiplied by an n*1 column 
vector, a circulant n*n matrix is multiplied by an n*1 column vector w times.  The result is 
a wn*1 column vector.

See inv_C for an explanation of why it is that C contains the diagonal representation of P^T.

See inv_C for an explanation of why it is that the conjugate of C contains the 
diagonal representation of P.

See multiply_H for an explanation of why it is that the operations below correspond 
to the multiplication of circulant matrices by column vectors.

---

It is also true that this function is used to compute PP^T u.

*/

void mul(int f, Cx **m, double mu, double *in, double *out, int n, int w)
{
  int i,j,k;  double s = 1./(double)n;  Cx t;

  In = in; Out = out; E = mu;
  for (i = 0; i < ((THREADS)-1); i++){
    sem_post(S+i);
  }

  //printf("after sem_post\n");

  for (j = 0; j < w; j++){
    if (j%(THREADS) == 0){
      for (i = 0; i < n; i++){
        Tc1[i][0] = IN(j,i)*s; Tc1[i][1] = 0.;                   // complexify input
        OUT(j,i) = mu*IN(j,i);                                   // mu*in
      }
      fft(n,reverse_cpy(n,Tc1,Tc0));                             // transform input

      //printf("after fft in mul\n");

      if (!j)
        for (i = 0; i < n; i++) c_mul(m[j][i],Tc0[i],Tc3[0][i]); // multiply m[j]
      else 
        for (i = 0; i < n; i++)
          c_add(Tc3[0][i],c_mul(m[j][i],Tc0[i],t),Tc3[0][i]);    // multiply m[j] and accumulate
    }
  }
  for (i = 0; i < ((THREADS)-1); i++){
    sem_wait(F+i);
  }

  //printf("after sem_wait\n");

  for (i = 0; i < n; i++){                                       // global accumulate
    for (k = 1; k < THREADS; k++){
      c_add(Tc3[0][i],Tc3[k][i],Tc3[0][i]);
    }
  }
  for (i = 0; i < ((THREADS)-1); i++){
    sem_post(S+i);
  }

  //printf("after sem_post\n");

  //sem_post(S+0); sem_post(S+1); sem_post(S+2);
  for (j = 0; j < w; j++){
    if (j%(THREADS) == 0){
      for (i = 0; i < n; i++) c_cjm(m[j][i],Tc3[0][i],Tc2[i]);   // multiply conjugate m[j]
      ifft(n,reverse_cpy(n,Tc2,Tc1));
      if (f)
        for (i = 0; i < n; i++) OUT(j,i) += Tc1[i][0];
      else
        for (i = 0; i < n; i++) OUT(j,i) -= Tc1[i][0];
    }
  }
  for (i = 0; i < ((THREADS)-1); i++){
    sem_wait(F+i);
  }

  //printf("after sem_wait\n");

  //sem_wait(F+0); sem_wait(F+1); sem_wait(F+2);                   // wait for threads
}

void m_t(int f, int k, int n, int w, Cx *t1, Cx *t2, Cx *t3, double mu, double *in, double *out)
{
  int h,i,j;  double s = 1./(double)n;  Cx t;

  //printf("%d in m_t\n",k);
  for (h = j = 0; j < w; j++){
    //printf("%d in for loop; j = %d\n",k,j);
    if (j%(THREADS) == k+1){
      for (i = 0; i < n; i++){
        t1[i][0] = IN(j,i)*s; t1[i][1] = 0.;                     // complexify input
        OUT(j,i) = mu*IN(j,i);                                   // mu*in
      }
      //printf("%d calling fft (%d)\n",k,n);
      fft(n,reverse_cpy(n,t1,t2));                               // transform input
      //printf("%d back from fft\n",k);
      if (!h++)
        for (i = 0; i < n; i++) c_mul(C[j][i],t2[i],t3[i]);      // multiply C[j]
      else
        for (i = 0; i < n; i++)
          c_add(t3[i],c_mul(C[j][i],t2[i],t),t3[i]);             // multiply C[j] and accumulate
    }
  }
  //printf("before %d posting F\n",k);
  sem_post(F+k);                                                 // signal finished
  //printf("after %d posting F\n",k);
  sem_wait(S+k);                                                 // wait to start
  //printf("after %d waiting S\n",k);
  for (j = 0; j < w; j++){
    if (j%(THREADS) == k+1){
      for (i = 0; i < n; i++) c_cjm(C[j][i],Tc3[0][i],t2[i]);    // multiply conjugate C[j]
      ifft(n,reverse_cpy(n,t2,t1));
      if (f)
        for (i = 0; i < n; i++) OUT(j,i) += t1[i][0];
      else
        for (i = 0; i < n; i++) OUT(j,i) -= t1[i][0];
    }
  }   
  sem_post(F+k);                                                 // signal finished
}

void *multiply_thread(void *a)
{
  long int k = (long int)a;  int n = GX_DIM*GY_DIM, w = WVDIM;  Cx *tc1,*tc2,*tc3;  int j;

  tc1            = malloc(n*sizeof(Cx));
  tc2            = malloc(n*sizeof(Cx));
  Tc3[k+1] = tc3 = malloc(n*sizeof(Cx));
 t:
  sem_wait(S+k);     

  //printf("%ld\n",k);

  pthread_mutex_lock(&Mutex);                                    // wait to start
  j = J;
  switch(J){
  case 0:
    pthread_mutex_unlock(&Mutex);
    //printf("%ld starting %d\n",k,j);
    multiply_A_thread(k,n,w,tc1,tc2,tc3,(double *)In,(double *)Out);
    break;
  case 1:
    pthread_mutex_unlock(&Mutex);
    //printf("%ld starting %d\n",k,j);
    m_t(0,k,n,w,tc1,tc2,tc3,E,(double *)In,(double *)Out);       // mul(0,C,e,in,out,n,w);
    break;
  case 2:
    pthread_mutex_unlock(&Mutex);
    //printf("%ld starting %d\n",k,j);
    m_t(1,k,n,w,tc1,tc2,tc3,0.,(double *)In,(double *)Out);      // mul(1,C,0,in,out,n,w);
    break;
  case 3:
    pthread_mutex_unlock(&Mutex);
    //printf("%ld starting %d\n",k,j);
    multiply_PTP_thread(k,n,w,Tc3[0],tc1,tc2,tc3);
    break;
  case 4:
    pthread_mutex_unlock(&Mutex);
    //printf("%ld starting %d\n",k,j);
    multiply_P_thread(k,n,w,Tc3[0],tc1,tc2);
    break;
  case 5:
    pthread_mutex_unlock(&Mutex);
    //printf("%ld starting %d\n",k,j);
    multiply_PT_thread(k,n,w,tc1,tc2,tc3);
    break;
  case 6:
    pthread_mutex_unlock(&Mutex);
    //printf("%ld starting %d\n",k,j);
    free(tc1); free(tc2); free(tc3); 
    pthread_exit(NULL);
  }
  //printf("%ld finished with %d\n",k,j);
  goto t;
  return NULL;
}

/*

This function does 2 things:

1) Either compute P^T and store it to disk or read it from disk 

2) Compute Mh from paper where M = (1/mu)I - PP^T

*/

void multiply_C(double *in, double *out, int n, int w) 
{
  char b[1<<14]; int i,k; static int e = 1;  FILE *f;

  // the first time this function is called, either compute P^T and store in on disk or read it from 
  // disk depending on how #if 1 is set
  if (e){

// ###
#if COMPUTE_PT                          
    inv_C(A,MU,C,n,w);                                           // 
    for (e = i = 0; i < w; i++){
      sprintf(b,"C_%d",i);
      if ((f = fopen(b,"w"))){
        k=write(fileno(f),C[i],n*sizeof(Cx)); // "diagonal" representation of circulant
        fclose(f);
      } else {
        error(1,"could not write %s\n",b);
      }
    }
#else                             /* read precomputed C */
    for (e = i = 0; i < w; i++){
      sprintf(b,"C_%d",i);
      if ((f = fopen(b,"r"))){
        k=read(fileno(f),C[i],n*sizeof(Cx));
        fclose(f);
      } else {
        error(1,"could not open %s\n",b);
      }
    }
#endif
  }
  pthread_mutex_lock(&Mutex);
  J = 1;
  pthread_mutex_unlock(&Mutex);

  //printf("eiei\n");
                                                                 // compute Mh from paper where M = (1/mu)I - PP^T
  mul(0,C,1./MU,in,out,n,w);

  //printf("eiei\n");
   
}

void multiply_PPT(double *in, double *out, int n, int w) 
{
  pthread_mutex_lock(&Mutex);
  J = 2;
  pthread_mutex_unlock(&Mutex);
   
  mul(1,C,0.,in,out,n,w);                                        // compute PP^T u
}

/*

This function reads f off the hard drive.  

It will be multiplied by H on the left to produce a 'synthetic' g = Hf

*/

/*

void pgm2f(double *f, char *name)
{
  char s[1<<14]; int i,k; FILE *file;

  //if ((WVDIM != 24)||(XDIM != 81)||(YDIM != 90)) error(1,"inappropriate geometry\n");
  if (!(file = fopen(name,"r"))) error(1,"could not open %s\n",s);
  for (i = 3; i--; fgets(s,1<<14,file));                         // header
  for (i = 0; i < WVDIM*XDIM*YDIM; i++){                         // problem size
    fscanf(file,"%d",&k);                                        // read one value
    f[i] = (double)k;                                            // populate f
  }    
  fclose(file);
}

*/

/* 

This function stores the reconstructed f to a pgm file such that each WVDIM-th of the image 
is a subimage corresponding a single wavelength

*/

void f2pgm(double *f, char *name)
{

  int h,i,j,k,m; FILE *file;

  if (!(file = fopen(name,"w"))) error(1,"could not open %s\n",name);
  // fprintf(file,"P3\n# %s\n%d %d\n255\n",name,4*XDIM,2*YDIM);
  // fprintf(file,"P2\n%d %d\n255\n",6*XDIM,4*YDIM);
  fprintf(file,"P2\n%d %d\n255\n",XDIM,YDIM);




// for (m = 0; m < 4; m++){  

  // for (k = 0; k < YDIM; k++){

    // for (j = 0; j < 6; j++){

      // for (i = 0; i < XDIM; i++){
      for (i = 0; i < XDIM*YDIM; i++){
        // if ((h = 0.5+f[j*XDIM*YDIM+k*XDIM+m*6*XDIM*YDIM+i]) < 0) h = 0;
        // if ((h = 0.5+f[j*XDIM*YDIM+k*XDIM+m*6*XDIM*YDIM+i]) < 0) h = 0;
        // if (h > 255) h = 255;
        // fprintf(file,"%d\n",h);
        fprintf(file,"%lf\n",f[i]);
      }

    // }

  // }
// }



  fclose(file);
}

void g2ppm(double *g)
{
  int h,i,j; FILE *file;  double v,m = 0;

  if (!(file = fopen("focal.ppm","w"))) error(1,"could not open focal\n");
  fprintf(file,"P2\n# %s\n%d %d\n255\n","focal.ppm",GX_DIM,GY_DIM);
  for (i = 0; i < GX_DIM*GY_DIM; i++){
    v = g[i]*3.80649075e-06;
    if (v > m) j = i, m = v;
    if ((h = 0.5+v) < 0) h = 0;
    if (h > 255) h = 255;
    fprintf(file,"%d\n",h);
  }
  fclose(file);
  printf("%d : max = %lf\n",j,m);
}

void f2ppm(double *f, char *name)
{
  int h,i; FILE *file;

  if ((WVDIM != 24)||(XDIM != 81)||(YDIM != 90)) error(1,"inappropriate geometry\n");
  if (!(file = fopen(name,"w"))) error(1,"could not open %s\n",name);
  fprintf(file,"P3\n# %s\n%d %d\n255\n",name,4*XDIM,2*YDIM);
  for (i = 0; i < WVDIM*XDIM*YDIM; i++){
    if ((h = 0.5+f[i]) < 0) h = 0;
    if (h > 255) h = 255;
    fprintf(file,"%d\n",h);
  }
  fclose(file);
}


/*
void f2pgm(double *f, char *dir)
{
  char s[1<<14]; int h,i,j,k; FILE *file;

  mkdir(dir,0777);
  for (k = 0; k < WVDIM; k++) {
    sprintf(s,"%s/f%d.pgm",dir,k+1);
    file = fopen(s,"w");
    fprintf(file,"P2\n# %s\n%d %d\n255\n",s,XDIM,YDIM);
    for (i = 0; i < YDIM; i++) {
      for (j = 0; j < XDIM; j++) {
        if ((h = 0.5+f[k*XDIM*YDIM+i*XDIM+j]) < 0) h = 0;
        if (h > 255) h = 255;
        fprintf(file,"%d\n",h);
      }
    }
    fclose(file);
  }
}
*/

int solve(double **z, double *y, int *p, int n)
  /*
    returns 1 if a invertible and then y[p[]] is solution to a[][]x[] = y[]
    returns 0 otherwise
  */
{
  int i,j,k,*q;  double m,t,**a;

  if (n == 1){
    p[0]  = 0;
    y[0] /= z[0][0];
    return (fabs(z[0][0]) > 1e-10);
  }
  q = malloc(n*sizeof(int));
  a = malloc(n*sizeof(double *));
  for (i = 0; i < n; i++) a[i] = malloc(n*sizeof(double));
  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++)
      a[i][j] = z[i][j];
  }
  for (i = 0; i < n; q[i] = i, i++);
  for (i = 0; i < n-1; i++){
    for (m = 0, k = i; k < n; k++)
      if ((t = fabs(a[i][k])) > m) m = t, j = k;                 // find pivot
    if (m < 1e-10){
      free(q);
      for (i = 0; i < n; i++) free(a[i]); free(a);
      return 0;                                                  // small determinant?
    }
    if (j > i){                                                  // swap
      k = q[j], q[j] = q[i], q[i] = k;
      for (k = 0; k < n; k++) t = a[k][j], a[k][j] = a[k][i], a[k][i] = t;
    }
    for (k = i+1; k < n; k++)
      if (0 != a[k][i]){
        for (m = a[k][i]/a[i][i], j = i+1; j < n; j++) a[k][j] -= m*a[i][j];
        y[k] -= m*y[i];
      }
  }
  for (i = n; i--; ){
    for (j = i+1; j < n; j++) y[i] -= a[i][j]*y[j];
    y[i] /= a[i][i];
    p[q[i]] = i;
  }
  free(q);
  for (i = 0; i < n; i++) free(a[i]); free(a);
  return 1;
}

double get_time() 
{
  struct timeval tv;  struct timezone tz;
  
  gettimeofday(&tv,&tz);
  return tv.tv_sec + tv.tv_usec/1000000.;
} 

/* 

This function reports the error at a given cg step.  It's input is:

 e: current cg 'error'
it: cg iteration number

It reports:

current cg 'error'
cg iteration number

*/

void show_error(double e, int it)
{
  printf("%3d : conjugate gradient error = %e\n",it,e);
}

#define ISNORMAL(x) (finite(x)&&(x != 0))

/*

This function adds poisson distributed noise with mean and variance Hf_i and 
gaussian noise with zero mean and variance sigma_s to g.

It enforces the condition worked out in the paper that 
g_i be >= sqrt(sigma_s^4 + sigma_s^2) - sigma_s^2.

*/

void add_noise(double *g, double s, int n)
{
  int i; 

  //double b = sqrt(pow(s,4.)+pow(s,2.))-pow(s,2.);

  for (i = 0; i < n; i++){
    g[i] = poisson(g[i]) + normal(0.,s);                         // poisson noise, gaussian (normal) noise <--- random.c
    //if (g[i] < b) g[i] = b;                                    // condition (3.5) from paper
  }
}

/*

This function writes g to disk.

There will be one float per line.

*/

void g2disk (double *g)
{
  int i; FILE *file;

  if (!(file = fopen(FOCAL_PLANE,"w"))) error(1,"could not open %s\n",FOCAL_PLANE);



      for (i = 0; i < GX_DIM*GY_DIM; i++){                       // problem size
        fprintf(file,"%f\n",g[i]);                               // print value
      }


	  
  fclose(file);
}

void disk2g (double *g, char *name)
{
  int i,k; FILE *file;

  if (!(file = fopen(name,"r"))) error(1,"could not open %s\n",name);



      for (i = 0; i < GX_DIM*GY_DIM; i++){                       // problem size
        k=fscanf(file,"%lf",&g[i]);                                 // print value
      }


	  
  fclose(file);
}

/*

This function 'removes noise'.  It produces the most probable noise free focal plane image x = Hf.

It does so by following the recipe given by (3.4) in the paper.

*/

void remove_noise(double *g, double s, int n)
{
  int i; double b = sqrt(pow(s,4.)+pow(s,2.))-pow(s,2.);

  for (i = 0; i < n; i++)
    if (g[i] < b) g[i] = b;                                      // condition (3.5) from paper
    g[i] = (sqrt(pow(1.+2.*pow(s,2.),2.)+4*(g[i]*(g[i]+2.*pow(s,2.))-pow(s,2.)))-(1.+2.*pow(s,2.)))/2.;
}

int main(int argc, char **argv)
{
  char buf[256]; int w,n,m,i,j,k,l,itr; double t0,t1,a,b,c,q,e,*f,*g,*h,*o,*p,*x,*y,*r,*s,*t,*u,*v,*z,*ts; 
  pthread_t thread[(THREADS)-1];
  signal(SIGINT,quit);

  printf("MU = %e\n",MU);

  R   = RESTART;                                // # of vose-horton steps
  Lim = LIM;                                    // # of cg steps
  itr = 0;                                      // current cg step

  w = WVDIM;                                    // # wavelengths spectrometer calibrated at
  n = GX_DIM*GY_DIM;                            // # rows in H
  m = w * XDIM*YDIM;                            // # columns in H

  init(n,m,w);                                  // malloc space; read diagonalized representation of system matrix off disk; read fftw wisdom off disk

printf("dkdk\n");  

  for (i = 0; i < ((THREADS)-1); i++)             // n threads
    pthread_create(&thread[i],NULL,multiply_thread,(void *)(long int)i);

  f = calloc(m,sizeof(double));                 // image

f[0]=255.0;

  h = malloc(sizeof(double)*n);                 // x = H f

  y  = calloc(w*n,sizeof(double));              // answer
  u  = malloc(sizeof(double)*w*n);              // residual
  s  = malloc(sizeof(double)*w*n);              // h = script_H^T x
  t  = malloc(sizeof(double)*w*n);
  z  = malloc(sizeof(double)*w*n);              // Mh

  multiply_C(u,z,n,w);                          // compute P^T and store on disk

printf("fjfj\n");  

  g = malloc(sizeof(double)*n);
  o = malloc(sizeof(double)*n);
  p = malloc(sizeof(double)*n);
  r = malloc(sizeof(double)*n);
  v = malloc(sizeof(double)*n);
  x = malloc(sizeof(double)*n);

#if GENERATE_G
  //pgm2f(f,SYNTHETIC_F);

f2pgm(f,PICTURE_OUT);
//exit(0);

printf("eiei\n");  

  multiply_H(f,h,n,m,w);                        // h = Hf

  //add_noise(h,2.,n);                            // add poisson distributed noise, normal distributed noise
  g2disk(h);                                    // write g to disk
  exit(0);
#else
  disk2g(h,FOCAL_PLANE);                        // get g from disk
#endif

printf("ghgh\n");  

#if REMOVE_NOISE                                
  remove_noise(h,2.,n);                         // remove noise via (3.4) from paper; enforce condition (3.5)
#endif

  t1 = 0.;
  t0 = get_time();                              // begin marking time
  printf("Start...\n");

  embed(h,s,n,w);                               // h = script_H^T x
  memcpy(u,s,sizeof(double)*w*n);               // u = residual
  j = 0;                                        // vose-horton step number 

  do {                                          // R = RESTART # of vose-horton steps

    multiply_C(u,z,n,w);                        // Mh <--- M = (1/mu) I - PP^T

    for (i = 0; i < w*n; i++){                  // multiply by Z'
      if (ZERO(i))
        t[i] = z[i];
      else
        t[i] = 0.;
    }

    multiply_PT(t,g,n,w);                       // multiply by P^T

    for (a = itr = i = 0; i < n; i++){          // set up conjugate gradient method to solve v = ((1/muu)I - P^T Z' P) u
      o[i] = x[i] = 0.;                         // old guess is o, new is x
      p[i] = r[i] = g[i];                       // g replaces b - Ax from paper
      a += r[i]*r[i];                           // current error
    }
    c = a;

    while ((e = a) > EPSILON){                  // stop if error below threshold
      printf("cg...\n");
      multiply_PTP(p,v,n,w);                    // matrix: (1/mu)I - P^T Z' P
      for (q = i = 0; i < n; i++){
        q += p[i]*v[i]; 
      }
      a = e/q;                                  // a = e/(p,v)
      if (!ISNORMAL(a)) break;                  // break if a is zero, infinite, or nan
      for (b = i = 0; i < n; i++){
        x[i]  = o[i] + a*p[i];                  // x = x + a*p
        r[i] -= a*v[i];                         // r = r - a*v
        b    += r[i]*r[i];                      // a = (r,r)
      }
      if ((b < c)||(Lim == LIM)){               // at least one step; stop at local minimum of error
        c = b;
        ts = x; x = o; o = ts;                  // exchange old guess with new

#if !REPORT
        if (--Lim < 1) break;                   // stop at LIM number of steps
#endif
      } else {
        break;
      }
      a = b;
      if (!ISNORMAL(a)) break;
      b /= e;
      if (!ISNORMAL(b)) break;
      for (i = 0; i < n; i++){
        p[i] = r[i] + b*p[i];                   // p = r + (a/e)*p
      }


#if REPORT                                      // report error for this conjugate gradient step 
      t1 += get_time() - t0;                    // collect time so far; stop recording time

      show_error(a,itr++);                      // report error: cg error, 'f', f, size of f, cg iteration #

      t0 = get_time();                          // reporting over; begin marking time
#else
      printf("%2d : %e\n",itr++,a);
#endif
#if REPORT
      if (--Lim < 1) break;                     // LIM number of cg steps
#endif
    }
    printf("restart...\n");

    multiply_P(o,t,n,w);                        // multiply by P

    for (i = 0; i < w*n; i++){                  // add to Mh; multiply by Z' on the left
      if (ZERO(i))
        t[i] += z[i];
      else
        t[i] = 0;
    }

#define IU(i) ((int *)u)[i]

    multiply_PPT(t,Vec,n,w);                    // multiply by P P^T

    for (k = i = 0; i < w*n; i++)               // perform the computation Z(u + MU v); next suck f out of y_infinity
      if (ZERO(i)){
        Vec[i] = 0.;
      } else {
        IU(k)  = i;                             // corresponding nonzero index
        t[k++] = Vec[i] = z[i] + MU*Vec[i];     // t is condensed version of Vec; i.e., f
      }

                                                // store the most recently computed f to disk
    sprintf(buf,"Vec%d",j); swap_out(buf,t,k*sizeof(double));

    multiply_A(Vec,z,n,w);                      // compute script_H^T script_H y_infinity
    for (M[j][j] = B[j] = i = 0; i < k; i++){
      t[i]     = z[IU(i)];                      // t is condensed version of z = A Vec; i.e., H^T H f
      B[j]    += Vec[IU(i)]*s[IU(i)];           // V^T x
      M[j][j] += Vec[IU(i)]*z[IU(i)];           // part of V^T V
    }

                                                // store the most recently computed H^T T f to disk
    sprintf(buf,"AVec%d",j); swap_out(buf,t,k*sizeof(double));

    
    for (l = j; l--; ){                         // compute the rest of V^T V
      sprintf(buf,"AVec%d",l); swap_in(buf,z,k*sizeof(double));
      for (M[j][l] = i = 0; i < k; i++) M[j][l] += Vec[IU(i)]*z[i];
    }
    for (i = 0; i < j; i++) M[i][j] = M[j][i];
    memcpy(Alpha,B,(j+1)*sizeof(double));

                                                // solve for eta = (V^T V)^{-1} V^T x 
    if (!solve(M,Alpha,(int *)t,j+1)) goto end;

                                                // compute sum eta_i y_infinity_i
    for (i = 0; i < w*n; i++) y[i] = Alpha[((int *)t)[j]]*Vec[i];
    for (l = j; l--; ){
      sprintf(buf,"Vec%d",l); swap_in(buf,Vec,k*sizeof(double)); 
      for (i = 0; i < k; i++)
        y[IU(i)] += Alpha[((int *)t)[l]]*Vec[i];
    }
#undef IU
    if ((R < 1)||(Lim == LIM)) break;           // 

    multiply_A(y,t,n,w);                        // compute script_H^T script_H 'y_infinity'

    for (i = 0; i < w*n; i++){                  // compute 'script_H^T (x - script_H sum eta_i y_infinity_i)'
      u[i] = s[i] - t[i];
    }
    Lim = LIM; j++;                             // reset lim
  } while (--R > 0);                            // R number of vose-horton steps
 end:

  for (k = i = 0; i < n*w; i++){                // suck f out of the final answer
    if (ZERO(i)) continue;
    Td4[k++] = y[i];
  }
  t1 += get_time() - t0;                        // compute time
  printf("time = %lf\n",t1);                    // print time
#if PICTURE
  f2pgm(Td4,PICTURE_OUT);                       // produce final ppm image
#endif
  if (k != m) error(1,"extraction error: k = %d, m = %d\n",k,m);
  putchar('\n');

  multiply_H(Td4,x,n,m,w);                      // multiply reconstructed image by H on the left
  for (e = q = i = 0; i < n; i++){

    if (ZERO(i)) continue;                      // ### ###
    b = fabs(x[i]-h[i]);                        // absolute difference between original focal plane and reconstructed focal plane
    if (e < b) e = b;                           // max absolute difference
    if (fabs(h[i]) > 1.) b /= fabs(h[i]);       // absolute difference / original
    if (q < b) q = b;                           // max absolute difference / original
  }

                                                // max absolute difference (max absolute difference / original)
  printf("error       = %e (%e %c)\n",e,q*100,'%');
  pthread_mutex_lock(&Mutex);
  J = 6; 
  pthread_mutex_unlock(&Mutex);
  
  for (i = 0; i < ((THREADS)-1); i++){
    sem_post(S+i);
  }
  
  for (i = 0; i < ((THREADS)-1); i++){
    pthread_join(thread[i],NULL);
  }
  
  
  free(f);
  free(h);
  free(y);
  free(u);
  free(s);
  free(t);
  free(z);
  free(g);
  free(o);
  free(p);
  free(r);
  free(v);
  free(x);

  free(Vec);
  free(Tc0);
  free(Tc1);
  free(Tc2);
  free(Tc3[0]);
  free(Td4);

#if NEW
  free(B);
  free(Alpha);
  for (i = 0; i < RESTART; i++) free(M[i]);
  free(M);
#endif

  for(i = 0; i < w; i++){
    free(C[i]);
    free(A[i]);
  }
  free(A);
  free(C);
  return 0;
}

