#define  TWO_32     (4294967296.0)
#define  TWO_LONGSZ TWO_32

typedef unsigned Ulong;

Ulong  rtab[55];
int    rndx,nflg;
double nv;

#define rndm   ((++rndx>54)?rtab[rndx=nrndm()]:rtab[rndx]) /* random 32-bit generator  */
#define U01    (rndm/TWO_LONGSZ)                           /* random in interval [0,1) */
#define U(x)   (U01*(x))                                   /* random in interval [0,x) */
#define rnd(n) ((Ulong)U(n))                               /* random from set {0..n-1} */

typedef struct dist
{ /* for generating random number i with probability p[i] */
  double *p;
  Ulong  *a;
  Ulong   n;
  double m1;
  double m2;
} distribution;

void sv_rnd(Ulong *v);                             // save state; also saves for normal, thus v[59]
void rst_rnd(Ulong *v);                            // restore state; also restores for normal, thus v[59]
int nrndm();
void initrand(Ulong j);                            // Initialize the 32 bit random number generator with seed j
distribution *allocdist(long n);                   // allocate a distribution over {0..n-1}
distribution *initdist(distribution *d, double s); // Initialize the distribution d
void freedist( distribution *d );                  // frees dist
Ulong drand(distribution *d);                      // Return element from {0..d->n-1} according to d
double normal(double m, double s);
