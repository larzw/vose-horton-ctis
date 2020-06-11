#include <math.h>
#include <fftw3.h>
#include "complex.h"

Cx c_ZERO = {0,0}, c_ONE = {1,0};

double *c_pow(Cx a, double p, Cx b)
{
  double r, theta;
  
  r = sqrt((a[0])*(a[0])+(a[1])*(a[1]));
  if (a[0] == 0.) {
    if (a[1] > 0.) {
      theta = M_PI/2.;
    } else {
      theta = 3. * M_PI / 2.;
    }
  } else {
    theta = atan((a[1])/(a[0]));
  }
  if ((a[0] < 0.) && (a[1] < 0.)) {
    theta = theta - M_PI;
  }
  if ((theta < 0.) && (a[0] < 0.) && (a[1] > 0.)) {
    theta = theta + M_PI;
  }

  r = pow(r, p);
  theta = theta * p;

  b[0] = r*cos(theta);
  b[1] = r*sin(theta);
  return b;
}  

double c_msq(Cx a)
{
  return a[0]*a[0] + a[1]*a[1];
}

double *c_cnj(Cx a, Cx r)
{
  r[0] =  a[0];
  r[1] = -a[1];
  return r;
}

double *c_add(Cx a, Cx b, Cx r)
{
  r[0] = a[0] + b[0];
  r[1] = a[1] + b[1];
  return r;
}

double *c_sub(Cx a, Cx b, Cx r)
{
  r[0] = a[0] - b[0];
  r[1] = a[1] - b[1];
  return r;
}

double *c_mul(Cx a, Cx b, Cx r)
{
  r[0] = a[0]*b[0] - a[1]*b[1];
  r[1] = a[0]*b[1] + b[0]*a[1];
  return r;
}

double *c_cjm(Cx a, Cx b, Cx r)
{
  r[0] = a[0]*b[0] + a[1]*b[1];
  r[1] = a[0]*b[1] - b[0]*a[1];
  return r;
}

double *c_div(Cx a, Cx b, Cx r)
{
  Cx c; double d = c_msq(b);

  c_mul(a,c_cnj(b,c),r);
  r[0] /= d;
  r[1] /= d;
  return r;
}

double *c_inv(Cx a, Cx r)
{
  Cx c = {1.,0.};
  return c_div(c,a,r);
}

double *c_set(Cx a, Cx r)
{
  r[0] = a[0]; r[1] = a[1];
  return r;
}

double *c_e(double x, Cx r)
{
  r[0] = cos(2.0*M_PI*x); r[1] = sin(2.0*M_PI*x);
  return r;
}
