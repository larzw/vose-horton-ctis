#ifndef COMPLEX_H
#define COMPLEX_H

#define PI 3.14159265358979323846

typedef double Cx[2];

double * c_cnj (Cx  , Cx        );
double * c_inv (Cx  , Cx        );
double * c_add (Cx  , Cx  , Cx);
double * c_cjm (Cx  , Cx  , Cx);
double * c_sub (Cx  , Cx  , Cx);
double * c_mul (Cx  , Cx  , Cx);
double * c_div (Cx  , Cx  , Cx);
double * c_e   (double, Cx        );
double * c_pow (Cx  , double, Cx);
double * c_set (Cx  , Cx        );
double c_msq (Cx                );

#endif


