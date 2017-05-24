/* ROCSTAT - random.h */

#ifndef _RANDOM_H
#define _RANDOM_H

double random_uniform(double a, double b);
void random_unit_vector(double *v);
void random_normal(double *n, double *v);

double random_uniform_r(unsigned short xsubi[3], double a, double b);
void random_unit_vector_r(unsigned short xsubi[3], double *v);
void random_normal_r(unsigned short xsubi[3], double *n, double *v);

#endif
