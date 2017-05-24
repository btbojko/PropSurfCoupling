/* ROCSTAT - random.c */

#include <math.h>
#include <stdlib.h>

double random_uniform(double a, double b){
	return drand48()*(b-a)+a;
}

void random_unit_vector(double *v){
	double z = random_uniform(-1.0, 1.0);
	double r = sqrt(1.0 - z*z);
	double t = random_uniform(0.0, 2*M_PI);

	v[0] = r * cos(t);
	v[1] = r * sin(t);
	v[2] = z;
}

void random_normal(double *n, double *v){
	double norm, x[3];

	do{
		random_unit_vector(x);
		n[0] = x[1]*v[2] - x[2]*v[1];
		n[1] = x[2]*v[0] - x[0]*v[2];
		n[2] = x[0]*v[1] - x[1]*v[0];
		norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
	} while(norm < 0.01);

	n[0] /= norm; n[1] /= norm; n[2] /= norm;
}

double random_uniform_r(unsigned short xsubi[3], double a, double b){
	return erand48(xsubi)*(b-a)+a;
}

void random_unit_vector_r(unsigned short xsubi[3], double *v){
	double z = random_uniform_r(xsubi, -1.0, 1.0);
	double r = sqrt(1.0 - z*z);
	double t = random_uniform_r(xsubi, 0.0, 2*M_PI);

	v[0] = r * cos(t);
	v[1] = r * sin(t);
	v[2] = z;
}

void random_normal_r(unsigned short xsubi[3], double *n, double *v){
	double norm, x[3];

	do{
		random_unit_vector_r(xsubi, x);
		n[0] = x[1]*v[2] - x[2]*v[1];
		n[1] = x[2]*v[0] - x[0]*v[2];
		n[2] = x[0]*v[1] - x[1]*v[0];
		norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
	} while(norm < 0.01);

	n[0] /= norm; n[1] /= norm; n[2] /= norm;
}

