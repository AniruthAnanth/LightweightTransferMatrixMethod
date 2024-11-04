

#include <stdio.h>
#include <math.h>
#include <complex.h>

const long double PI = 3.1415;

double angle(complex double z) {
	double angle = cimagl(z)/creall(z);
	return angle;
}

complex double complex_modulo(double complex z1, double complex z2) {
	double complex remainder = abs(z1) % abs(z2);

	return remainder * cexp(I*angle(z1));
}

complex double sin_complex(int v) {
	double complex a = fmod(creall(v), (2*PI));
	double complex b = complex_modulo(cimagl(v),2*PI*I);

	return (csin(a) * ccosh(b) + I*(ccos(a)*csinh(b)));
}

complex double cos_complex(int v) {
	double complex a = fmod(creall(v), (2*PI));
	double complex b = complex_modulo(cimagl(v),2*PI*I);

	return (ccos(a) * ccosh(b) - I*(csin(a)*csinh(b)));
}

complex double* transfer_matrix(double k_0, double n, double d, double theta) {
	complex double k_z = k_0*n*ccos(theta);

	complex double q_1 = cos_complex(k_z * d);
	complex double q_2 = I*sin_complex(k_z * d);
	complex double n_cos_th = n*ccos(theta);
	complex double array[] = {{q_1, q_2/n_cos_th},{n_cos_th*q_2,q_1}};

	return array;
}

int main() {
}