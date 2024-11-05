#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <omp.h>

complex double complex_modulo(complex double z1, complex double z2) {
    double abs_z1 = cabs(z1);
    double abs_z2 = cabs(z2);
    double remainder = fmod(abs_z1, abs_z2);
    double angle_z1 = carg(z1);
    complex double result = remainder * cexp(I * angle_z1);
    return result;
}

complex double sin_complex(complex double v) {
    double a = fmod(creal(v), 2 * M_PI);
    complex double b = complex_modulo(cimag(v) * I, 2 * M_PI * I);
    return (sin(a) * cosh(cimag(b))) + I * (cos(a) * sinh(cimag(b)));
}

complex double cos_complex(complex double v) {
    double a = fmod(creal(v), 2 * M_PI);
    complex double b = complex_modulo(cimag(v) * I, 2 * M_PI * I);
    return cos(a) * cosh(cimag(b)) - I * sin(a) * sinh(cimag(b));
}

void transfer_matrix(complex double result[2][2], double k_0, complex double n, double d, double theta) {
    complex double k_z = k_0 * n * cos(theta);
    complex double q_1 = cos_complex(k_z * d);
    complex double q_2 = I * sin_complex(k_z * d);
    complex double n_cos_th = n * cos(theta);
    result[0][0] = q_1;
    result[0][1] = q_2 / n_cos_th;
    result[1][0] = n_cos_th * q_2;
    result[1][1] = q_1;
}

void solve_tmm(double* R, double* T, complex double layers[][2], int num_layers, double wavelength, double theta) {
    double k_0 = (2 * M_PI) / wavelength;
    complex double M[2][2] = { {1, 0}, {0, 1} };
    complex double n_0 = layers[0][0];
    complex double n_l = layers[num_layers - 1][0];

    #pragma omp parallel for
    for (int i = 1; i < num_layers - 1; ++i) {
        complex double n = layers[i][0];
        double d = creal(layers[i][1]);
        double theta_i = asin(creal(n_0) * sin(theta) / creal(n));
        complex double T_i[2][2];
        transfer_matrix(T_i, k_0, n, d, theta_i);

        #pragma omp critical
        {
            complex double temp[2][2];
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    temp[j][k] = 0;
                    for (int l = 0; l < 2; ++l) {
                        temp[j][k] += M[j][l] * T_i[l][k];
                    }
                }
            }
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    M[j][k] = temp[j][k];
                }
            }
        }
    }

    complex double q_1 = n_0 * cos(theta);
    double theta_l = asin(creal(n_0) * sin(theta) / creal(n_l));
    complex double q_2 = n_l * cos(theta_l);
    complex double M_in[2][2] = {
        {1, 1},
        {q_1, -q_1}
    };
    complex double M_out[2][2] = {
        {1, 1},
        {q_2, -q_2}
    };

    complex double M_total[2][2];
    complex double M_temp[2][2];
    complex double M_inv[2][2];

    complex double det = M_in[0][0] * M_in[1][1] - M_in[0][1] * M_in[1][0];
    M_inv[0][0] = M_in[1][1] / det;
    M_inv[0][1] = -M_in[0][1] / det;
    M_inv[1][0] = -M_in[1][0] / det;
    M_inv[1][1] = M_in[0][0] / det;

    for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 2; ++k) {
            M_temp[j][k] = 0;
            for (int l = 0; l < 2; ++l) {
                M_temp[j][k] += M[j][l] * M_out[l][k];
            }
        }
    }

    for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 2; ++k) {
            M_total[j][k] = 0;
            for (int l = 0; l < 2; ++l) {
                M_total[j][k] += M_inv[j][l] * M_temp[l][k];
            }
        }
    }

    complex double r = M_total[1][0] / M_total[0][0];
    complex double t = 1 / M_total[0][0];
    *R = pow(cabs(r), 2);
    *T = pow(cabs(t), 2) * creal(n_l * cos(theta_l)) / creal(n_0 * cos(theta));
}

int main() {
    clock_t begin = clock();
    double R, T;
    complex double layers[3][2] = {
        {1.0 + 0.0 * I, 0.0},
        {1.5 + 0.0 * I, 100.0},
        {1.0 + 0.0 * I, 0.0}
    };

    #pragma omp parallel
    {
        #pragma omp single
        solve_tmm(&R, &T, layers, 3, 500.0, 0.0);
    }

    printf("Reflectance: %f\n", R);
    printf("Transmittance: %f\n", T);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("%f", time_spent);
}