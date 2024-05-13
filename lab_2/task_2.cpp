#define _USE_MATH_DEFINES

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <iomanip>

using namespace std;

double F(double x, double t) {
    return 0.;
}

double *progonka(int N, double *A, double *B, double *C, double *f, double *y) {
    double *v, *u;
    v = new double[N + 1];
    u = new double[N + 1];
    v[0] = 0;
    u[0] = 0;

    for (int i = 1; i < N + 1; i++) {
        v[i] = -A[i - 1] / (B[i - 1] + C[i - 1] * v[i - 1]);
        u[i] = (f[i - 1] - C[i - 1] * u[i - 1]) / (B[i - 1] + C[i - 1] * u[i - 1]);
    }

    y[N] = u[N];
    for (int i = N - 1; i >= 0; i--) {
        y[i] = v[i + 1] * y[i + 1] + u[i + 1];
    }
    return y;
}

void IFDS(int N, double h, double tau) {
    double *A = new double[N + 1];
    double *B = new double[N + 1];
    double *C = new double[N + 1];
    double *f = new double[N + 1];
    double *x = new double[N + 1];
    double *uPrev = new double[N + 1];
    double *uNext = new double[N + 1];
    double *t = new double[N + 1];

    for (int i = 0; i < N + 1; i++) {
        x[i] = i * h;
        uPrev[i] = sinh(x[i]);
        uNext[i] = 0;
    }

    for (int i = 0; i < N + 1; i++) {
        t[i] = i * tau;
    }
    double pogr = 0;

    for (int tt = 0; tt < N + 1; tt++) {
        A[0] = 0.;
        B[0] = -1.;
        C[0] = 1.;

        A[N] = 0.;
        B[N] = 1.;
        C[N] = 0.;

        f[0] = 1. * h;
        f[N] = 1 - t[tt + 1];

        for (int j = 1; j < N; j++) {
            A[j] = tau / (2 * h * h);
            B[j] = 1 + (tau / h * h);
            C[j] = tau / (2 * h * h);

            f[j] = A[j] * uPrev[j + 1] + B[j] * uPrev[j] + C[j] * uPrev[j - 1];
        }

        progonka(N + 1, A, B, C, f, uNext);

        for (int i = 0; i < N + 1; i++) {
            uPrev[i] = uNext[i];
            double analitik = exp(t[tt]) * sinh(x[i]);
            double delta = abs(analitik - uNext[i]);
            if (delta > pogr) {
                pogr = delta;
            }
        }
    }

    cout << "IFDS Error = " << setprecision(10) << pogr << endl;
}

void CrankNicolson(int N, double h, double tau) {
    double *A = new double[N + 1];
    double *B = new double[N + 1];
    double *C = new double[N + 1];

    double *f = new double[N + 1];
    double *x = new double[N + 1];

    double *uPrev = new double[N + 1];
    double *uNext = new double[N + 1];
    double *t = new double[N + 1];

    for (int i = 0; i < N + 1; i++) {
        x[i] = i * h;
        t[i] = i * tau;
        uPrev[i] = sinh(x[i]);
        uNext[i] = 0;
    }

    double pogr = 0;


    for (int tt = 0; tt < N + 1; tt++) {

        A[0] = 0.;
        B[0] = -1.;
        C[0] = 1.;

        A[N] = 0.;
        B[N] = 1.;
        C[N] = 0.;

        f[0] = 1. * h;
        f[N] = 1 - t[tt + 1];

        for (int j = 1; j < N; j++) {
            A[j] = -(tau / 2 * h * h);
            B[j] = 1 + (tau / h * h);
            C[j] = -(tau / 2 * h * h);

            f[j] = A[j] * uPrev[j + 1] + B[j] * uPrev[j] + C[j] * uPrev[j - 1] + F(x[j], t[tt]) * h * h;
        }

        progonka(N + 1, A, B, C, f, uNext);

        for (int i = 0; i < N + 1; i++) {
            uPrev[i] = uNext[i];
            double analitik;
            analitik = exp(t[tt]) * sinh(x[i]);

            double delta = abs(analitik - uNext[i]);
            if (delta > pogr) {
                pogr = delta;
            }
        }
    }
    cout << "CrankNicolson Error = " << setprecision(10) << pogr << endl;
}

int main(int argc, char *argv[]) {
    int N = 100;
    double h, tau;
    double a = 0., b = 1.;

    if (argc == 2) {
        N = atoi(argv[1]);
    }

    h = (double) ((b - a) / double(N));
    tau = (h * h) / 2;

    cout << "N = " << N << ", h = " << h << ", tau = " << tau << endl;

    IFDS(N, h, tau);

    CrankNicolson(N, h, tau);

    return 0;
}