#define _USE_MATH_DEFINES

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <iomanip>

using namespace std;

#define eps 0.01

double fi(double x, double t) {
    return sinh(x);
}

double F(double x, double t) {
    return 0.;
}

double psiOne(double x, double t) {
    return exp(t);
}

double psiTwo(double x, double t) {
    return exp(t + 1);
}

double K(double x, double t) {
    return 1.;
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

void CS(int N, double h, double tau) {
    cout << "\n > Conservative Scheme" << endl;
    double *A = new double[N + 1];
    double *B = new double[N + 1];
    double *C = new double[N + 1];
    double *f = new double[N + 1];
    double *x = new double[N + 1];
    double *uPrev = new double[N + 1];
    double *uNext = new double[N + 1];
    double *t = new double[N + 1];
    double coef = 1.;
    double ke, kw;

    for (int i = 0; i < N + 1; i++) {
        x[i] = i * h;
        t[i] = i * tau;
        uPrev[i] = sinh(x[i]);
        uNext[i] = 0;
    }

    double pogr = 0;

    float timeStart = clock() / (float) CLOCKS_PER_SEC;

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
            kw = (uPrev[j - 1] + uPrev[j]) / 2.;
            ke = (uPrev[j + 1] + uPrev[j]) / 2.;

            A[j] = (-coef / (2 * h * h)) * kw;
            B[j] = (1. / tau + coef * ke / (h * h) + coef * kw / (h * h));
            C[j] = (-coef / (2 * h * h)) * ke;

            f[j] = uPrev[j] * (1. / tau - (1. - coef) * (ke) / (h * h) - (1. - coef) * (kw) / (h * h)) +
                   uPrev[j + 1] * (1. - coef) * (ke) / (h * h) - uPrev[j - 1] * (1. - coef) * (kw) / (h * h) +
                   uPrev[j] * uPrev[j] * F(x[j], t[tt + 1]);
            //cout << A[j] << ", " << B[j] << ", "  << C[j] << ", "  << f[j] << endl;
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

    float timeStop = clock() / (float) CLOCKS_PER_SEC;
    cout << "Time = " << timeStop - timeStart << endl;
    cout << "Error = " << setprecision(10) << pogr << endl;
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
    tau = 5e-09;

    cout << "N = " << N << ", h = " << h << ", tau = " << tau << endl;

    CS(N, h, tau);

    return 0;
}