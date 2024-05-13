#define _USE_MATH_DEFINES

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <iomanip>

using namespace std;

double F(double x, double t) {
    return 0;
}

double EFDS(int N, double h, double tau) {
    cout << "> EFDS started..." << endl;
    double *X = new double[N + 1];
    double *T = new double[N + 1];
    double **u = new double *[N + 1];

    for (int i = 0; i < N + 1; i++) {
        X[i] = i * h;
        T[i] = i * tau;
    }

    for (int t = 0; t < N + 1; t++) {
        u[t] = new double[N + 1];
    }

    for (int i = 0; i < N + 1; i++) {
        for (int j = 0; j < N + 1; j++) {
            u[i][j] = 0.0;
        }
    }

    for (int j = 0; j < N + 1; j++) {
        u[0][j] = sinh(X[j]);
    }

    double k, FF;
    k = 1;

    for (int t = 0; t < N; t++) {
        for (int x = 1; x < N; x++) {
            FF = F(X[x], T[t]);
            u[t + 1][x] =
                    (tau * k / (h * h)) * (u[t][x - 1] + u[t][x + 1]) - tau * (2. / (h * h) - 1. / tau) * u[t][x] +
                    tau * FF;
        }
        u[t + 1][0] = u[t + 1][1] - h;
        u[t + 1][N] = 1 - T[t + 1];

    }

    double pogr = 0;
    ofstream fout1("EFDS-" + to_string(N) + ".txt");
    ofstream fout2("an-" + to_string(N) + ".txt");

    for (int m = 0; m < N + 1; m++) {
        for (int j = 0; j < N + 1; j++) {
            double analitik;
            analitik = exp(T[m]) * sinh(X[j]);
            double delta = abs(analitik - u[m][j]);

            if (delta > pogr) {
                pogr = delta;
            }
            if (m == N - 1)
                cout << "  x = " << setprecision(10) << X[j] << " numsol = " << u[m][j] << " analitik = " << analitik
                     << endl;
        }
    }

    fout1.close();
    fout2.close();

    cout << "max_pogr = " << setprecision(10) << pogr << endl;
    return 0;
}


int main(int argc, char *argv[]) {
    int N = 100;
    double h, tau;
    double a = 0., b = 1.;

    if (argc == 2) {
        N = atoi(argv[1]);
    }

    h = (double) ((b - a) / double(N));
    tau = h * h / 2.;

    cout << "N = " << N << "; h = " << h << "; tau = " << tau << endl;

    EFDS(N, h, tau);

    return 0;
}