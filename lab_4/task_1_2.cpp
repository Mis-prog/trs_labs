
#define _USE_MATH_DEFINES

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <iomanip>

using namespace std;

const double C = 1.;

double phi(double x, double t) {
    return (1/6) * sin(5*x) + sin(x);
}

double F(double x, double t) {
    return (1/6) * (cos(t+5*x) - 6*cos(t-x)) + cos(t-x) + (5/6)*cos(t+5*x);
}

double psi(double x, double t) {
    return (1/6) * sin(t) + sin(-t);
}

double ExplicitLeftCorner(int N, int M, double h, double tau) {
    double *X = new double[N + 1];
    double *T = new double[M + 1];
    double **u = new double*[M + 1];

    double gamma = tau / h;

    for (int i = 0; i < N + 1; i++)
    {
        X[i] = i * h;
    }

    for (int i = 0; i < M + 1; i++)
    {
        T[i] = i * tau;
    }

    for (int t = 0; t < M + 1; t++) {
        u[t] = new double[N + 1];
    }

    for (int i = 0; i < M + 1; i++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            u[i][j] = 0.0;
        }
    }

    for (int j = 0; j < N + 1; j++)
    {
        u[0][j] = phi(j*h, 0);
    }

    for (int t = 0; t < M; t++)
    {
        for (int x = 1; x < N+1; x++)
        {
            u[t + 1][x] = (1. - gamma)*u[t][x] + gamma * u[t][x - 1] + F(X[x], T[t])*tau;
        }
        u[t + 1][0] = psi(0, (t+1)*tau);
    }

    double pogr = 0.;

    ofstream fout("ELC-"+to_string(N)+"-"+to_string(M)+".txt");
    for (int t = 0; t < M + 1; t++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            double exact;
            exact = (1/6)*sin(T[t]+5*X[j]) + sin(X[j]-T[t]);

            double delta = abs(exact - u[t][j]);

            if (delta >= pogr) {
                pogr = delta;
            }
            if (t == M)
                fout << X[j] << " " << delta << endl;
        }
    }
    fout.close();
    cout << "Error = " << pogr << endl;
    return 0;
}

double ImplicitLeftCorner(int N, int M, double h, double tau) {
    double *X = new double[N + 1];
    double *T = new double[M + 1];
    double **u = new double*[M + 1];

    double gamma = tau / h;

    for (int i = 0; i < N + 1; i++)
    {
        X[i] = i * h;
    }

    for (int i = 0; i < M + 1; i++)
    {
        T[i] = i * tau;
    }

    for (int t = 0; t < M + 1; t++) {
        u[t] = new double[N + 1];
    }

    for (int i = 0; i < M + 1; i++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            u[i][j] = 0.0;
        }
    }

    for (int j = 0; j < N + 1; j++)
    {
        u[0][j] = phi(j*h, 0);
    }

    double Func;

    for (int t = 0; t < M; t++)
    {
        u[t + 1][0] = psi(0, (t + 1)*tau);
        for (int x = 1; x < N+1; x++)
        {
            u[t + 1][x] = (u[t][x] + (tau / h) * u[t + 1][x - 1] + F(X[x], T[t+1])*tau) / (1. + tau / h);
        }
        u[t + 1][0] = psi(N, (t+1)*tau);
    }

    double pogr = 0;

    ofstream fout("ILC-"+to_string(N)+"-"+to_string(M)+".txt");
    for (int t = 0; t < M + 1; t++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            double exact;
            exact = (1/6)*sin(T[t]+5*X[j]) + sin(X[j]-T[t]);

            double delta = abs(exact - u[t][j]);

            if (delta >= pogr) {
                pogr = delta;
            }
            if (t == M)
                fout << X[j] << " " << T[t] << " " << u[t][j] << endl;
        }
    }
    fout.close();
    cout << "Error = " << pogr << endl;
    return 0;
}

int main(int argc, char* argv[])
{
    int N, M;
    double h, tau;
    double a = 0., b = 2., T = 10.;

    if (argc == 3)
    {
        N = atoi(argv[1]);
        M = atoi(argv[2]);
    }
    else
    {
        N = 100;
        M = 500;
    }

    h = (double)((b - a) / double(N));
    tau = (double)T / (double)M;

    if (tau/h > 1) {
        cout << "[!] May cause unstable solution. 'M' should be greater than " << (int)(T/h) << endl;
        //return 0;
    }

    cout << "N = " << N << "; M = " << M << endl;
    cout << "h = " << h << "; tau = " << tau << endl;

    //ExplicitLeftCorner(N, M, h, tau);
    ImplicitLeftCorner(N, M, h, tau);

    return 0;
}