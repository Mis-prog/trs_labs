#define _USE_MATH_DEFINES

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <iomanip>

#define eps 0.01
using namespace std;

double phi(double x, double t) {
    return 0.;
}

double F(double x, double t) {
    return sin(2*t);
}

double uExact(double x, double t)
{
    return -(sin(2*x + 2*t) + sin(2*x - 2*t) - 2*sin(2*x)) / 8;
}

double psi0(double x, double t)
{
    return -sin(2*t);
}

double psi1(double x, double t)
{
    return cos(2) * sin(t) * sin(t);
}

double ExplicitCross(int N, int M, double h, double tau) {
    double *X = new double[N + 1];
    double *T = new double[M + 1];
    double **u = new double*[M + 1];

    for (int i = 0; i < N + 1; i++)
    {
        X[i] = i * h;
    }

    for (int i = 0; i < M + 1; i++)
    {
        T[i] = i * tau;
        u[i] = new double[N + 1];
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
        u[0][j] = phi(X[j], 0);
        u[M][j] = phi(X[j], 0);
    }

    for (int t = 1; t < M; t++)
    {
        for (int x = 1; x < N+1; x++)
        {
            u[t+1][x] = 2*u[t][x] - u[t-1][x] + tau*tau * (u[t][x-1]/h/h - 2 * u[t][x]/h/h + u[t][x+1]/h/h + F(X[x], T[t]));
        }
        u[t + 1][0] = -h * psi0(0, (t + 1)*tau) + u[t+1][1];
        u[t + 1][N] = psi1(N, (t + 1)*tau);
    }

    double pogr = 0;

    ofstream fout("CrossILC-"+to_string(N)+"-"+to_string(M)+".txt");
    for (int t = 0; t < M + 1; t++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            fout << X[j] << " " << T[t] << " " << u[t][j] << endl;
        }
    }
    fout.close();
    return 0;
}

int main(int argc, char* argv[])
{
    int N, M;
    double h, tau;
    double a = 0., b = 1., T = 5.;

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

    ExplicitCross(N, M, h, tau);

    return 0;
}

