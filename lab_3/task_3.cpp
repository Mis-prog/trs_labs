#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;
#define PI 3.14159

const double epsilon = 1.0e-10;

double f(double x, double y) {
    return -PI * PI * sin(PI * y) * (5. * cos(PI * x) * cos(PI * x) - 3.);
}

inline double Acc(double *arr, int N, int i, int j) {
    double temp = 0;
    int K = sqrt(N);
    if ((abs(i - j) == K) && (j > i)) {
        temp = arr[i * 5 + 4];
        return temp;
    }
    if ((abs(i - j) == K) && (j < i)) {
        temp = arr[i * 5 + 0];
        return temp;
    }
    if (abs(i - j) <= 1) {
        temp = arr[i * 5 + j - i + 2];
        return temp;
    }
    return temp;
}

double VecScalarProduct(double *a, double *b, int N) {
    double result = 0;
    for (int i = 0; i < N; i++)
        result += a[i] * b[i];
    return result;
}

double *VecSum(double *a, double *b, double a_coef, double b_coef, int N) {
    double *c = new double[N];
    for (int i = 0; i < N; i++)
        c[i] = a_coef * a[i] + b_coef * b[i];
    return c;
}

double VecNorm(double *a, int N) {
    double result = abs(a[0]);
    for (int i = 1; i < N; i++) {
        if (abs(a[i]) > result)
            result = abs(a[i]);
    }
    return result;
}

double *MatrixVectorProduct(double *A, double *v, int N) {
    double *b = new double[N];
    double temp = 0;
    double elemA = 0;
    for (int i = 0; i < N; i++) {
        temp = 0;
        for (int k = 0; k < N; k++) {
            elemA = Acc(A, N, i, k);
            temp += elemA * v[k];
        }
        b[i] = temp;
    }
    return b;
}

double *CGM(double *A, double *b, double *x_0, int N, int maxiter, double acc) {
    double *x_k = new double[N];
    double *x_k1 = new double[N];
    double *p_k = new double[N];
    double *p_k1 = new double[N];
    double *r_k = new double[N];
    double *r_k1 = new double[N];
    double *q_k;
    double alpha_k = 0, beta_k = 0;
    double delta = 0;
    int iter = 0;

    for (int i = 0; i < N; i++)
        x_k[i] = x_0[i];

    r_k = VecSum(b, MatrixVectorProduct(A, x_k, N), 1.0, -1.0, N);
    for (int i = 0; i < N; i++)
        p_k[i] = r_k[i];

    do {
        iter++;
        q_k = MatrixVectorProduct(A, p_k, N);
        alpha_k = -VecScalarProduct(r_k, r_k, N) / VecScalarProduct(p_k, q_k, N);
        x_k1 = VecSum(x_k, p_k, 1.0, -alpha_k, N);
        r_k1 = VecSum(r_k, q_k, 1.0, alpha_k, N);
        delta = VecScalarProduct(r_k1, r_k1, N);
        //cout << delta << endl;
        if (delta < acc)
            break;
        beta_k = VecScalarProduct(r_k1, r_k1, N) / VecScalarProduct(r_k, r_k, N);
        p_k1 = VecSum(r_k1, p_k, 1.0, beta_k, N);
        for (int i = 0; i < N; i++) {
            x_k[i] = x_k1[i];
            r_k[i] = r_k1[i];
            p_k[i] = p_k1[i];
        }
    } while (iter < maxiter);
    cout << iter << " Iterations" << endl;
    return x_k1;
}

double Solution(double x, double y) {
    //return (exp(x*y)-1)*(1-x)*(1-y);
    return sin(PI * x) * sin(PI * x) * sin(PI * y);
}

double ChebNormErr(double *y, double lx, double ly, int N) {
    double hx = lx / N;
    double hy = (ly * hx) / lx;
    double max = -1.0;
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            if (abs(y[j * (N + 1) + i] - Solution(hx * i, hy * j)) >= max)
                max = abs(y[j * (N + 1) + i] - Solution(hx * i, hy * j));
        }
    }
    return max;
}

int main(int argc, char *argv[]) {
    int N = 50;

    if (argc == 2) {
        N = atoi(argv[1]);
    }
    double lx = 1, ly = 1;
    double hx = lx / N;
    double hy = (ly * hx) / lx;
    double *A = new double[(N - 1) * (N - 1) * 5];
    double *y_0 = new double[(N - 1) * (N - 1)];
    double *y = new double[(N + 1) * (N + 1)];
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < (N - 1) * (N - 1); j++) {

            if (i == 0) {
                if (j <= N - 2) A[j * 5 + i] = 0;
                else A[j * 5 + i] = 1 / pow(hx, 2);
            }
            if (i == 1) {
                if ((j == 0) || (j % (N - 1) == 0)) A[j * 5 + i] = 0;
                else A[j * 5 + i] = 1 / pow(hy, 2);
            }
            if (i == 2) A[j * 5 + i] = -2 * (1 / pow(hx, 2) + 1 / pow(hy, 2));
            if (i == 3) {
                if ((j == (N - 1) * (N - 1)) || ((j + 1) % (N - 1) == 0)) A[j * 5 + i] = 0;
                else A[j * 5 + i] = 1 / pow(hy, 2);
            }
            if (i == 4) {
                if (j >= (N - 2) * (N - 1)) A[j * 5 + i] = 0;
                else A[j * 5 + i] = 1 / pow(hx, 2);
            }
        }
    }
    double *b = new double[(N - 1) * (N - 1)];
    for (int i = 0; i < (N - 1); i++) {
        for (int j = 0; j < (N - 1); j++) {
            b[j * (N - 1) + i] = -f((i + 1) * hx, (j + 1) * hy);
            y_0[j * (N - 1) + i] = 0;
        }
    }

    y_0 = CGM(A, b, y_0, (N - 1) * (N - 1), 1000, epsilon);
    for (int i = 0; i < (N + 1); i++) {
        for (int j = 0; j < (N + 1); j++) {
            if ((i == 0) || (i == N) || (j == 0) || (j == N))
                y[j * (N + 1) + i] = 0;
            else
                y[j * (N + 1) + i] = y_0[(j - 1) * (N - 1) + i - 1];
        }
    }
    ofstream filewrite("CGM.txt");
    for (int i = 0; i < N + 1; i++) {
        for (int j = 0; j < N + 1; j++)
            filewrite << i * hx << " " << j * hy << " " << y[j * (N + 1) + i] << endl;
    }

    cout << ChebNormErr(y, lx, ly, N) << endl;
    return 0;
}