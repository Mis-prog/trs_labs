#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <math.h>


using namespace std;


double F(double u) {
    return u * u * u;
    //return 1;
}

double func(double x, double t) {
    return ((1. - x) * (1. - x) - x * sin(t) - 2. * t);
}

double k1(double u) {
    return u * u + 1;
    //return 1;
}

vector<double> TripleDiag(vector<double> A, vector<double> B, vector<double> C, vector<double> D) {
    vector<double> P(A.size());
    vector<double> Q(A.size());
    vector<double> X(A.size());
    int N = A.size();
    P[0] = C[0] / B[0];
    Q[0] = D[0] / B[0];
    for (int i = 1; i < N; ++i) {
        if (i < N - 1)
            P[i] = C[i] / (B[i] - A[i] * P[i - 1]);
        Q[i] = (D[i] - A[i] * Q[i - 1]) / (B[i] - A[i] * P[i - 1]);
    }

    // backward
    X[N - 1] = Q[N - 1];
    for (int i = N - 2; i >= 0; --i) {
        X[i] = Q[i] - P[i] * X[i + 1];
    }

    return X;
}

double u1_func(double x, double t) {
    return exp(t) * sinh(x);
}
void z3() {
    vector<double> nX{10, 450}; // по x
    vector<double> nT{10, 450};// по t
    for (int k = 0; k < nX.size(); k++) {
        ofstream result_file("../../../trs_labs_/output/lab_2/z3_u" + to_string(k) + ".txt");
        ofstream result_file_error("../../../trs_labs_/output/lab_2/z3_u_error" + to_string(k) + ".txt");

        auto begin = std::chrono::steady_clock::now();
        double ht = 1. / (nT[k] - 1);
        double hx = 1. / (nX[k] - 1);

        vector<vector<double>> u(nT[k], vector<double>(nX[k]));
        for (int i = 0; i < nX[k]; i++) {
            u[0][i] = sinh(i * hx); // началка
        }

        for (int j = 1; j < nT[k]; j++) {
            vector<double> A(nX[k]), B(nX[k]), C(nX[k]), D(nX[k]);
            //граничные
            A[0] = 0;
            B[0] = hx - 1;
            C[0] = 1.;
            D[0] = hx * exp(j * ht);
            A[nX[k] - 1] = -1.;
            B[nX[k] - 1] = hx + 1;
            C[nX[k] - 1] = 0;
            D[nX[k] - 1] = hx * exp(j * ht + 1);
            //
            double g = hx * hx / ht;
            for (int i = 1; i < nX[k] - 1; i++) {
                double ke = k1((u[j - 1][i + 1] + u[j - 1][i]) / 2);
                double kw = k1((u[j - 1][i] + u[j - 1][i - 1]) / 2);
                A[i] = -kw / hx / hx;
                B[i] = 1. / ht + (kw + ke) / hx / hx;
                C[i] = -ke / hx / hx;
                D[i] = u[j - 1][i] / ht;

            }
            u[j] = TripleDiag(A, B, C, D);
        }

        auto end = std::chrono::steady_clock::now();
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);

        for (int i = 0; i < nX[k]; i++) {
            for (int j = 0; j < nT[k]; j++) {
                result_file << u[j][i] << " ";
            }
            result_file << endl;
        }
        double dif;
        double max_dif;
        for (int i = 0; i < nT[k]; i++) {
            for (int j = 0; j < nX[k]; j++) {
                max_dif = 0;
                dif = abs(u1_func(j * hx, i * ht) - u[i][j]);
                result_file_error << dif << " ";
                if (dif > max_dif) {
                    max_dif = dif;
                }
            }
            result_file_error << endl;
        }

        cout << "Step x: " << hx << setw(13 + nX.size()) << "Step T: " << ht
             << setw(8 + nX.size()) << "Error: " << max_dif << setw(8 + nX.size()) << "Time, ms: "
             << elapsed_ms.count() * 1e-3 << endl;
    }
}

int main() {
    z3();

    return 0;
}