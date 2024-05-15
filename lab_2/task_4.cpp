#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <fstream>
#include <math.h>
#include <iomanip>

using namespace std;

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

double F(double u) {
    return u * u * u;
    //return 1;
}

double func(double x, double t) {
    return ((1. - x) * (1. - x) - x * sin(t) - 2. * t);
}

double k1(double u) {
    return sin(u);
    //return 1;
}

void z4(double eps) {
    vector<double> nX{10, 100}; // по x
    vector<double> nT = nX;
    double delta = 0;
    for (int k = 0; k < nX.size(); k++) {
        ofstream result_file("../../../trs_labs_/output/lab_2/z4_u" + to_string(k) + ".txt");

        auto begin = std::chrono::steady_clock::now();

        double ht = 1 / (nT[k] - 1);
        double hx = 1 / (nX[k] - 1);

        vector<vector<double>> u(nT[k], vector<double>(nX[k]));
        vector<double> previous;
        vector<double> A(nX[k]), B(nX[k]), C(nX[k]), D(nX[k]);

        for (int i = 0; i < nX[k]; i++) {
            u[0][i] = i * hx; // началка
        }

        for (int j = 1; j < nT[k]; j++) {
            //граничные
            A[0] = 0;
            B[0] = 1;
            C[0] = 0;
            D[0] = j * ht;
            A[nX[k] - 1] = -1. / hx;
            B[nX[k] - 1] = 1. + 1. / hx;
            C[nX[k] - 1] = 0;
            D[nX[k] - 1] = 2 * cos(j * ht);
            //
            double g = hx * hx / ht;
            for (int i = 1; i < nX[k] - 1; i++) {
                double ke = k1((u[j - 1][i + 1] + u[j - 1][i]) / 2);
                double kw = k1((u[j - 1][i] + u[j - 1][i - 1]) / 2);
                A[i] = -kw / hx / hx;
                B[i] = kw / hx / hx + ke / hx / hx + 1. / ht;
                C[i] = -ke / hx / hx;
                D[i] = func(i * hx, j * ht) * F(u[j - 1][i]) / hx + 1. / ht * u[j - 1][i];
            }
            u[j] = TripleDiag(A, B, C, D);


            int iter = 0;
            do {
                delta = 0;
                iter++;
                previous.assign(u[j].begin(), u[j].end());
                for (int i = 1; i < nX[k] - 1; i++) {
                    double ke = k1((u[j][i + 1] + u[j][i]) / 2);
                    double kw = k1((u[j][i] + u[j][i - 1]) / 2);
                    A[i] = -kw / hx / hx;
                    B[i] = kw / hx / hx + ke / hx / hx + 1. / ht;
                    C[i] = -ke / hx / hx;
                }
                u[j] = TripleDiag(A, B, C, D);

                for (int y = 1; y < nX[k]; y++) {
                    double diff = abs(k1(u[j][y] / 2 + u[j][y - 1] / 2) - k1(previous[y] / 2 + previous[y - 1] / 2));
                    if (diff > delta) {
                        delta = diff;
                    }
                }
                previous.clear();
            } while (delta > eps);
            cout << "step x,t: " << hx << setw(13 + nX.size()) << "iter: " << iter << " error: " << delta << endl;
        }

        auto end = std::chrono::steady_clock::now();

        auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);

        for (int i = 0; i < nX[k]; i++) {
            for (int j = 0; j < nT[k]; j++) {
                result_file << u[j][i] << " ";
            }
            result_file << endl;
        }
        cout << "Step x: " << hx << setw(13 + nX.size()) << "Step T: " << ht
             << setw(8 + nX.size()) << "Error: " << "..." << setw(8 + nX.size()) << "Time, ms: "
             << elapsed_ms.count() * 1e-3 << endl;

    }
}


int main() {
    z4(1e-6);

    return 0;
}