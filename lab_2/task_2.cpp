#include <iostream>
#include <vector>
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

double u1_func(double x, double t) {
    return exp(t) * sinh(x);
//    return x * cos(t) + t * (1 - x) * (1 - x);
}


void z2_p1() {
    double ht, hx;
    vector<int> nX{ 10,100,130 };
    vector<int> nT = nX;
    for (int i = 0; i < nX.size(); i++) {

        clock_t start = clock();
        hx = 1. / (nX[i] - 1);
//        nT[i] = 2 * (nX[i]) * (nX[i]);
        ht = 1. / (nT[i] - 1);
        vector<vector<double>> M(nT[i], vector<double>(nX[i])); // t , x ;
        ofstream file_to_cout("../../../trs_labs_/output/lab_2/z2_p1" + to_string(i + 1) + ".txt");
        vector<double> A(nX[i]), B(nX[i]), C(nX[i]), D(nX[i]);

        for (int k = 0; k < nX[i]; k++) {
            M[0][k] = sinh(k*hx);
        }
        A[0] = 0;
        B[0] = hx-1;
        C[0] = 1.;
        for (int k = 1; k < nX[i] - 1; k++) {
            A[k] = -1. / hx / hx;
            B[k] = 1. / ht + 2. / hx / hx;
            C[k] = -1. / hx / hx;
        }
        A[nX[i] - 1] = -1.;
        B[nX[i] - 1] = hx+1;
        C[nX[i] - 1] = 0;
        for (int j = 1; j < nT[i]; j++) {
            D[0] = hx * exp(j * ht);
            for (int k = 1; k < nX[i] - 1; k++) {
                D[k] = M[j - 1][k] / ht;
            }
            D[nX[i] - 1] =hx * exp(j * ht + 1);

            M[j] = TripleDiag(A, B, C, D);
        }

        clock_t end = clock();
        double time = (double)(end - start) / CLOCKS_PER_SEC;

        for (int j = 0; j < nT[i]; j++) {//cout
            for (int k = 0; k < nX[i]; k++) {
                file_to_cout << M[j][k] << " ";
            }
            file_to_cout << endl;
        }

        double dif;
        double max_dif;
        for (int j = 0; j < nT[i]; j++)
        {
            for (int k = 1; k < nX[i]; k++)
            {
                max_dif = 0;
                dif = abs(u1_func(k * hx, j * hx) - M[j][k]);
                if (dif > max_dif)
                {
                    max_dif = dif;
                }
            }
        }

        cout << "Step x: " << hx << setw(13 + nX.size()) << "Step T: " << ht
             << setw(8 + nX.size()) << "Error: " << max_dif << setw(8 + nX.size()) << "Time: " << time << endl;
    }
}

void z2_p2() {
    double ht, hx;
    vector<double> nX{ 10,100,400 };
    vector<double> nT = nX;
    for (int i = 0; i < nX.size(); i++) {
        clock_t start = clock();
        hx = 1. / (nX[i] - 1);
        //nT[i] = 2 * (nX[i]) * (nX[i]);
        ht = 1. / (nT[i] - 1);
        vector<vector<double>> M(nT[i], vector<double>(nX[i])); // t , x ;
        ofstream file_to_cout("../../../trs_labs_/output/lab_2/z2_p2" + to_string(i + 1) + ".txt");
        vector<double> A(nX[i]), B(nX[i]), C(nX[i]), D(nX[i]);

        for (int k = 0; k < nX[i]; k++) {
            M[0][k] = sinh(k*hx);
        }
        for (int j = 1; j < nT[i]; j++) {
            A[0] = 0;
            B[0] = hx-1;
            C[0] = 1.;
            D[0] = hx * exp(j * ht);
            for (int k = 1; k < nX[i] - 1; k++) {
                A[k] = -1. / 2. / hx / hx;
                B[k] = 1. / ht + 1. / hx / hx;
                C[k] = -1. / 2. / hx / hx;
                D[k] = 1. / 2. / hx / hx * (M[j - 1][k - 1] + M[j - 1][k + 1] - 2 * M[j - 1][k]) + 1. / ht * M[j - 1][k];
            }
            A[nX[i] - 1] = -1.;
            B[nX[i] - 1] = hx+1;
            C[nX[i] - 1] = 0;
            D[nX[i] - 1] = hx * exp(j * ht + 1);
            M[j] = TripleDiag(A, B, C, D);
        }

        clock_t end = clock();
        double time = (double)(end - start) / CLOCKS_PER_SEC;


        for (int j = 0; j < nT[i]; j++) {//cout
            for (int k = 0; k < nX[i]; k++) {
                file_to_cout << M[j][k] << " ";
            }
            file_to_cout << endl;
        }

        double dif;
        double max_dif;
        for (int j = 0; j < nT[i]; j++)
        {
            for (int k = 1; k < nX[i]; k++)
            {
                max_dif = 0;
                dif = abs(u1_func(k * hx, j * hx) - M[j][k]);
                if (dif > max_dif)
                {
                    max_dif = dif;
                }
            }
        }

        cout << "Step x: " << hx << setw(13 + nX.size()) << "Step T: " << ht
             << setw(8 + nX.size()) << "Error: " << max_dif << setw(8 + nX.size()) << "Time: " << time << endl;
    }
}

int main(){
    z2_p1();
    z2_p2();

    return 0;
}