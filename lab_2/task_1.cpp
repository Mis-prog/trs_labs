#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>
#include <chrono>
#include <math.h>

using namespace std;



double u1_func(double x, double t) {
    return exp(t) * sinh(x);
}

void z1() {
    double ht, hx;

    double cur_t = 0;

    ofstream error("../../../trs_labs_/output/lab_2/z1_error.txt");
    error << "count" << " " << "error" << " " << "time" << endl;
    vector<int> nX{10,100,130};
    vector<int> nT(nX.size());

    for (int i = 0; i < nX.size(); i++) {
        clock_t start = clock();
        hx = 1. / (nX[i] - 1);
        nT[i] = 2 * (nX[i]) * (nX[i]);
        ht = 1. / (nT[i] - 1);
        vector<vector<double>> M(nT[i], vector<double>(nX[i], 0)); // t , x ;
        ofstream file_to_cout("../../../trs_labs_/output/lab_2/z1_" + to_string(i + 1) + ".txt");
        ofstream file_to_cout_2("../../../trs_labs_/output/lab_2/z1_error" + to_string(i + 1) + ".txt");

        for (int k = 0; k < nX[i]; k++) {
            M[0][k] = sinh(k * hx);
        }
        double gamma = ht / hx / hx;
        for (int j = 1; j < nT[i]; j++) {

            for (int k = 0; k < nX[i]; k++) {
                if (k == 0) {
                    cur_t = j * ht;
                    M[j][k + 1] = gamma * (M[j - 1][k] - 2 * M[j - 1][k + 1] + M[j - 1][k + 2]) + M[j - 1][k + 1];
                    M[j][k] = (hx * exp(cur_t) - M[j][k + 1]) / (hx - 1);
//                    cout << M[j][k] << endl;
                    continue;
                }
                if (k == (nX[i] - 1)) {
                    cur_t = j * ht;
//                    M[j][k - 1] = gamma * (M[j - 1][k - 2] - 2 * M[j - 1][k - 1] + M[j - 1][k]) + M[j - 1][k - 1];
                    M[j][k] = (hx * exp(cur_t + 1) + M[j][k - 1]) / (hx + 1);
                    continue;
                }

                M[j][k] = gamma * (M[j - 1][k - 1] - 2 * M[j - 1][k] + M[j - 1][k + 1]) + M[j - 1][k];
            }
        }

        clock_t end = clock();
        double time = (double) (end - start) / CLOCKS_PER_SEC;

        for (int j = 0; j < nT[i]; j++) {
            for (int k = 0; k < nX[i]; k++) {
                file_to_cout << M[j][k] << " ";
            }
            file_to_cout << endl;
        }

        double dif;
        double max_dif;
        for (int j = 0; j < nT[i]; j++) {
            for (int k = 1; k < nX[i]; k++) {
                max_dif = -1000;
                dif = abs(abs(u1_func(k * hx, j * ht)) - abs(M[j][k]));
                file_to_cout_2 << dif << " ";
//                cout << dif << " ";
                if (dif > max_dif) {
                    max_dif = dif;
                }
            }
            file_to_cout_2 << endl;
//            cout << endl;
        }

        error << nX[i] << " " << max_dif << " " << time << endl;
        cout << "Step x: " << hx << setw(13 + nX.size()) << "Step T: " << ht
             << setw(8 + nX.size()) << "Error: " << max_dif << setw(8 + nX.size()) << "Time: " << time << endl;
    }

    error.close();
}

int main(){
     z1();

     return 0;
}