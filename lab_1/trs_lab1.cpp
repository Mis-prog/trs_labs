#include <iostream>
#include "task5.h"
#include "task4.h"

int main() {
    int n = 1e2;
    double x_step;
    t4(n, -10, 10, "../../../trs_labs_/output/lab_1/4.csv", x_step);
    n = 1e4;
    t5(n, -10, 10, "../../../trs_labs_/output/lab_1/5.csv", x_step);

    return 0;
}