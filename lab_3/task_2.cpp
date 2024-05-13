#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;
#define PI 3.14159

const double epsilon = 1.0e-10;

double f(double x, double y)
{
    return -PI * PI * sin(PI * y) * (5. * cos(PI * x) * cos(PI * x) - 3.);
}

double phi(double x, double y)
{
    return 0;
}

double Solution(double x, double y)
{
    return sin(PI * x) * sin(PI * x) * sin(PI * y);
}

double ChebNormErr(double *y, double lx, double ly, int N)
{
    double hx = lx/N;
    double hy = (ly*hx)/lx;
    double max=-1.0;
    for (int i=0; i<=N; i++)
    {
        for (int j=0; j<=N; j++)
        {
            if (abs(y[j*(N+1)+i]-Solution(hx*i,hy*j))>=max)
                max=abs(y[j*(N+1)+i]-Solution(hx*i,hy*j));
        }
    }
    return max;
}

double delta(double *y_0, double *y, double N)
{
    double max=-1;
    double temp=0;
    for (int i=0; i<(N+1)*(N+1); i++)
    {
        temp=abs(y_0[i]-y[i]);
        if (temp>max)
            max=temp;
    }
    return max;
}

void SOR_Puasson(double *y, double lx, double ly, int N, double w)
{
    double hx = lx/N;
    double hy = (ly*hx)/lx;
    double *y_0 = new double[(N+1)*(N+1)];
    double D = 0;
    int iter = 0;
    for (int i=0; i<=N; i++)
    {
        for (int j=0; j<=N; j++)
        {
            if ((i==0)||(i==N)||(j==0)||(j==N))
            {
                y_0[j*(N+1)+i]=phi(i*hx,j*hy);
                y[j*(N+1)+i]=y_0[j*(N+1)+i];
            }
            else
                y_0[j*(N+1)+i]=0;
        }
    }
    do
    {
        for (int i=1; i<N; i++)
        {
            for (int j=1; j<N; j++)
            {
                y[j*(N+1)+i]=(1-w)*y_0[j*(N+1)+i]+w*
                                                  ((y_0[j*(N+1)+i+1]+y[j*(N+1)+i-1])/pow(hx,2)+
                                                   (y_0[(j+1)*(N+1)+i]+y[(j-1)*(N+1)+i])/pow(hy,2)+f(i*hx,j*hy))/(2.0/pow(hx,2)+2.0/pow(hy,2));
            }
        }
        D = delta(y_0,y,N);
        for (int i=0; i<(N+1)*(N+1); i++)
            y_0[i]=y[i];
        iter++;
    }
    while (D>epsilon);
    cout << iter << endl;
}
int main(int argc, char* argv[])
{
    int N=100;

    double lx=1.0;
    double ly=1.0;

    if (argc == 2)
    {
        N = atoi(argv[1]);
    }

    double hx = 1. / (double)(N);
    double hy = 1. / (double)(N);

    double *y = new double[(N+1)*(N+1)];
    SOR_Puasson(y,lx,ly,N,1.94);
    ofstream filewrite("SOR.txt");
    for (int i=0; i<=N; i++)
    {
        for (int j=0; j<=N; j++)
            filewrite << i*hx << " " << j*hy << " " << y[j*(N+1)+i] << endl;
    }

    cout << ChebNormErr(y,lx,ly,N) << endl;
    return 0;
}