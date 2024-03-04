#define _USE_MATH_DEFINES
#define  _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cstdio>
#include <cmath>

#include <vector>

using namespace std;

double Y0(double x, double y0, double y1) {  return  y1;}
double Y1(double x, double y0, double y1) { return 5.0 + 7.0 * y1 - 12.0 * y0; }
double W(double x, double y0, double y1) { return 7.0 * y1 - 12.0 * y0; }

void Save(char* fname, std::vector<double> &x, std::vector<double> &data) {
	FILE *f = fopen(fname, "w");
	for (size_t i = 0; i < data.size(); i++) {
		fprintf(f, "%15f %20.10lf\n", x[i], data[i]);
	}
	fclose(f);
}


double maple_result(double x)
{
	return ((double)1.0 / (double)3.0)*exp(3 * x) + exp(4 * x) * (double)1.0 / (double)4.0 + (double)5.0 / (double)12.0;
}

void Euler(std::vector<double> &x, std::vector<double> &y0, std::vector<double> &y1) {
	size_t n = y0.size();
	y0[0] = 1; 
	y1[0] = 2;
	for (size_t i = 1; i < n; i++) { 
		y0[i] = y0[i - 1] +(x[i]-x[i-1])*Y0(x[i-1], y0[i-1], y1[i - 1]);
			y1[i] = y1[i - 1] + (x[i] - x[i - 1])*Y1(x[i - 1], y0[i - 1], y1[i - 1]);
	}
	
}
		
void Runge_Kutta(std::vector<double> &x, std::vector<double> &y0, std::vector<double> &y1)
{
	double k1, k2, k3, k4, k11, k22, k33, k44;
	double h =0;
	size_t n = y0.size();
	

	for (size_t i = 0; i < n-1; i++) {
		h = abs(x[i + 1] - x[i]);
		k1 = Y0(x[i], y0[i],y1[i]);
		k11 = Y1(x[i], y0[i], y1[i]);
		k2 = Y0(x[i] + h / 2, y0[i] + (h*k1) / 2, y1[i] + (h*k11) / 2);
		k22 = Y1(x[i] + h / 2, y0[i] + (h*k1) / 2, y1[i] + (h*k11) / 2);
		k3 = Y0(x[i] + h / 2, y0[i] + (h*k2) / 2, y1[i] + (h*k22) / 2);
		k33 = Y1(x[i] + h / 2, y0[i] + (h*k2) / 2, y1[i] + (h*k22) / 2);
		k4 = Y0(x[i] + h, y0[i] + h*k3, y1[i] + h*k33);
		k44 = Y1(x[i] + h, y0[i] + h*k3, y1[i] + h*k33);
		y0[i + 1] = y0[i] + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
		y1[i + 1] = y1[i] + (h / 6)*(k11 + 2 * k22 + 2 * k33 + k44);
		
	}
	
}

void AdamsBasfort(std::vector<double> &x, std::vector<double> &y0, std::vector<double> &y1)
{
	double k1, k2, k3, k4, k11, k22, k33, k44;
	double h = 0;
	size_t n = y0.size();
	y0[0] = 1;
	y1[0] = 2;
	for (size_t i = 0; i < n - 1; i++) {
		h = abs(x[i + 1] - x[i]);
		if (i <= 2)
		{
			k1 = Y0(x[i], y0[i], y1[i]);
			k11 = Y1(x[i], y0[i], y1[i]);;
			k2 = Y0(x[i] + h / 2, y0[i] + (h*k1) / 2, y1[i] + (h*k11) / 2);
			k22 = Y1(x[i] + h / 2, y0[i] + (h*k1) / 2, y1[i] + (h*k11) / 2);
			k3 = Y0(x[i] + h / 2, y0[i] + (h*k2) / 2, y1[i] + (h*k22) / 2);
			k33 = Y1(x[i] + h / 2, y0[i] + (h*k2) / 2, y1[i] + (h*k22) / 2);
			k4 = Y0(x[i] + h, y0[i] + h*k3, y1[i] + h*k33);
			k44 = Y1(x[i] + h, y0[i] + h*k3, y1[i] + h*k33);
			y0[i + 1] = y0[i] + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
			y1[i + 1] = y1[i] + (h / 6)*(k11 + 2 * k22 + 2 * k33 + k44);
		}
		if (i > 2) {
			y0[i + 1] = y0[i] + h*(55 * Y0(x[i], y0[i], y1[i]) - 59 * Y0(x[i - 1], y0[i - 1], y1[i - 1]) + 37 * Y0(x[i - 2], y0[i - 2], y1[i - 2]) - 9 * Y0(x[i - 3], y0[i - 3], y1[i - 3])) / 24.0;
			y1[i + 1] = y1[i] + h*(55 * Y1(x[i], y0[i], y1[i]) - 59 * Y1(x[i - 1], y0[i - 1], y1[i - 1]) + 37 * Y1(x[i - 2], y0[i - 2], y1[i - 2]) - 9 * Y1(x[i - 3], y0[i - 3], y1[i - 3])) / 24.0;
		}
	}
	
}


void Plot(int a) {
	if (a==2)
	system("\"C:\\Program Files\\gnuplot\\bin\\wgnuplot.exe\" plot1.plt");
	if(a==1)
		system("\"C:\\Program Files\\gnuplot\\bin\\wgnuplot.exe\" plot2.plt");
	if (a == 3)
		system("\"C:\\Program Files\\gnuplot\\bin\\wgnuplot.exe\" plot3.plt");
	if (a == 4)
		system("\"C:\\Program Files\\gnuplot\\bin\\wgnuplot.exe\" plot4.plt");
}

void Delta()
{
	int size = 10;
	std::vector<double> h(size+1);
	
	std::vector<double> E_norma(size+1);
	std::vector<double> RK_norma(size+1);
	std::vector<double> AB_norma(size+1);

	for (size_t i = 0; i < E_norma.size(); i++) {
		int n = 100000/pow(2,i);
		
		h[i]= (double)1 / (n - 1);
		std::vector<double> delta(n);
		std::vector<double> y0(n);
		std::vector<double> y1(n);
		std::vector<double> x(n);

		for (int j = 0; j<n; j++)
		{
			x[j] = j*h[i];
		}
		y0[0] = 1;
		y1[0] = 2;
		Euler(x, y0, y1);
		double Sum=0;
		for (size_t j = 0; j <n; j++)
		{
			delta[j] = abs(maple_result(x[j]) - y0[j]);
			if (Sum < delta[j])Sum = delta[j];
		}
		E_norma[i] = log(Sum);

		y0[0] = 1;
		y1[0] = 2;
		Runge_Kutta(x, y0, y1);
		double Sum1=0;
		for (size_t j = 0; j <n; j++)
		{
			delta[j] = abs(abs(maple_result(x[j])) - abs(y0[j]));
			if (Sum1 < delta[j])Sum1 = delta[j];
		}
		RK_norma[i] = log(Sum1);

		y0[0] = 1;
		y1[0] = 2;
		AdamsBasfort(x, y0, y1);
		double Sum2=0;
		for (size_t j = 0; j <n; j++)
		{
			delta[j] = abs(maple_result(x[j])- y0[j]);
			if (Sum2 < delta[j])Sum2 = delta[j];
		}
		AB_norma[i] = log(Sum2);

		h[i] = log(h[i]);
	}
	
	Save("DelEuler.txt", h,E_norma );
	Save("DelRK.txt", h, RK_norma);
	Save("DelAB.txt", h, AB_norma);
}



bool PraviloRunge(std::vector<double> &x, std::vector<double> &y0, std::vector<double> &y02h, int p)
{
	size_t n = y0.size();
	std::vector<double> delta(n);
	std::vector<double> h(n);
	std::vector<double> raznost(n);
	double norma = 0;
	bool a=true;
	for (size_t i = 0; i < n; i++) {
		delta[i] = abs(abs(maple_result(x[i])) - abs(y0[i]));
		raznost[i] = abs(y0[i] - y02h[i * 2]) / (pow(2, p) - 1);
		if (delta[i] < raznost[i])a = false;
	
	}
	
	return a;
}



void Shooting(std::vector<double> &x)
{
	size_t n = x.size();
	double k1, k2, k3, k4, k11, k22, k33, k44, K;
	double h = 0;
	std::vector<double> v(n);
	std::vector<double> dv(n);
	std::vector<double> w(n);
	std::vector<double> dw(n);
	std::vector<double> u(n);
	std::vector<double> Del(n);
	v[0] = 1;
	dv[0] = 0;
	Runge_Kutta(x, v, dv);
	w[0] = 0;
	dw[0] = 1;

	for (size_t i = 0; i < n - 1; i++) {
		h = (x[i + 1] - x[i]);
		k1 = Y0(x[i], w[i], dw[i]);
		k11 = W(x[i], w[i], dw[i]);
		k2 = Y0(x[i] + h / 2, w[i] + (h*k1) / 2, dw[i] + (h*k11) / 2);
		k22 = W(x[i] + h / 2, w[i] + (h*k1) / 2, dw[i] + (h*k11) / 2);
		k3 = Y0(x[i] + h / 2, w[i] + (h*k2) / 2, dw[i] + (h*k22) / 2);
		k33 = W(x[i] + h / 2, w[i] + (h*k2) / 2, dw[i] + (h*k22) / 2);
		k4 = Y0(x[i] + h, w[i] + h*k3, dw[i] + h*k33);
		k44 = W(x[i] + h, w[i] + h*k3, dw[i] + h*k33);
		w[i + 1] = w[i] + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
		dw[i + 1] = dw[i] + (h / 6)*(k11 + 2 * k22 + 2 * k33 + k44);
	}
	
	K = (20.7613831493 - v[n-1]) / w[n-1];
	
	for (size_t i = 0; i < n ; i++) {
		u[i] = v[i] + K*w[i];
		Del[i]= log(abs(u[i] - maple_result(x[i]))) ;
	}

	
	Save("SolKray.txt", x, Del);
	Save("Vshoot.txt", x, v);
	Save("Wshoot.txt", x, w);
	Save("Ushoot.txt", x, u);
}



int main() {
	int n = 100;
	std::vector<double> y0(n);
	std::vector<double> y1(n);
	std::vector<double> x(n);
	std::vector<double> y02h(n*2);
	std::vector<double> y12h(n*2);
	std::vector<double> x2h(n*2);
	double h = (double)1/(n-1);
	for (int i = 0; i<n;i++)
	{
		x[i] = i*h;
		
	}
	h = h/2.0;

	for (int i = 0; i<n*2; i++)
	{
		x2h[i] = i*h;
	}

	
	
	for (size_t i = 0; i < y0.size(); i++) {
		y0[i] = maple_result(x[i]);
		
	}
	Save("Mpl.txt", x, y0);

	for (size_t i = 0; i < n; i++)
	{
		y0[i] = 0;
		y1[i] = 0;
	}

	
	Euler(x,y0,y1);
	Euler(x2h, y02h, y12h);
	cout<<"Pravilo Runge Euler:"<<boolalpha <<PraviloRunge(x, y0, y02h,1)<<endl;
	Save("Euler.txt", x, y0);
	

	for (size_t i = 0; i < n; i++)
	{
	y0[i]=0;
	y1[i] = 0;
	}

	y0[0] = 1;
	y1[0] = 2; 
	y02h[0] = 1;
	y12h[0] = 2;
	Runge_Kutta(x, y0, y1);
	Runge_Kutta(x2h, y02h, y12h);
	cout << "Pravilo Runge RK:" << boolalpha << PraviloRunge(x, y0, y02h, 4)<<endl;
	Save("RK.txt", x, y0);
	

	for (size_t i = 0; i < n; i++)
	{
		y0[i] = 0;
		y1[i] = 0;
	}



	AdamsBasfort(x, y0, y1);
	AdamsBasfort(x2h, y02h, y12h);
	cout << "Pravilo Runge AB:" << boolalpha << PraviloRunge(x, y0, y02h, 4)<<endl;
	Save("AB.txt", x, y0);
	
	Plot(2);
	Delta();
	Plot(1);
	Shooting(x);
	Plot(3);
	Plot(4);
	return 0;
}