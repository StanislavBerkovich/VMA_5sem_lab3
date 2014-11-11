#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

ofstream out("Output.txt");

double function(double x)
{
	return cos(x*x);
}

void print_vector(double* vector, int n)
{
	for(int i=0;i<n;++i)
		out << setiosflags(ios::left) << setw(10) << setprecision(5) << vector[i];
	out << endl;
}

double* fill_x(double start, double end, int n)
{
	double* x = new double[n+1];
	double h = (end - start)/n;
	for(int i=0;i<=n;++i)
		x[i] = start + i*h;
	return x;
}

double* fill_f(double* x, int x_number)
{
	double* f = new double[x_number];
	for(int i=0;i<x_number;++i)
		f[i] = function(x[i]);
	return f;
}

double* countSimpleSweep(int number, double* a, double* b, double* c, double* F)
{
	double* result = new double[number];
	double* alpha = new double[number];
	double* beta = new double[number];
	alpha[0] = b[0]/c[0];
	beta[0] = F[0]/c[0];
	double del;
	for(int k=0;k<number-1;++k)
	{
		del = c[k+1] - alpha[k]*a[k];
		alpha[k+1] = b[k+1]/del;
		beta[k+1] = (F[k+1] + beta[k]*a[k])/del;
	}
	result[number - 1] = beta[number-1];
	for(int i=number-2;i>=0;--i)
		result[i] = alpha[i]*result[i+1] + beta[i];
	return result;
}

int main()
{
	int n = 20;
	double a = 1.;
	double b = 3.;
	double* x = fill_x(a, b, n);
	double* f = fill_f(x, n+1);
	print_vector(x, n+1);
	print_vector(f, n+1);
	return 0;
}