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