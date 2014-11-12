#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
using namespace std;

ofstream out("Output.txt");

double function(double x)
{
	return cos(x*x);
}

template <typename T>
void print_vector(vector<T> &vect)
{
	for(int i=0;i<vect.size();++i)
		out << setiosflags(ios::left) << setw(10) << setprecision(5) << vect[i];
	out << endl;
}

vector<double> fill_x(double start, double end, int n)
{
	vector<double> x = vector<double>(n+1);
	double h = (end - start)/n;
	for(int i=0;i<x.size();++i)
		x[i] = start + i*h;
	return x;
}

vector<double> fill_f(vector<double> &x)
{
	vector<double> f = vector<double>(x.size());
	for(int i=0;i<x.size();++i)
		f[i] = function(x[i]);
	return f;
}

vector<double> countSimpleSweep(int number, vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &F)
{
	vector<double> result = vector<double>(number);
	vector<double> alpha = vector<double>(number);
	vector<double> beta = vector<double>(number);
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

double h(vector<double> &x, int i)
{
    return x[i] - x[i-1];
}

void fill_sweep_vectors(vector<double> &a, vector<double> &b, vector<double> &F, vector<double> x, vector<double> f)
{
    double temp = 0.;
    int N = b.size();
    for(int i=1;i<=N;++i)
    {
        temp = 2*(h(x, i) + h(x, i+1));
        a[i-1] = -h(x, i)/temp;
        b[i] = -h(x, i+1)/temp;
        F[i] = 6/temp*((f[i+1] - f[i])/h(x, i+1) - (f[i] - f[i-1])/h(x, i));
    }
    b[0] = 0.;
    F[0] = 0.;
    a[N-1] = 0.;
    F[N] = 0.;
}

double P3(int i, double val, vector<double> &x, vector<double> &f, vector<double> &M)
{
    double result = 0.;
    result = M[i-1]*powf(x[i] - val, 3.)/(6*h(x, i));
    result += M[i]*powf(val - x[i-1], 3.)/(6*h(x, i));
    result +=(f[i-1] - powf(h(x, i), 2.)/6*M[i-1])*(x[i] - val)/h(x, i);
    result += (f[i] - M[i]*powf(h(x, i), 2.)/6*M[i])*(val - x[i-1])/h(x, i);
    return result;
}

int main()
{
	int n = 20;
	double start = 1.;
	double end = 3.;
	vector<double> x = fill_x(start, end, n);
	vector<double> f = fill_f(x);
	print_vector(x);
	print_vector(f);
	vector<double> a = vector<double>(n-1);
	vector<double> b = vector<double>(n-1);
	vector<double> c = vector<double>(n, 1);
	vector<double> F = vector<double>(n);
	fill_sweep_vectors(a, b, F, x, f);
	vector<double> M = countSimpleSweep(n, a, b, c, F);
	print_vector(M);
	vector<double> tmp = vector<double>(n);
	for(int i=1;i<=n;++i)
        tmp[i-1] = P3(i, x[i], x, f, M);
	print_vector(tmp);
	return 0;
}
