#include "nr3.h"
#include "interp_1d.h"
#include "quadrature.h"
#include "romberg.h"
#include <iostream>
#include <math.h>

using namespace std;

double integrand(double x)
{
	double val;
	val = pow(x, 4)*log(x + sqrt(pow(x, 2) + 1));
	return val;
}

int addition(int a, int b)
{
	int r;
	r = a + b;
	return r;
}

int main()
{
	double result;
	result = qromb(integrand, 0, 2);

	cout << "The result is " << result;
}