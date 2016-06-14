#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <stdio.h>
#include "nr3.h"
#include "interp_1d.h"
#include "quadrature.h"
#include "romberg.h"

using namespace std;

const Doub TORAD = M_PI / 180.0;

float integrand(float x)
{
	float val;
	val = pow(x, 4)*log(x + sqrt(pow(x, 2) + 1));
	return val;
}

int addition(int a, int b)
{
	int r;
	r = a + b;
	return r;
}

Doub thetaPrime(const Doub r, const Doub thv, const Doub sig, const Doub phi) {
	Doub numer = r*pow(pow(cos(thv), 2) - 0.25*pow(sin(2.0*thv), 2)*pow(cos(phi), 2), 0.5);
	Doub denom = 1.0 + 0.5*r*sin(2.0*thv)*cos(phi);
	return numer / denom;
}

Doub energyProfile(const Doub thp, const Doub sig, const Doub kap) {
	return exp2(-pow(thp / sig, 2.0*kap));
}

struct RootFunc
{
	Doub r0, kap, thv, sig, phi;
	RootFunc(const Doub R0, const Doub KAP, const Doub THV, const Doub SIG, const Doub PHI) : r0(R0), kap(KAP), thv(THV*TORAD), sig(SIG), phi(PHI*TORAD) {}
	Doub f(const Doub r) {
		const Doub thp = thetaPrime(r, thv, sig, phi);
		const Doub eng = energyProfile(thp, sig, kap);
		const Doub lhs = (pow(r, 2) + 2.0*r*tan(thv)*cos(phi) + pow(tan(thv), 2))*eng;
		const Doub thp0 = thetaPrime(r, thv, sig, 0.0);
		const Doub eng0 = energyProfile(thp0, sig, kap);
		const Doub rhs = pow(r0 + tan(thv), 2)*eng0;
		return lhs - rhs;
	}
	Doub df(const Doub r) {
		const Doub thp = thetaPrime(r, thv, sig, phi);
		const Doub first = r + tan(thv)*cos(phi);
		const Doub second = pow(r, 2) + 2 * r*tan(thv)*cos(phi) + pow(tan(thv), 2);
		const Doub frac = (kap*log(2.0)*pow(thp / sig, 2.0*kap)) / (r*(1.0 + 0.5*r*sin(2.0*thv)*cos(phi)));
		const Doub exponent = 2.0*energyProfile(thp, sig, kap);
		return (first - second*frac)*exponent;
	}
};

template <class T>
Doub rtnewt(T &RootFunc, const Doub g, const Doub xacc) {
	const Int JMAX = 20;
	Doub rtn = g;
	for (Int j = 0; j < JMAX; j++) {
		Doub f = RootFunc.f(rtn);
		Doub df = RootFunc.df(rtn);
		Doub dx = f / df;
		//cout << f << ", " << df << ", " << dx;
		rtn -= dx;
		if (abs(dx) < xacc) {
			cout << "Convergance in " << j << " steps.";
			return rtn;
		}
	}
	throw("Maximum number of iterations exceeded in rtnewt");
}

int main(void)
{
	/*
	RootFunc fx(0.1, 0.0, 0.0, 2.0, 0.0);
	cout << "function at 0.1: " << fx(0.1) << ", derivative at 0.1: " << fx.df(0.1);*/
	Doub R0, KAPPA, THETA_V, SIGMA, PHI, G;
	cout << "Enter values for r0, kappa, thetaV, sigma, phi, and initial r guess.";
	cin >> R0 >> KAPPA >> THETA_V >> SIGMA >> PHI >> G;
	RootFunc fx(R0, KAPPA, THETA_V, SIGMA, PHI);
	Doub root = rtnewt(fx, G, 1.0e-11);
	cout << "The root value is: " << root;
	return 0;
}