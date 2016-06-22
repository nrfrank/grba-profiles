#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "nr3.h"
#include "interp_1d.h"
#include "quadrature.h"
//#include "romberg.h"

using namespace std;

const Doub TORAD = M_PI / 180.0;

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
	const Doub r0, kap, thv, sig, phi;
	RootFunc(const Doub R0, const Doub KAP, const Doub THV, const Doub SIG, const Doub PHI) : 
		r0(R0), kap(KAP), thv(THV*TORAD), sig(SIG), phi(PHI*TORAD) {}
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
Doub rtnewt(T &func, const Doub g, const Doub xacc) {
	const Int JMAX = 20;
	Doub rtn = g;
	for (Int j = 0; j < JMAX; j++) {
		Doub f = func.f(rtn);
		Doub df = func.df(rtn);
		Doub dx = f / df;
		//cout << f << ", " << df << ", " << dx;
		rtn -= dx;
		if (abs(dx) < xacc) {
			//cout << "Convergance in " << j << " steps.";
			//pair <Doub, Int> rtnPair(rtn, j);
			return rtn;
		}
	}
	throw("Maximum number of iterations exceeded in rtnewt");
}

struct paramsPhi
{
	const Doub r0, kap, thv, sig;
	paramsPhi(const Doub R0, const Doub KAP, const Doub THV, const Doub SIG, const Doub PHI) :
		r0(R0), kap(KAP), thv(THV*TORAD), sig(SIG) {}
};

struct integrandPhi
{
	const Doub r0, kap, thv, sig, phi;
	integrandPhi(const Doub R0, const Doub KAP, const Doub THV, const Doub SIG, const Doub PHI) :
		r0(R0), kap(KAP), thv(THV*TORAD), sig(SIG), phi(PHI*TORAD) {}
	Doub operator() (const Doub x) {
		return pow(sin(x), 2);
	}
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

struct QuadraturePhi {
	Int n;
	virtual Doub next() = 0;
};

template<class T>
struct TrapzdPhi : QuadraturePhi {
	Doub a, b, s;
	T &func;
	TrapzdPhi() {};
	TrapzdPhi(T &funcc, const Doub aa, const Doub bb) :
		func(funcc), a(aa), b(bb) {
		n = 0;
	}
	Doub next() {
		Doub x, tnm, sum, del;
		Int it, j;
		n++;
		FILE * ofile;
		ofile = fopen("phi-integration-test.txt", "a+");
		if (n == 1) {
			return (s = 0.5*(b - a)*(func(a) + func(b)));
		}
		else {
			for (it = 1, j = 1; j<n - 1; j++) it <<= 1;
			tnm = it;
			del = (b - a) / tnm;
			x = a + 0.5*del;
			for (sum = 0.0, j = 0; j < it; j++, x += del) {
				sum += func(x);
				fprintf(ofile, "%d\t%d\t%f\t%f\n", n, j, x, sum);
			}
			fclose(ofile);
			s = 0.5*(s + (b - a)*sum / tnm);
			return s;
		}
	}
};

template<class T>
Doub qsimpPhi(T &func, const Doub a, const Doub b, const Doub eps = 1.0e-10) {
	const Int JMAX = 20;
	Doub s, st, ost = 0.0, os = 0.0;
	TrapzdPhi<T> t(func, a, b);
	for (Int j = 0; j<JMAX; j++) {
		st = t.next();
		s = (4.0*st - ost) / 3.0;
		if (j > 5)
			if (abs(s - os) < eps*abs(os) ||
				(s == 0.0 && os == 0.0)) return s;
		os = s;
		ost = st;
	}
	throw("Too many steps in routine qsimp");
}

int main(void)
{
	
	Doub R0, KAPPA, THETA_V, SIGMA, PHI, G;
	R0 = 0.1;
	KAPPA = 1.0;
	THETA_V = 6.0;
	SIGMA = 2.0;
	PHI = 5.0;
	G = R0;
	/*
	char filename[50];
	FILE * ofile;
	Doub YVAL = 0.5, RMAX;
	cout << "Enter values for R0, KAPPA, and THETA_V: \n";
	cin >> R0 >> KAPPA >> THETA_V;
	sprintf(filename, "phiRoot_r0=%f_kap=%f_thv=%f.txt", R0, KAPPA, THETA_V);
	ofile = fopen(filename, "w");
	fprintf(ofile, "PHI\tR\tNSTEPS\n");
	for (int p = 0; p <= 500; p++) {
		PHI = 360.0*p / 500.0;
		RootFunc fx(R0, KAPPA, THETA_V, SIGMA, PHI);
		pair <Doub, Int> root = rtnewt(fx, G, 1.0e-11);
		fprintf(ofile, "%f\t%f\t%d\n", PHI, root.first, root.second);
		G = root.first;
	}
	fclose(ofile);
	*/

	// Test out modifications to standard Trapzd and qsimp routines.
	FILE * ofile;
	ofile = fopen("phi-integration-test.txt", "w");
	fprintf(ofile, "n\tj\tx\ts\n");
	fclose(ofile);
	integrandPhi intFunc(R0, KAPPA, THETA_V, SIGMA, PHI);
	Doub intVal = qsimpPhi(intFunc, 0.0, 360.0);
	
	//RootFunc rFunc(R0, KAPPA, THETA_V, SIGMA, PHI);
	Doub rootVal = rtnewt(intFunc, G, 1.0e-11);
	cout << intVal << "\t" << rootVal;

	// Run a test to better understand the Trapzd method.
	/*Int n, it, j;
	Doub x, tnm, del, a = 0.0, b = 360.0;
	cout << "Enter maximum refinement level (n): \n";
	cin >> n;
	for (it = 1, j = 1; j < n - 1; j++) {
		it <<= 1;
		printf("j = %d, it = %d\n", j, it);
	}
	tnm = it;
	del = (b - a) / tnm;
	printf("tnm = %f, del = %f\n", tnm, del);
	x = a + 0.5*del;
	for (j = 0; j < it; j++, x += del) {
		printf("j = %d, x = %f\n", j, x);
	}*/

	return 0;
}