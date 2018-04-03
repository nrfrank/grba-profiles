#ifndef GRBA_PROFILES_BASE
#define GRBA_PROFILES_BASE

class GrbaProfiles {
public:
    GrbaProfiles(double THV, double KAP, double SIG, double K, double P, double GA);

    double ThetaPrime(double phi, double r);

    double EnergyProfile(double phi, double r);

    const double thv, kap, sig, k, p, ga, gk, bg, tan_thv;

protected:
    const double tan_thv_sq, sin_2thv, cos_thv, sin_thv, chi_exp, y_exp;
};

#endif