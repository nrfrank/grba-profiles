#include "base.h"

#include <iostream>
#include <cmath>


GrbaProfiles::GrbaProfiles(double THV, double KAP, double SIG, double K, double P,
                           double GA) : thv(THV), kap(KAP), sig(SIG), k(K), p(P), ga(GA), gk((4.0 - k) * ga * ga),
                                        bg((1.0 - p) / 2.0), tan_thv(tan(thv)), tan_thv_sq(tan(thv) * tan(thv)),
                                        sin_2thv(sin(2.0 * thv)), cos_thv(cos(thv)), sin_thv(sin(thv)),
                                        chi_exp((7.0 * k - 23.0 + bg * (13.0 + k)) / (6.0 * (4.0 - k))),
                                        y_exp(0.5 * (bg * (4.0 - k) + 4.0 - 3.0 * k)) {

}

double GrbaProfiles::ThetaPrime(double phi, double r) {
    double cos_phi = cos(phi);
    double numer = r * pow(pow(cos_thv, 2) - 0.25 * pow(sin_2thv, 2) * pow(cos_phi, 2), 0.5);
    double denom = 1.0 + 0.5 * r * sin_2thv * cos_phi;
    return numer / denom;
}

double GrbaProfiles::EnergyProfile(double phi, double r) {
    double thp = ThetaPrime(phi, r);
    return exp2(-pow(thp / sig, 2.0 * kap));
}
