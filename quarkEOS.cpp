#include "quarkEOS.hpp"
#include "NumMethods.hpp"
#include <cmath>
#include <iostream>

using namespace std;

const double pi = 4.0*atan(1.0);
const double mE = 0.511;            // Mass of electron (MeV)
const double mU = 2.3;              // Mass of up quark (MeV)
const double mD = 4.8;              // Mass of down quark (MeV)
const double mS = 95.0;            // Mass of strange quark (MeV)
const double mMU = 105.6583745;     // Mass of muon (MeV)

double quark :: pr(double cp, double mass, double deg) {
    double z, pr;
    z = mass/cp;
    if (mass > cp) {
        return 0;
    } else {
        pr = deg*pow(cp,4)/(24.0*pow(pi,2))*( sqrt(1.0-pow(z,2))*(1.0-5.0/2.0*pow(z,2)) + 3.0/2.0*pow(z,4)*log((1.0+sqrt(1.0-pow(z,2)))/z) );
        return pr;
    }
}

double quark :: en(double cp, double mass, double deg) {
    double z, en;
    z = mass/cp;
    if (mass > cp) {
        return 0;
    } else {
        en = deg*pow(cp,4)/(8.0*pow(pi,2))*( sqrt(1.0-pow(z,2))*(1.0-0.5*pow(z,2)) - 0.5*pow(z,4)*log((1.0+sqrt(1.0-pow(z,2)))/z) );
        return en;
    }
}

double quark :: nb(double cp, double mass, double deg) {
    double z, nb;
    z = mass/cp;
    if (mass > cp) {
        return 0;
    } else {
        nb = deg*pow(cp,3)/(6.0*pow(pi,2))*pow((1.0-pow(z,2)),3.0/2.0);
        return nb;
    }
}

double quark :: NRpr(double cp, double mass, double deg) {
    double pr;
    if (mass > cp) {
        return 0;
    } else {
        pr = deg*pow(cp-mass,5.0/2.0)*pow(2.0*mass,3.0/2.0)/(15*pow(pi,2));
        return pr;
    }
}

double quark :: NRen(double cp, double mass, double deg) {
    double en;
    if (mass > cp) {
        return 0;
    } else {
        en = deg*mass*sqrt(2.0*mass)*pow(cp-mass,3.0/2.0)*(3.0*cp+2.0*mass)/(15*pow(pi,2));
        return en;
    }
}

double quark :: NRnb(double cp, double mass, double deg) {
    double nb;
    if (mass > cp) {
        return 0;
    } else {
        nb = deg*mass*sqrt(2.0*mass)*pow(cp-mass,3.0/2.0)/(3*pow(pi,2));
        return nb;
    }
}
