#include "quarkEOS.hpp"
#include "NumMethods.hpp"
#include "nucEOS.hpp"
#include "Conversions.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

const double pi = 4.0*atan(1.0);
const double mE = 0.511;            // Mass of electron (MeV)
const double mU = 2.3;              // Mass of up quark (MeV)
const double mD = 4.8;              // Mass of down quark (MeV)
const double mS = 95.0;            // Mass of strange quark (MeV)
const double mMU = 105.6583745;     // Mass of muon (MeV)
const double mP = 939.0; //938.27231;    // Mass of proton (MeV)
const double mN = 939.0; //939.56542052;
const double mNuc = (mP+mN)/2.0;

quark qk;
Convert conv1;
data2 dm1;
nummeth nm2;

// Relativistic Nuclear Field Theory 
// FSU Model
/*
void EOS :: whitedwarf(double nucNum) {
    double en, pr, mun, nb, mue, ke;
    //double K;
    mun = mN+mE+0.2;
    ofstream out("whitedwarf" + to_string(int(nucNum)) + ".txt");
    double mev4 = conv1.energyCONV(6,0); // mev4 to mev/fm3
    while (mun<1500) {
        mue = 0.5*(pow(mE,2.0) - pow(mN,2.0) + pow(mun,2.0))/mun;
        nb = qk.nb(mun,mN,2.0) + qk.nb(mun-mue,mP,2.0);
        pr = qk.pr(mue,mE,2.0);
        en = qk.nb(mue,mE,2)*mNuc*nucNum + qk.en(mue,mE,2);
        //cout << mun << "  " << mue << "  " << pr << endl;
        out << scientific << setprecision(20) << nb*mev4 << "  " << en*mev4 << "  " << pr*mev4 << endl;
        mun = mun + 0.001;
    }
}
*/

// Get the EOS for a white dwarf
void EOS :: whitedwarf(double nucNum) {
    double en, pr, mue;
    //double K;
    mue = mE+1e-3;
    ofstream out("whitedwarf" + to_string(int(nucNum)) + ".txt");
    double mev4 = conv1.energyCONV(6,0); // mev4 to mev/fm3
    while (mue<25) {
        en = qk.nb(mue,mE,2)*mNuc*nucNum + qk.en(mue,mE,2);
        pr = qk.pr(mue,mE,2);

        // Polytropes

        //K = pow(3*pow(pi,2)/(mNuc*nucNum),polyindex)/(12*pow(pi,2));
        //en = qk.nb(mue,mE,2)*mNuc*nucNum + qk.en(mue,mE,2);
        //pr = K*pow(en,polyindex);

        //cout << mue << "  " << en << "  " << pr << endl;
        out << scientific << setprecision(20) << mue << "  " << qk.nb(mue,mE,2)*mev4 << "  " << en*mev4 << "  " << pr*mev4 << endl;
        mue = mue +1e-3;
    }
}

/*
void EOS :: fermigas() {
    double en, pr, mub;
    mub = mN+0.1;
    ofstream out("fermigas.txt");
    double mev4 = conv1.energyCONV(6,0); // mev4 to mev/fm3
    while (mub<2000) {
        en = qk.en(mub,mN,2);
        pr = qk.pr(mub,mN,2);

        //cout << mue << "  " << en << "  " << pr << endl;
        out << scientific << setprecision(10) << mub << "  " << qk.nb(mub,mN,2)*mev4 << "  " << en*mev4 << "  " << pr*mev4 << endl;
        mub = mub + 1e-1;
    }
}
*/

// Get the energy density for the FSU model
double EOS :: get_en(double kf, double t, double gss, double gsoms2, double gww, double gwomw2, double gpp, double gpomp2, double b, double c, double h, double lambda) {
    double integralp = 0; double integraln = 0; double res,mss2,mww2,mpp2,mstrp,mstrn,kfn,kfp,rn,rp;

    mss2 = pow(gss,2.0)/gsoms2;  // (m_sigma*sigma)^2
    mww2 = pow(gww,2.0)/gwomw2;  // (m_omgea*omega)^2
    mpp2 = pow(gpp,2.0)/gpomp2;  // (m_rho*rho)^2
    mstrp = mP-gss; mstrn = mN-gss;   // effective masses
    kfn = kf*pow(1.0+t,1.0/3.0); kfp = kf*pow(1.0-t,1.0/3.0);     // define the fermi momentum for the proton and neutron given t
    rn = sqrt(pow(kfn,2) + pow(mstrn,2)); rp = sqrt(pow(kfp,2) + pow(mstrp,2));       // convenient quanitity to define since used a lot

    // integrals for the energy density
    integralp = kfp*pow(rp,3)/4.0 - pow(mstrp,2)*kfp*rp/8.0 - pow(mstrp,4)/8.0*log(kfp+rp) + pow(mstrp,4)/8.0*log(abs(mstrp));
    integraln = kfn*pow(rn,3)/4.0 - pow(mstrn,2)*kfn*rn/8.0 - pow(mstrn,4)/8.0*log(kfn+rn) + pow(mstrn,4)/8.0*log(abs(mstrn));
    
    // add free contributions
    res = 0.5*mss2 + 0.5*mww2 + +0.5*mpp2 + 1.0/3.0*b*mNuc*pow(gss,3.0) + 0.25*c*pow(gss,4.0) + 0.75*h*pow(gww,4.0) + 1.5*lambda*pow(gpp,2.0)*pow(gww,2.0)
            + 1.0/pow(pi,2)*integralp + 1.0/pow(pi,2)*integraln;
    return res;
}

// Get the pressure for the FSU model
double EOS :: get_pr(double kf, double t, double gss, double gsoms2, double gww, double gwomw2, double gpp, double gpomp2, double b, double c, double h, double lambda) {
    double integralp = 0; double integraln = 0; double res,mss2,mww2,mpp2,mstrp,mstrn,kfn,kfp,rn,rp;

    mss2 = pow(gss,2.0)/gsoms2;  // (m_sigma*sigma)^2
    mww2 = pow(gww,2.0)/gwomw2;  // (m_omgea*omega)^2
    mpp2 = pow(gpp,2.0)/gpomp2;  // (m_rho*rho)^2
    mstrp = mP-gss; mstrn = mN-gss;   // effective masses
    kfn = kf*pow(1.0+t,1.0/3.0); kfp = kf*pow(1.0-t,1.0/3.0);     // define the fermi momentum for the proton and neutron given t
    rn = sqrt(pow(kfn,2) + pow(mstrn,2)); rp = sqrt(pow(kfp,2) + pow(mstrp,2));   // convenient quanitity to define since used a lot

    // integrals for the pressure
    integralp = pow(kfp,3)*rp - 3.0/4.0*kfp*pow(rp,3) + 3.0/8.0*pow(mstrp,2)*kfp*rp + 3.0/8.0*pow(mstrp,4)*log(kfp+rp) - 3.0/8.0*pow(mstrp,4)*log(abs(mstrp));
    integraln = pow(kfn,3)*rn - 3.0/4.0*kfn*pow(rn,3) + 3.0/8.0*pow(mstrn,2)*kfn*rn + 3.0/8.0*pow(mstrn,4)*log(kfn+rn) - 3.0/8.0*pow(mstrn,4)*log(abs(mstrn));
    
    // add free contributions
    res = -0.5*mss2 - 1.0/3.0*b*mNuc*pow(gss,3.0) - 0.25*c*pow(gss,4.0) + 0.5*mww2 + 0.5*mpp2 + 0.25*h*pow(gww,4.0) + 0.5*lambda*pow(gpp,2.0)*pow(gww,2.0)
        + 1.0/(3.0*pow(pi,2))*integralp + 1.0/(3.0*pow(pi,2))*integraln;
    return res;
}

// Calculate the compressibility 
double EOS :: get_K(double kf, double gss, double gwomw2, double gsoms2, double b, double c, double h, double gww) {
    double integralp = 0; double integraln = 0; double integral = 0; double res, mstrp,mstrn,rn,rp;

    mstrp = mP-gss; mstrn = mN-gss;   // effective masses
    rn = sqrt(pow(kf,2) + pow(mstrn,2)); rp = sqrt(pow(kf,2) + pow(mstrp,2));     // convenient quanitity to define since used a lot

    integralp = (0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstrp,2.0))/rp + 3.0*pow(mstrp,2.0)/2.0*log(abs(mstrp/(kf+rp)));
    integraln = (0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstrn,2.0))/rn + 3.0*pow(mstrn,2.0)/2.0*log(abs(mstrn/(kf+rn)));
    integral = integralp + integraln;
    res = gwomw2/(1.0 + 3.0*gwomw2*h*pow(gww,2.0))*6.0*pow(kf,3)/pow(pi,2) + 3.0/2.0*pow(kf,2.0)*(1.0/rp + 1.0/rn) 
            - 3.0/2.0*gsoms2*pow(kf,3.0)/pow(pi,2.0)*pow(mstrp/rp + mstrn/rn,2.0)*pow(1.0+gsoms2*(2.0*b*mNuc*gss + 3.0*c*pow(gss,2.0) + 1.0/pow(pi,2.0)*integral),-1.0);
    return res;
}

// Calculate the symmetry energy at saturation (J)
double EOS :: get_J(double kf, double gss, double gpomp2, double gww, double lambda) {
    double res,mstarp,mstarn,rn,rp;
    
    mstarp = mP - gss; mstarn = mN - gss; // effective masses
    rn = sqrt(pow(kf,2) + pow(mstarn,2)); rp = sqrt(pow(kf,2) + pow(mstarp,2));   // convenient quanitity to define since used a lot
    res = pow(kf,2.0)/12.0*(1.0/rp + 1.0/rn) + 1.0/(12.0*pow(pi,2.0))*gpomp2/(1.0+gpomp2*pow(gww,2.0)*lambda)*pow(kf,3.0);
    return res;
}

// Calculate the drivative of the symmetry energy at saturatiob (L)
double EOS :: get_L(double kf, double gss, double gsoms2, double gww, double gwomw2, double gpomp2, double h, double lambda, double b, double c) {
    double res, mstarp,mstarn,rn,rp,intp,intn,integral,gsdsdk,gwdwdk;

    mstarp = mP - gss; mstarn = mN - gss; // effective masses
    rn = sqrt(pow(kf,2) + pow(mstarn,2)); rp = sqrt(pow(kf,2) + pow(mstarp,2));   // convenient quanitity to define since used a lot
    intp = (0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstarp,2.0))/rp + 3.0*pow(mstarp,2.0)/2.0*log(abs(mstarp/(kf+rp)));    // proton integral
    intn = (0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstarn,2.0))/rn + 3.0*pow(mstarn,2.0)/2.0*log(abs(mstarn/(kf+rn)));    // neutron integral
    integral = intp + intn;

    // g_sigma*dsigma/dk
    gsdsdk = 2.0/pow(pi,2.0)*pow(kf,2.0)* gsoms2/2.0*(mstarp/rp + mstarn/rn)*pow(1.0+ gsoms2*(2.0*b*mNuc*gss + 3.0*c*pow(gss,2.0) + 1.0/pow(pi,2.0)*integral),-1.0);
    
    // g_omega*domega/dk
    gwdwdk = gwomw2*2.0*pow(kf/pi,2.0)/(1.0 + 3.0*gwomw2*h*pow(gww,2.0));
    
    res = pow(kf,2.0)/6.0*(1.0/rp+1.0/rn) + pow(kf,3.0)/12.0*( (mstarp*gsdsdk - kf)/pow(rp,3.0) + (mstarn*gsdsdk - kf)/pow(rn,3.0) ) 
        + pow(kf,3.0)/(4.0*pow(pi,2.0))*gpomp2/(1.0+gpomp2*pow(gww,2.0)*lambda) - pow(kf,4.0)/(6.0*pow(pi,2.0))*pow(gpomp2,2.0)*gww*lambda*gwdwdk/pow(1.0 + gpomp2*lambda*pow(gww,2.0),2.0);
    return res;
}

// Calculate (g_rho/m_rho)^2 Analytically
double get_gpomp2(double kf, double asym, double L, double gss, double gsoms2, double gww, double gwomw2, double h, double b, double c) {
    double gpomp2,z1,z2, mstarp,mstarn,rn,rp,intp,intn,integral,gsdsdk,gwdwdk,alpha,beta,ap,an,gamma;
    mstarp = mP - gss; mstarn = mN - gss;   // effective masses

    // convenient quantities to define (makes expressions easier to read)
    rn = sqrt(pow(kf,2) + pow(mstarn,2)); rp = sqrt(pow(kf,2) + pow(mstarp,2));
    z2 = pow(kf,3.0)/(12.0*pow(pi,2.0)); z1 = pow(kf,2.0)/12.0*(1.0/rp + 1.0/rn);

    // integrals for proton and neutron
    intp = (0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstarp,2.0))/rp + 3.0*pow(mstarp,2.0)/2.0*log(abs(mstarp/(kf+rp)));
    intn = (0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstarn,2.0))/rn + 3.0*pow(mstarn,2.0)/2.0*log(abs(mstarn/(kf+rn)));
    integral = intp + intn;

    // g_sigma*dsigma/dk and g_omega*d_omega/dk
    gsdsdk = 2.0/pow(pi,2.0)*pow(kf,2.0)* gsoms2/2.0*(mstarp/rp + mstarn/rn)*pow(1.0+ gsoms2*(2.0*b*mNuc*gss + 3.0*c*pow(gss,2.0) + 1.0/pow(pi,2.0)*integral),-1.0);
    gwdwdk = gwomw2*2.0*pow(kf/pi,2.0)/(1.0 + 3.0*gwomw2*h*pow(gww,2.0));

    // more convenient quantities to define (no physical meaning)
    alpha = 1.0/z2*pow(rn*rp,3.0)*gwdwdk*pow(kf,4.0)*pow(asym-z1,2.0);
    beta = asym*pow(rp*rn*kf,3.0)*(-1.5*gww+gwdwdk*kf);
    ap = pow(rp,3.0)*gww*z2*pow(pi*kf,2.0)*(-pow(rn,2.0) + 0.5*pow(kf,2.0) - 0.5*kf*mstarn*gsdsdk);
    an = pow(rn,3.0)*gww*z2*pow(pi*kf,2.0)*(-pow(rp,2.0) + 0.5*pow(kf,2.0) - 0.5*kf*mstarp*gsdsdk);
    gamma = pow(rp*rn,3.0)*(1.5*gww*pow(kf,3.0)*z1 - gwdwdk*pow(kf,4.0)*z1 + 6.0*pow(pi,2.0)*gww*L*z2);
    
    gpomp2 = alpha/(beta+ap+an+gamma);
    return gpomp2;
}

// Field Equation for g_sigma sigma (used for bisection method to solve for g_sigma*sigma)
double gss_FE(double kf, double gss, double gsoms2, double b, double c, double t) {
    double integral = 0;
    double res,mstrp,mstrn,kfn,kfp,rn,rp;
    mstrp = mP-gss; mstrn = mN-gss;     // effective masses
    kfn = kf*pow(1.0+t,1.0/3.0); kfp = kf*pow(1.0-t,1.0/3.0);   // fermi momenta of neutron and proton
    rn = sqrt(pow(kfn,2) + pow(mstrn,2)); rp = sqrt(pow(kfp,2) + pow(mstrp,2));     // convenient quantities to define
    integral = gsoms2/pow(pi,2)*  (mstrp*(kfp*rp/2.0 - pow(mstrp,2)/2.0*log(kfp+rp) + pow(mstrp,2)/2.0*log(abs(mstrp)) ) +
                 mstrn*(kfn*rn/2.0 - pow(mstrn,2)/2.0*log(kfn+rn) +pow(mstrn,2)/2.0*log(abs(mstrn)) )  );
    
    res = integral - b*mNuc*gsoms2*pow(gss,2.0) - c*gsoms2*pow(gss,3.0) - gss;
    return res;
}

// derivative of the sigma field equation (for use in newtons method)
double gssFE_derivative(double kf, double gss, double gsoms2, double b, double c, double t) {
    double integrald = 0;
    double res, mstr, kfn,kfp, rn,rp;

    mstr = mNuc-gss;    // get effective mass
    kfn = kf*pow(1.0+t,1.0/3.0); kfp = kf*pow(1.0-t,1.0/3.0);   // get fermi momenta
    rn = sqrt(pow(kfn,2) + pow(mstr,2)); rp = sqrt(pow(kfp,2) + pow(mstr,2));   // covenient quantities to define
    integrald = gsoms2/pow(pi,2)*(-0.5*kfn*rn-0.5*kfp*rp+pow(mstr,2.0)* (0.5*log(kfn+rn) + 0.5*log(kfp+rp) - 3.0*log(abs(mstr)) - 1.0
                - 0.5*kfn/rn + 0.5*pow(mstr,2.0)/(kfn*rn+rn*rn) + log(kfn+rn) - 0.5*kfp/rp + 0.5*pow(mstr,2.0)/(kfp*rp+rp*rp) + log(kfp+rp) )); // derivative of the integral expression
    
    res = integrald - 2.0*b*mNuc*gsoms2*gss - 3.0*c*gsoms2*pow(gss,2.0) - 1.0;
    return res;
}

// obtain coupling constants given bulk properties with optional terminal output (print=true)
// flag is for MCMC calculations that tells you if the corresponding bulk properties yield an unrealistic EOS since gw/mw becomes imaginary
// Fills a coupling array of the form [ (gs/ms)^2 , (gw/mw)^2 , (gp/mp)^2, b , c , h , lambda ]
int EOS :: get_parameters(double BA, double p0, double J, double mstar, double K, double L, double h, double params[7], bool print, bool flag) {
    double a1,a2,a3,b1,b2,b3,c1,c2,c3,g1,g2,g3,z1,z2,m1,m2,m3,m4,n1,n2,n3,n4;
    double kf,gss,mstarp,mstarn,gww,gwomw2,rn,rp,sdensn,sdensp,sdens,gsoms2,b,c,gpomp2,gpp,lambda,pintegralK,nintegralK,integralK;

    p0 = p0*pow(197.32698,3); // convert density from 1/fm3 to MeV^3;
    kf = pow(3.0/2.0*pow(pi,2)*p0,1.0/3.0); // get fermi momentum

    gss = mNuc - mstar; // given mstar get gss
    mstarp = mP - gss;  mstarn = mN - gss;  // get effective masses (technically redundant since mN = mP)
    gww = mNuc + BA - 0.5*sqrt(pow(kf,2) + pow(mstarp,2)) - 0.5*sqrt(pow(kf,2) + pow(mstarn,2));    // get gww at saturation
    rn = sqrt(pow(kf,2) + pow(mstarn,2)); rp = sqrt(pow(kf,2) + pow(mstarp,2)); // convenient quantities to define
    
    // check to see if the coupling constants are realistic 
    if (p0/pow(gww,3.0) < h) {
        cout << "unrealistic: limit is: " << p0/pow(gww,3.0) << " input is: " << h << "  " << gww << endl;
        if (flag == true) {
            exit(0);
        }
    }

    gwomw2 = gww/(p0-h*pow(gww,3.0));   // get (gw/mw)^2

    // proton and neutron integrals
    pintegralK = (0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstarp,2.0))/rp + 3.0*pow(mstarp,2.0)/2.0*log(abs(mstarp/(kf+rp)));
    nintegralK = (0.5*pow(kf,3.0) + 3.0/2.0*kf*pow(mstarn,2.0))/rn + 3.0*pow(mstarn,2.0)/2.0*log(abs(mstarn/(kf+rn)));
    integralK = pintegralK + nintegralK;

    // algebraically convenient expressions
    n1 = 6.0*pow(kf,3.0)/pow(pi,2.0)*gwomw2/(1.0+ 3.0*gwomw2*h*pow(gww,2.0));
    n2 = 3.0/2.0*pow(kf,2.0)*(1.0/rp + 1.0/rn);
    n3 = 3.0/2.0*pow(kf,3.0)/pow(pi,2.0)*pow(mstarp/rp + mstarn/rn,2.0);    
    n4 = 1.0/pow(pi,2.0)*integralK;
    a1 = gss;
    a2 = mNuc*pow(gss,2.0);
    a3 = pow(gss,3.0);
    b1 = K - n1 - n2;
    b2 = 2.0*mNuc*gss*(K-n1-n2);
    b3 = 3.0*pow(gss,2.0)*(K-n1-n2);
    g1 = 0.5*pow(gss,2.0);
    g2 = 1.0/3.0*mNuc*pow(gss,3.0);
    g3 = 0.25*pow(gss,4.0);
    m1 = 1.0/pow(pi,2.0)* (kf*pow(rp,3)/4.0 - pow(mstarp,2)*kf*rp/8.0 - pow(mstarp,4)/8.0*log(kf+rp) + pow(mstarp,4)/8.0*log(abs(mstarp)));
    m2 = 1.0/pow(pi,2.0)* (kf*pow(rn,3)/4.0 - pow(mstarn,2)*kf*rn/8.0 - pow(mstarn,4)/8.0*log(kf+rn) + pow(mstarn,4)/8.0*log(abs(mstarn)));
    m3 = 0.5*pow(gwomw2,-1.0)*pow(gww,2.0);
    m4 = 0.75*h*pow(gww,4.0);

    // scalar densities
    sdensp = 1.0/pow(pi,2)*mstarp*(kf*rp/2.0 - pow(mstarp,2)/2.0*log(abs(kf+rp)) + pow(mstarp,2)/2.0*log(abs(mstarp)) );
    sdensn = 1.0/pow(pi,2)*mstarn*(kf*rn/2.0 - pow(mstarn,2)/2.0*log(abs(kf+rn)) + pow(mstarn,2)/2.0*log(abs(mstarn)) );
    sdens = sdensp + sdensn;

    // algebraically convenient expressions continued
    c1 = sdens;
    c2 = -n3 - n4*(K-n1-n2);
    c3 = p0*(mNuc+BA)-m1-m2-m3-m4;

    // algebraic expressions for the sigma coupling constants
    gsoms2 = (a3*b2*g1-a2*b3*g1-a3*b1*g2+a1*b3*g2+a2*b1*g3-a1*b2*g3)/(c3*a3*b2-c3*a2*b3-c2*a3*g2+c1*b3*g2+c2*a2*g3-c1*b2*g3);
    b = (c3*a3*b1-c3*a1*b3-c2*a3*g1+c1*b3*g1+c2*a1*g3-c1*b1*g3)/(-a3*b2*g1+a2*b3*g1+a3*b1*g2-a1*b3*g2-a2*b1*g3+a1*b2*g3);
    c = (c3*a2*b1-c3*a1*b2-c2*a2*g1+c1*b2*g1+c2*a1*g2-c1*b1*g2)/(a3*b2*g1-a2*b3*g1-a3*b1*g2+a1*b3*g2+a2*b1*g3-a1*b2*g3);
    
    // get(gp/mp)^2
    gpomp2 = get_gpomp2(kf,J,L,gss,gsoms2,gww,gwomw2,h,b,c);
    gpp = 0.0;   // valid at saturation

    z2 = pow(kf,3.0)/(12.0*pow(pi,2.0));
    z1 = pow(kf,2.0)/12.0*(1.0/rp + 1.0/rn);
    lambda = 1.0/(gpomp2*pow(gww,2.0))*(z2*gpomp2/(J-z1) - 1.0);  // calculate last remaining coupling

    if (print == true) {
        cout << "gsoms2: " << setprecision(10) << gsoms2 << endl;
        cout << "gwomw2: " << setprecision(10) << gwomw2 << endl;
        cout << "gpmop2: " << setprecision(10) << gpomp2 << endl;
        cout << "b: " << setprecision(10) << b << endl;
        cout << "c: " << setprecision(10) << c << endl;
        cout << "h: " << setprecision(10) << h << endl;
        cout << "lambda: " << setprecision(10) << lambda << endl;
        cout << "mstar/m: " << mstar/mNuc << endl;
        cout << "B/A: " << get_en(kf,0,gss,gsoms2,gww,gwomw2,gpp,gpomp2,b,c,h,lambda)/p0 - mNuc << endl;
        cout << "J: " << get_J(kf,gss,gpomp2,gww,lambda) << endl;
        cout << "L: " << get_L(kf,gss,gsoms2,gww,gwomw2,gpomp2,h,lambda,b,c) << endl;;
        cout << "K: " << get_K(kf,gss,gwomw2,gsoms2,b,c,h,gww) << endl; // create function to calulcate K using derivatives;
    }

    params[0] = gsoms2;
    params[1] = gwomw2;
    params[2] = gpomp2;
    params[3] = b;
    params[4] = c;
    params[5] = h;
    params[6] = lambda;

    return 0;
}

// get the EOS for a given number of points and set of couplings
// coupling array of the form [ (gs/ms)^2 , (gw/mw)^2 , (gp/mp)^2, b , c , h , lambda ]
// will output an array of the form [nb fm-3, en MeV/fm3, pr MeV/fm3, dpde, mub MeV, dP/drho, Keff MeV]
int EOS :: get_EOS(double params[7], double** &eos, int npoints, bool print, bool unstable) {
    double k,en,pr,dens,mue,gww,gss,gpp,mstarp,mstarn,kp,kn,mun,mup,t,check,conv_mev4,p0f,ssize,gsoms2,gwomw2,gpomp2,b,c,h,lambda;
    double dydx1,dydx2;
    conv_mev4 = conv1.energyCONV(6,0); // mev4 to mev/fm3
    k = 40.0;   // initial fermi momentum
    p0f = 8.0; // N TIMES SATURATION DENSITY

    // step size that increases with increasing density
    ssize = (-10.4687 + 55.5001*pow(p0f,7.0/25))/npoints;

    // Import the coupling constants from array
    gsoms2 = params[0];
    gwomw2 = params[1];
    gpomp2 = params[2];
    b = params[3];
    c = params[4];
    h = params[5];
    lambda = params[6];
    cout << gpomp2 << "  " << lambda << endl;

    // create EOS array
    eos = new double*[npoints];
    for (int i=0; i<npoints; i++) {
        eos[i] = new double[7];
    }
    
    // return the EOS for a specified number of points
    for (int i=0; i<npoints; ++i) {
        dens = 2.0*pow(k,3)/(3.0*pow(pi,2));    // density at given fermi momentum
        t = get_t_betaeq(k,gsoms2,gwomw2,gpomp2,b,c,lambda,h); // get t from charge neutrality and beta equil
        kp = k*pow(1.0-t,1.0/3.0); kn = k*pow(1.0+t,1.0/3.0);   // get proton and neutron fermi momenta
        gss = FSUgssnewton(dens,gsoms2,b,c,t,1e-7);     // get the sigma field
        gww = FSUgwwnewton(gwomw2,gpomp2,lambda,h,dens,t);  // get omega field
        gpp = -0.5*gpomp2*dens*t/(1.0+lambda*gpomp2*pow(gww,2.0));  // get rho field
        //cout << setprecision(10) << dens << "  " << t << "  " << 0 << "  " << gss << "  " << gww << "  " << gpp << endl;
        mstarp = mP - gss; mstarn = mN - gss; // get efffecitve masses
        
        // get the chemical potentials
        mup = sqrt(pow(kp,2) + pow(mstarp,2)) + gww + 0.5*gpp;  
        mun = sqrt(pow(kn,2) + pow(mstarn,2)) + gww - 0.5*gpp;
        mue = mun - mup;

        // get the EOS
        en = get_en(k,t,gss,gsoms2,gww,gwomw2,gpp,gpomp2,b,c,h,lambda) + qk.en(mue,mE,2.0) + qk.en(mue,mMU,2.0);
        pr = get_pr(k,t,gss,gsoms2,gww,gwomw2,gpp,gpomp2,b,c,h,lambda) + qk.pr(mue,mE,2.0) + qk.pr(mue,mMU,2.0);
        eos[i][0] = dens*conv_mev4;
        eos[i][1] = conv_mev4*en;
        eos[i][2] = conv_mev4*pr;
        eos[i][4] = mun;
        eos[i][5] = 0.5*sqrt(pow(kp,2.0)+pow(mstarp,2.0))*pow(kp/k,2.0)*pow(1.0-t,1.0/3.0) + 0.5*sqrt(pow(k,2.0)+pow(mstarp,2.0))*pow(kn/k,2.0)*pow(1.0+t,1.0/3.0) + gww;   // de/drho
        eos[i][6] = effK2(k,gss,gsoms2,gww,gwomw2,gpomp2,gpp,h,lambda,b,c,t);
        
        // check for thermodynamic stability
        check = mun*dens - en - pr;
        if (abs(check)>100.0) {
            cout << "consistency fail for k: " << k << " with value: " << check << endl;
            exit(0);
        }

        k = k + ssize*log(1.0+k);
    }

    // get the speed of sound by taking derivatives
    eos[0][3] = (eos[1][2] - eos[0][2])/(eos[1][1] - eos[0][1]);
    eos[1][3] = (eos[2][2] - eos[0][2])/(eos[2][1] - eos[0][1]);
    eos[0][5] = eos[0][3]*eos[0][5];
    eos[1][5] = eos[1][3]*eos[1][5];
    
    for (int i=2; i<(npoints-2); ++i) {
        dydx1 = (eos[i+1][2] - eos[i-1][2])/(eos[i+1][1] - eos[i-1][1]);
        dydx2 = (eos[i+2][2] - eos[i-2][2])/(eos[i+2][1] - eos[i-2][1]);
        eos[i][3] = 0.5*(dydx1 + dydx2);
        eos[i][5] = eos[i][3]*eos[i][5];
        if (eos[i][4] < 0) {
            unstable = true;
        }
    }
    
    eos[npoints-2][3] = (eos[npoints-1][2] - eos[npoints-3][2])/(eos[npoints-1][1] - eos[npoints-3][1]);
    eos[npoints-1][3] = (eos[npoints-1][2] - eos[npoints-2][2])/(eos[npoints-1][1] - eos[npoints-2][1]);
    eos[npoints-2][5] = eos[npoints-2][3]*eos[npoints-2][5];
    eos[npoints-1][5] = eos[npoints-1][3]*eos[npoints-1][5];

    // print in case needed
    if (print == true) {
        dm1.print(eos,npoints,7,true,"FSUEOS.txt");
    }
    return 0;
}

// interpolate EOS for use in high speed calculations
// npoints specifies the number of points between each original array point to get
int EOS :: splineeos(double** eosarr, int nrows, int npoints, double** &neweos, bool print) {
    double** spline1; double** spline2; double** spline3; double** spline4; double** spline5; double** spline6;
    double a,b,c,d,y1,x,y2,y3,y4,y5,y6;
    int n,k;
    
    // initialize the cubic spline arrays
    dm1.cubicspline(eosarr,spline1,nrows,0,1); dm1.cubicspline(eosarr,spline2,nrows,0,2);
    dm1.cubicspline(eosarr,spline3,nrows,0,3); dm1.cubicspline(eosarr,spline4,nrows,0,4);
    dm1.cubicspline(eosarr,spline5,nrows,0,5); dm1.cubicspline(eosarr,spline6,nrows,0,6);
    
    // calculate the new array number of rows and create the new array
    n = (nrows-1)*npoints + nrows;
    neweos = new double*[n];
    for (int i=0; i<n; i++) {
        neweos[i] = new double[7];
    }
    
    k = 0;  // row index starts at zero
    for (int i=0; i<(nrows-1); ++i) {
        // fill the array at every nth point with the known value
        neweos[k][0] = spline1[i][0]; neweos[k][1] = spline1[i][1];
        neweos[k][2] = spline2[i][1]; neweos[k][3] = spline3[i][1];
        neweos[k][4] = spline4[i][1]; neweos[k][5] = spline5[i][1];
        neweos[k][6] = spline6[i][1];
        k = k + 1;  // iterate

        // compute the splines
        for (int j=1; j<(npoints+1); ++j) {
            x = spline1[i][0] + (spline1[i+1][0] - spline1[i][0])*j/(npoints+1);
            a = (spline1[i+1][0] - x)/(spline1[i+1][0] - spline1[i][0]);
            b = 1.0-a;
            c = 1.0/6.0*(pow(a,3.0)-a)*pow(spline1[i+1][0]-spline1[i][0],2.0);
            d = 1.0/6.0*(pow(b,3.0)-b)*pow(spline1[i+1][0]-spline1[i][0],2.0);
            y1 = a*spline1[i][1] + b*spline1[i+1][1] + c*spline1[i][2] + d*spline1[i+1][2];
            
            a = (spline2[i+1][0] - x)/(spline2[i+1][0] - spline2[i][0]);
            b = 1.0-a;
            c = 1.0/6.0*(pow(a,3.0)-a)*pow(spline2[i+1][0]-spline2[i][0],2.0);
            d = 1.0/6.0*(pow(b,3.0)-b)*pow(spline2[i+1][0]-spline2[i][0],2.0);
            y2 = a*spline2[i][1] + b*spline2[i+1][1] + c*spline2[i][2] + d*spline2[i+1][2];
            
            a = (spline3[i+1][0] - x)/(spline3[i+1][0] - spline3[i][0]);
            b = 1.0-a;
            c = 1.0/6.0*(pow(a,3.0)-a)*pow(spline3[i+1][0]-spline3[i][0],2.0);
            d = 1.0/6.0*(pow(b,3.0)-b)*pow(spline3[i+1][0]-spline3[i][0],2.0);
            y3 = a*spline3[i][1] + b*spline3[i+1][1] + c*spline3[i][2] + d*spline3[i+1][2];

            a = (spline4[i+1][0] - x)/(spline4[i+1][0] - spline4[i][0]);
            b = 1.0-a;
            c = 1.0/6.0*(pow(a,3.0)-a)*pow(spline4[i+1][0]-spline4[i][0],2.0);
            d = 1.0/6.0*(pow(b,3.0)-b)*pow(spline4[i+1][0]-spline4[i][0],2.0);
            y4 = a*spline4[i][1] + b*spline4[i+1][1] + c*spline4[i][2] + d*spline4[i+1][2];

            a = (spline5[i+1][0] - x)/(spline5[i+1][0] - spline5[i][0]);
            b = 1.0-a;
            c = 1.0/6.0*(pow(a,3.0)-a)*pow(spline5[i+1][0]-spline5[i][0],2.0);
            d = 1.0/6.0*(pow(b,3.0)-b)*pow(spline5[i+1][0]-spline5[i][0],2.0);
            y5 = a*spline5[i][1] + b*spline5[i+1][1] + c*spline5[i][2] + d*spline5[i+1][2];

            a = (spline6[i+1][0] - x)/(spline6[i+1][0] - spline6[i][0]);
            b = 1.0-a;
            c = 1.0/6.0*(pow(a,3.0)-a)*pow(spline6[i+1][0]-spline6[i][0],2.0);
            d = 1.0/6.0*(pow(b,3.0)-b)*pow(spline6[i+1][0]-spline6[i][0],2.0);
            y6 = a*spline6[i][1] + b*spline6[i+1][1] + c*spline6[i][2] + d*spline6[i+1][2];
            
            neweos[k][0] = x; neweos[k][1] = y1; neweos[k][2] = y2;
            neweos[k][3] = y3; neweos[k][4] = y4; neweos[k][5] = y5; neweos[k][6] = y6;
            k = k + 1;  // iterate over the npoints to interpolate
        }
    }

    // fill last point
    neweos[k][0] = spline1[nrows-1][0]; neweos[k][1] = spline1[nrows-1][1];
    neweos[k][2] = spline2[nrows-1][1]; neweos[k][3] = spline3[nrows-1][1];
    neweos[k][4] = spline4[nrows-1][1]; neweos[k][5] = spline5[nrows-1][1];
    neweos[k][6] = spline6[nrows-1][1];

    // can print if needed
    if (print == true) {
        dm1.print(neweos, n, 7, print, "splineresult.txt");
    }

    // cleanup the arrays
    dm1.cleanup(spline1,nrows); dm1.cleanup(spline2,nrows); dm1.cleanup(spline3,nrows);
    dm1.cleanup(spline4,nrows); dm1.cleanup(spline5,nrows); dm1.cleanup(spline6,nrows);
    
    return n;   // return the number of new rows
}

// NEEDS TO BE UPDATED SINCE EOS DOESNT OUTPUT BA AND MSTAR AND THE EOS ALREADY COMES IN ARRAY FORM
// get bulk properties from the couplings (needs the EOS for SNM)
void EOS :: get_bulkproperties(string eosSNM, int nbcol, int encol, int prcol, int BAcol, int mstrcol, double gsoms2, double gwomw2, double gpomp2, double b, double c, double h, double lambda) {
    double **arrSNM;
    double BA, p0, mstr, K, gss, J, kf, L, gww;
    int nrows, ncols;

    // import SNM EOS
    dm1.importdata(eosSNM, arrSNM);
    nrows = dm1.rowcount(eosSNM);
    ncols = dm1.colcount(eosSNM);

    // get saturation properties from array
    BA = dm1.findmin(arrSNM,BAcol,nrows,ncols);
    p0 = arrSNM[dm1.findvalue(arrSNM,nrows,ncols,BA,BAcol,0.01)][nbcol];
    p0 = p0*pow(197.32698,3);
    kf = pow(3.0/2.0*pow(pi,2)*p0,1.0/3.0)/197.32698;
    mstr = arrSNM[dm1.findvalue(arrSNM,nrows,ncols,BA,BAcol,0.01)][mstrcol];
    cout << "BA: " << BA << endl;
    cout << "p0: " << p0/pow(197.32698,3) << endl;
    cout << "kf: " << kf << endl;
    cout << "mstr/m: " << mstr << endl;

    // get properties from analytic expressions
    gss = mNuc - mstr*mNuc;
    gww = mNuc + BA - sqrt(pow(kf,2) + pow(mstr*mNuc,2));
    K = get_K(kf,gss,gwomw2,gsoms2,b,c,h,gww);
    J = get_J(kf,gss,gpomp2,gww,lambda);
    p0 = p0/pow(197.32698,3);
    L = get_L(kf,gss,gsoms2,gww,gwomw2,gpomp2,h,lambda,b,c);
    cout << "K: " << K << endl;
    cout << "J: " << J << endl;
    cout << "L: " << L << endl;

    dm1.cleanup(arrSNM,nrows);  // clear the SNM array
}

//########################################################################################################
//########################################################################################################
// Beta Equilibrium

// newtons method for solving for the sigma field
// runs newtons method until unviable and then bisects
double EOS :: FSUgssnewton(double dens, double gsoms2, double b, double c, double t, double eps) {
    double dy, newt, ymin,error,y,min,max,sol,k;
    int ib, in, MAXIT;
    
    error = eps;           // set min error
    y = error*2;           // initialize the bisection
    min = 1e-5; max = mNuc+1.0; // set bounds on solution
    sol = (min+max)/2.0;    // solution starts as midpoint
    k = pow(3.0/2.0*pow(pi,2)*dens,1.0/3.0);    // convert density to fermi momentum
    
    ib=0;   // counter for bisection
    in = 0; // counter for newtons
    MAXIT = 100;    // maximum number of iterations for newtons allowed
    bool flag = false;

    while (abs(y)>error) {
        y = gss_FE(k,sol,gsoms2,b,c,t);
        dy = gssFE_derivative(k,sol,gsoms2,b,c,t);
        newt = sol - y/dy;

        // if newtons method yields a value outside the solution bounds or its reached its maximum iterations allowed then use bisection
        if (newt > max || newt<min || flag == true) {
            ymin = gss_FE(k,min,gsoms2,b,c,t);
            if (y*ymin<0) {
                max = sol;
            } else {
                min = sol;
            }
            sol = (min+max)/2.0;    // new midpoint
            ib = ib+1;              // bisection count
        } else {
            sol = newt;             // newtons solution
            in = in+1;              // newtons count
        }
        
        // flag newton as unviable if it has had too many counts
        if (in > MAXIT) {
            flag = true;
        }
    }
    return sol;
}

// charge neutrality equation for use in bisection
double EOS :: chneutralityeq(double kf, double gsoms2, double gwomw2, double gpomp2, double b, double c, double lambda, double h, double t) {
    double kp, kn, gss, gww, gpp, munmmup, y, mstar, np,eps,dens,ne,nmu;
    
    eps = 1e-8; // set error for newtons method
    dens = 2.0*pow(kf,3)/(3.0*pow(pi,2));   // get density from fermi momentum
    
    kp = kf*pow(1.0-t,1.0/3.0); kn = kf*pow(1.0+t,1.0/3.0); // get nucleon fermi momenta
    np = pow(kp,3.0)/(3.0*pow(pi,2));   // proton density

    // get meson fields
    gss = FSUgssnewton(dens,gsoms2,b,c,t,eps);
    mstar = mNuc-gss;
    gww = FSUgwwnewton(gwomw2,gpomp2,lambda,h,dens,t);         // try newton
    gpp = -0.5*gpomp2*dens*t/(1.0+lambda*gpomp2*pow(gww,2.0));

    munmmup = sqrt(pow(kn,2) + pow(mstar,2)) - sqrt(pow(kp,2) + pow(mstar,2)) - gpp;    // mun - mup = mue
    ne = qk.nb(munmmup,mE,2.0);
    nmu = qk.nb(munmmup,mMU,2.0);
    y = np - ne - nmu;
    return y;
    
}

// perform newtons method and bisection to get t for a given fermi momentum kf
double EOS :: get_t_betaeq(double kf, double gsoms2, double gwomw2, double gpomp2, double b, double c, double lambda, double h) {
    double error,y,min,max,eE,eMU,at, bt, kp, kn, rn, rp, mstar, dy, newt, gss, gww, gpp, dgssdt, dgppdt, munmmup,dens,t;
    int MAXIT, in, ib;
    bool flag = false;  // flag for marking when newtons method fails

    error = 1e-7;           // set min error
    y = error*2;            // initialize bisection
    min = 0.001; max = 1.0; // set bounds on solution
    MAXIT = 1000;           // max iterations for newtons method
    dens = 2.0*pow(kf,3)/(3.0*pow(pi,2));   // get density from kf
    t = (min+max)/2.0;      // set solution as midpoint
    in =0;  // counter for newtons
    ib =0;  // counter for bisection

    while (abs(y)>error) {
        kp = kf*pow(1.0-t,1.0/3.0); kn = kf*pow(1.0+t,1.0/3.0); // nucleon fermi momenta
        gss = FSUgssnewton(dens,gsoms2,b,c,t,1e-8); // get sigma field
        mstar = mNuc-gss;   // effective mass
        rp = sqrt(pow(kp,2.0)+pow(mstar,2.0)); rn = sqrt(pow(kn,2.0)+pow(mstar,2.0));   // convenient quantity to define
        gww = FSUgwwnewton(gwomw2,gpomp2,lambda,h,dens,t);  // get omega field
        gpp = -0.5*gpomp2*dens*t/(1.0+lambda*gpomp2*pow(gww,2.0));  // get rho field
        munmmup = sqrt(pow(kn,2) + pow(mstar,2)) - sqrt(pow(kp,2) + pow(mstar,2)) - gpp;    // mue = mun-mup

        // convenient quantities to define
        at = kp*rp*0.5 - pow(mstar,2.0)*0.5*log(kp+rp) + pow(mstar,2.0)*log(mstar) + 0.5*kn*rn - pow(mstar,2.0)*0.5*log(kn+rn);
        bt = -kp*mstar*0.5/rp - kn*mstar*0.5/rn + mstar*log(kp+rp) + mstar*log(kn+rn) + pow(mstar,3.0)*0.5/(rp*kp+rp*rp) + pow(mstar,3.0)*0.5/(rn*kn+rn*rn) - mstar*(1.0-2.0*log(mstar));
        
        // derivatives (g_sigma dsigma/dt   and    g_rho drho/dt)
        dgssdt = -gsoms2*mstar/(3.0*pow(pi,2.0))*pow(kf,3.0)*((rp-rn)/(rp*rn)) /(-1.0/pow(pi,2.0)*gsoms2*at +gsoms2*mstar/pow(pi,2.0)*bt - 1.0 - gsoms2*(2.0*b*mNuc*gss + 3.0*c*gss*gss) );
        dgppdt = -0.5*gpomp2*dens*(1.0+gwomw2*lambda*pow(gpp,2.0)+3.0*gwomw2*h*pow(gww,2.0))/(lambda*pow(gpp,2.0)*gwomw2*(1.0-3.0*gpomp2*lambda*gww*gww) + (1.0*3.0*gww*gww*h*gwomw2)*(1.0+gpomp2*lambda*gww*gww) );
        
        // fermi momenta zero if sqrt is negative
        eE = sqrt(pow(munmmup,2.0)-pow(mE,2.0));
        eMU = sqrt(pow(munmmup,2.0)-pow(mMU,2.0));
        if (munmmup < mE) {
            eE = 0;
        }
        if (munmmup < mMU) {
            eMU = 0;
        }
        
        // Netwons method with bisection 
        dy = -pow(kf,3.0) - 3.0*munmmup*(eE + eMU)*(pow(kf,3.0)/(3.0*kn*rn) + pow(kf,3.0)/(3.0*kp*rp) + dgssdt*mstar*((rn-rp)/(rp*rn)) - dgppdt);
        y = chneutralityeq(kf, gsoms2, gwomw2, gpomp2, b, c, lambda, h, t);
        dy = dy/(3.0*pow(pi,2.0));
        newt = t - y/dy;

        // if newtons method yields a value outside the solution bounds or its reached its maximum iterations allowed then use bisection
        if (newt >= max || newt<min || flag == true) {
            if (y*chneutralityeq(kf, gsoms2, gwomw2, gpomp2, b, c, lambda, h, min)<0) {
                max = t;
            } else {
                min = t;
            }
            t = (min+max)/2.0;  // new midpoint
            ib = ib+1;  // bisection count
        } else {
            t = newt;
            y = chneutralityeq(kf, gsoms2, gwomw2, gpomp2, b, c, lambda, h, t);
            in = in+1;  // newton count
        }
        
        // flag newtons method if it takes too long
        if (in>MAXIT) {
            flag = true;
        }

        // give bisection more iterations and exit if it fails
        if (ib>(MAXIT+100)) {
            cout << "bisection beta eq failed for k: " << kf << endl;
            exit(0);
        }
    }
    return t;
}

// newton and bisection method to get omega field
double EOS :: FSUgwwnewton(double gwomw2, double gpomp2, double lambda, double h, double dens, double t) {
    double dy, gpp, newt, ymin, dgpp, error,y,min,max,x;
    int MAXIT, in,ib;
    
    error = 1e-8;           // set min error
    y = error*2;            // initialize bisection
    min = 1e-5; max = 4000; // set bounds on solution
    MAXIT = 500;   // maximum number of iterations for newtons method
    in = 0; // count for newtons
    ib = 0; // count for bisection
    bool flag = false;  // flag for newtons method failing
    x = (min+max)/2.0;   // solution starts at midpoint
    
    while (abs(y)>error) {
        // newtons method
        gpp = -0.5*gpomp2*dens*t/(1.0+lambda*gpomp2*pow(x,2.0));
        dgpp = 0.5*gpomp2*dens*t*pow(1.0+lambda*gpomp2*pow(x,2.0),-2.0)*lambda*gpomp2*2.0*x;
        y = x + lambda*pow(gpp,2.0)*gwomw2*x + h*gwomw2*pow(x,3.0) - gwomw2*dens;
        dy = 1.0 + lambda*gwomw2*pow(gpp,2.0) + lambda*2.0*gpp*dgpp*gwomw2*x + 3.0*h*gwomw2*pow(x,2.0);
        newt = x - y/dy;
        
        // if newtons method yields a value outside the solution bounds or its reached its maximum iterations allowed then use bisection
        if (newt > max || newt<min || flag == true) {
            gpp = -0.5*gpomp2*dens*t/(1.0+lambda*gpomp2*pow(min,2.0));
            ymin = min + lambda*pow(gpp,2.0)*gwomw2*min + h*gwomw2*pow(min,3.0) - gwomw2*dens;
            if (y*ymin<0) {
                max = x;
            } else {
                min = x;
            }
            x = (min+max)/2.0;  // new midpoint
            ib = ib+1;  // bisection count
        } else {
            x = newt;   // newtons solution
            in = in+1;  // newton count
        }
        
        // flag newton as unviable and try bisection
        if (in > MAXIT) {
            flag = true;
        }
        // if bisection also fails then exit
        if (ib > MAXIT) {
            cout << "bisection failed: " << endl;
            cout << "bisect: " << ib << "  " << y << "  " << x << "  " << gpp << "  " << dens << "  " << t << endl;
            exit(0);
        }
    }
    return x;
}

// Use the thermodynamic stability method to find the crust core transition point and glue the two together
// has optional ability to print the new EOS
// crust EOS is of the form (nb, en, p, dpde, mub)
int EOS :: ThermalCrust(double** crust, double** core, double** &neweos, int nrowscore, int nrowscrust, bool print, int nbcol, int prcol, int Keffcol) {
    int nrows;
    int ncols_core = 7;
    double trdens = 0;
    double x = core[nrowscore-1][Keffcol];
    double xp;

    // find where the compressibility becomes negative
    for (int i=(nrowscore-1);i>0;i=i-1) {
        xp = core[i][Keffcol];  // store next effK
        if (x>0) {          // check if effK starts stable
            if (x*xp<0) {       // check if effK becomes stable
                trdens = core[i][0];    // mark tr dens
                break;
            }
        }
        x = xp;     // store old effK
    }

    // Get the row number for the crust and core transitions
    int crustrow = dm1.findvalue(crust, nrowscrust,5,trdens,nbcol,1.0);
    int corerow = dm1.findvalue(core,nrowscore,ncols_core,trdens,nbcol,1.0);
    //cout << "Transition at: " << trdens << " 1/fm3" << endl; 
    while (core[corerow][nbcol] < crust[crustrow][nbcol]) {     // make sure the density is monotonic
        corerow = corerow + 1;
        if (corerow > (nrowscore-1)) {
            cout << "Thermal crust invalid nb for core" << endl;
            exit(0);
        }
    }

    while (core[corerow][prcol] < crust[crustrow][prcol]) {     // make sure the pressure is monotonic
        corerow = corerow + 1;
        if (corerow > (nrowscore-1)) {
            cout << "Thermal crust invalid nb for core" << endl;
            exit(0);
        }
    }

    // stitch the EOS
    nrows = crustrow + 1 + (nrowscore-corerow-1);
    neweos = new double*[nrows];
    for (int i=0; i<nrows; i++) {
        neweos[i] = new double[5];
    }

    for (int i=0; i<(crustrow+1); ++i) {
        neweos[i][0] = crust[i][0];
        neweos[i][1] = crust[i][1];
        neweos[i][2] = crust[i][2];
        neweos[i][3] = crust[i][3];
        neweos[i][4] = crust[i][4];
    }
        
    for (int i=(corerow+1); i<nrowscore; ++i) {
        neweos[crustrow+1+i-corerow-1][0] = core[i][0];
        neweos[crustrow+1+i-corerow-1][1] = core[i][1];
        neweos[crustrow+1+i-corerow-1][2] = core[i][2];
        neweos[crustrow+1+i-corerow-1][3] = core[i][3];
        neweos[crustrow+1+i-corerow-1][4] = core[i][4];
    }
    
    if (print == true) {
        dm1.print(neweos,nrows,5,true,"FSUEOSC.txt");
    }
    return nrows;
}

// get the effective compressibility for determining the crust core transition
double EOS :: effK2(double kf, double gss, double gsoms2, double gww, double gwomw2, double gpomp2, double gpp, double h, double lambda, double b, double c, double t) {
    double res, mstar, kp,kn,rp,rn, dkndpn,dkpdpp,dens,Yp,termNn,termPp,termPn,integral,gsdsdkn,gsdsdkp,dmunpn,dmuppp,dmuppn;
    mstar = mNuc - gss; // get the effective mass
    kp = kf*pow(1.0-t,1.0/3.0); kn = kf*pow(1.0+t,1.0/3.0); // get the fermi momenta
    rp = sqrt(pow(kp,2) + pow(mstar,2)); rn = sqrt(pow(kn,2) + pow(mstar,2));  // convenient quantities to define
    
    // stops errors for 1/0
    if (kp == 0) {
        kp = 1e-15;
    }
    dkndpn = pow(pi/kn,2.0);    // dk_n/drho_n
    dkpdpp = pow(pi/kp,2.0);    // dk_p/drho_p
    dens = 1.0/(3.0*pow(pi,2.0))*pow(kn,3.0) + 1.0/(3.0*pow(pi,2.0))*pow(kp,3.0);   // total baryon density
    Yp = 0.5*(1.0-t);   // proton fraction

    // algebraic terms
    termNn = 0.25*(4.0*gwomw2 + gpomp2 + lambda*gwomw2*gpomp2*pow(gpp,2.0) + 8.0*gwomw2*gpomp2*lambda*gpp*gww + (3.0*h+4.0*lambda)*gpomp2*gwomw2*pow(gww,2.0))/
            (lambda*gwomw2*pow(gpp,2.0)*(1.0 - 3.0*gpomp2*lambda*pow(gww,2.0)) + (1.0 + 3.0*h*gwomw2*pow(gww,2.0))*(1.0+gpomp2*lambda*pow(gww,2.0)));
    termPp = 0.25*(4.0*gwomw2 + gpomp2 + lambda*gwomw2*gpomp2*pow(gpp,2.0) - 8.0*gwomw2*gpomp2*lambda*gpp*gww + (3.0*h+4.0*lambda)*gpomp2*gwomw2*pow(gww,2.0))/
            (lambda*gwomw2*pow(gpp,2.0)*(1.0 - 3.0*gpomp2*lambda*pow(gww,2.0)) + (1.0 + 3.0*h*gwomw2*pow(gww,2.0))*(1.0+gpomp2*lambda*pow(gww,2.0)));
    termPn = 0.25*(4.0*gwomw2 - gpomp2 - lambda*gwomw2*gpomp2*pow(gpp,2.0) + (-3.0*h+4.0*lambda)*gpomp2*gwomw2*pow(gww,2.0))/
            (lambda*gwomw2*pow(gpp,2.0)*(1.0 - 3.0*gpomp2*lambda*pow(gww,2.0)) + (1.0 + 3.0*h*gwomw2*pow(gww,2.0))*(1.0+gpomp2*lambda*pow(gww,2.0)));
    integral = gsoms2/(2.0*pow(pi,2.0))*(pow(mstar,2.0)*(2.0 + kn/rn + kp/rp + 6.0*log(mstar) - 3.0*log(kn+rn) - 3.0*log(kp+rp))
                        - pow(mstar,4.0)*(1.0/(kp*rp+rp*rp) + 1.0/(kn*rn+rn*rn)) + kn*rn + kp*rp);
    
    // get derivatives of sigma field
    gsdsdkn = gsoms2*pow(kn,2.0)*mstar/(rn*pow(pi,2.0))*pow(1.0 + 3.0*c*gsoms2*pow(gss,2.0) + 2.0*b*mNuc*gsoms2*gss + integral,-1.0);   // g_sigma dsigma/dk_n
    gsdsdkp = gsoms2*pow(kp,2.0)*mstar/(rp*pow(pi,2.0))*pow(1.0 + 3.0*c*gsoms2*pow(gss,2.0) + 2.0*b*mNuc*gsoms2*gss + integral,-1.0);    // g_sigma dsigma_dk_p
    
    // derivatives of chemical potentials
    dmunpn = dkndpn*(kn-mstar*gsdsdkn)/rn + termNn; // dmu_n/drho_n
    dmuppp = dkpdpp*(kp-mstar*gsdsdkp)/rp + termPp; // dmu_p/drho_p
    dmuppn = dkndpn*(-mstar*gsdsdkn)/rp + termPn;   // dmu_p/drho_n

    res = dens/4.0*( (dmunpn + 2.0*dmuppn + dmuppp) + 2.0*(1.0-2.0*Yp)*(dmunpn - dmuppp) + pow(1.0-2.0*Yp,2.0)*(dmunpn - 2.0*dmuppn + dmuppp)
        - pow(dmunpn - dmuppp + (1.0-2.0*Yp)*(dmunpn - 2.0*dmuppn + dmuppp),2.0)/(dmunpn - 2.0*dmuppn + dmuppp) );
    cout << kf << "  " << res << endl;
    return res;
}

// get EOS for a given t (t=1 for PNM and t=0 for SNM)
// output is (nb fm-3, en MeV/fm3, pr MeV/fm3, B/A MeV, mstar MeV) 
int EOS :: get_PNM_SNM_EOS(double params[7], double** &eos, int npoints, bool print, double t) {
    double k, en, pr, dens, gww, gss, gpp, mstarn, kn, mun, check, conv_mev4, p0f;
    double ssize, gsoms2,gwomw2,gpomp2,b,c,h,lambda;

    conv_mev4 = conv1.energyCONV(6,0); // mev4 to mev/fm3 conversion
    k = 40.0;   // initial fermi momentum  
    p0f = 20.0; // N TIMES SATURATION DENSITY
    ssize = (-10.4687 + 55.5001*pow(p0f,7.0/25))/npoints;   // optimal setp size 
    
    // import the couplings
    gsoms2 = params[0];
    gwomw2 = params[1];
    gpomp2 = params[2];
    b = params[3];
    c = params[4];
    h = params[5];
    lambda = params[6];

    // initialize array for the EOS
    eos = new double*[npoints];
    for (int i=0; i<npoints; i++) {
        eos[i] = new double[5];
    }
    
    for (int i=0; i<npoints; ++i) {
        dens = 2.0*pow(k,3)/(3.0*pow(pi,2));    // density from kf
        kn = k*pow(1.0+t,1.0/3.0);  // neutron fermi momentum

        // get meson fields
        gss = FSUgssnewton(dens,gsoms2,b,c,t,1e-7);
        gww = FSUgwwnewton(gwomw2,gpomp2,lambda,h,dens,t);
        gpp = -0.5*gpomp2*dens*t/(1.0+lambda*gpomp2*pow(gww,2.0));

        mstarn = mN - gss;  // effective mass

        // get EOS
        mun = sqrt(pow(kn,2) + pow(mstarn,2)) + gww - 0.5*gpp;
        en = get_en(k,t,gss,gsoms2,gww,gwomw2,gpp,gpomp2,b,c,h,lambda);
        pr = get_pr(k,t,gss,gsoms2,gww,gwomw2,gpp,gpomp2,b,c,h,lambda);
        eos[i][0] = dens*conv_mev4; // density
        eos[i][1] = en*conv_mev4;   // energy density
        eos[i][2] = pr*conv_mev4;   // pressure
        eos[i][3] = en/dens - mNuc; // binding energy
        eos[i][4] = mstarn;         // effective mass

        check = mun*dens - en - pr; // check for thermodynamic consistence
        if (abs(check)>1.0) {
            cout << "consistency fail: " << check << endl;
            exit(0);
        }

        k = k + ssize*log(1.0+k);
    }

    // optional print
    if (print == true) {
        dm1.print(eos,npoints,5,true,"FSUEOSNM.txt");
    }
    return npoints;
}