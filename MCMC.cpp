#include "NumMethods.hpp"
#include "nucEOS.hpp"
#include "Conversions.hpp"
#include "MCMC.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <chrono>

using namespace std;
data2 dm2;
EOS nuc1;
nummeth nm1;
Convert conv2;
statistics stats;

const double pi = 4.0*atan(1.0);
const double mP = 939; //938.27231;    // Mass of proton (MeV)
const double mN = 939; //939.56542052;
const double mNuc = (mP+mN)/2.0;

//Box muller method
// returns a value from the normal distribution centered at some mean and stddev
double statistics :: rand_normal(double mean, double stddev) {
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

// returns a value from the uniform distribution centered at cntr with +- width w
// must be seeded
double statistics :: rand_uniform(double cntr, double w) {
    double x,y;
    x = rand()*1.0/RAND_MAX;
    y = cntr + 2.0*w*x - w;
    return y;
}

//######################################################################################################
//######################################################################################################

// target function for sampling a distribution given the covariance matrix
// takes in the proposed bulk properties and the old bulk properties to compute a chi square value
double RMFtarget(double params[7], double olds[7], double** cov, double** INIT, bool flag) {
    double chisqp = 0; double chisq0 = 0;   // set initial chi squares to zero
    double vec1[7]; double vec2[7]; double newparams[7];
    double BA,kf,J,mstar,K,L,h,p0,res;

    // initialize the chi square vectors
    for (int i=0; i<7; ++i) {
        vec1[i] = 0;
        vec2[i] = 0;
    }

    // (x-mu)^T COV^-1 (x-mu)
    for (int i=0; i<7; ++i) {
        for (int j=0; j<7; ++j) {
            vec1[i] = vec1[i] + cov[i][j]*(olds[j]-INIT[j][2]);
            vec2[i] = vec2[i] + cov[i][j]*(params[j]-INIT[j][2]);
        }
        chisq0 = chisq0 + vec1[i]*(olds[i]-INIT[i][2]);
        chisqp = chisqp + vec2[i]*(params[i]-INIT[i][2]);
    }

    // rescale the bulk properties to normal values;
    BA = params[0]*INIT[0][3]; kf = params[1]*INIT[1][3]; mstar = params[2]*INIT[2][3]*mNuc; 
    K = params[3]*INIT[3][3]; J = params[4]*INIT[4][3]; L = params[5]*INIT[5][3]; h = params[6]*INIT[6][3];

    // get the couplings to see if the EOS is unrealistic. If it is give it a negative probability so it doesnt get represented in the distribution
    // only happens for extreme EOS
    p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
    nuc1.get_parameters(BA,p0,J,mstar,K,L,h,newparams,false,flag);
    if (newparams[6] < 0) {
        return -1.0;
    }
    
    res = exp(chisq0/2.0 - chisqp/2.0);
    return res;
}

// Get a sample of a distribution given the covariance matrix and an initialization vector
// INIT vector contains suggested initial values and std deviations as well as means and scaling factors
int MCMC :: RMF(int npoints, double** INIT, string covdata) {
    srand(time(0)); // random seed
    double olds[7]; double cands[7]; int counts[7]; double stds[7]; double arate[7]; double agoal[7]; double params[7];
    double BA,kf,J,mstar,K,L,h,p0,r,a;
    ofstream out("prior.txt");  // output of the results

    // import cov matrix
    // specify acceptance rate goals, starting point for bulk properties and their widths, acceptance counter, averages
    double** cov; // cov matrix
    dm2.importdata(covdata,cov);
    for (int i=0; i<7; ++i) {
        agoal[i] = 0.5;
        counts[i] = 0;
        olds[i] = INIT[i][0];
        stds[i] = INIT[i][1];
    }
    
    // get the couplings to see if the EOS is unrealistic. Should never get flagged unless the first guess is not good since RMFtarget() already checks
    // just extra safety measure
    BA = olds[0]*INIT[0][3]; kf = olds[1]*INIT[1][3]; mstar = olds[2]*INIT[2][3]*mNuc; 
    K = olds[3]*INIT[3][3]; J = olds[4]*INIT[4][3]; L = olds[5]*INIT[5][3]; h = olds[6]*INIT[6][3];
    p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
    nuc1.get_parameters(BA,p0,J,mstar,K,L,h,params,false,true);
    if (params[6] < 0) {
        cout << "problem" << endl;
        exit(0);
    } 

    // BURN IN PHASE
    for (int i=0; i<1000000; ++i) {
        for (int j=0; j<7; ++j) {
            for (int k=0; k<7; ++k) {
                cands[k] = olds[k]; 
            }

            cands[j] = stats.rand_uniform(olds[j],stds[j]); // generate candidate point centered at old point within +- stdav
            a = RMFtarget(cands,olds,cov,INIT,true);    // get probability of acceptance from the likelihood
            if (a>1.0) {    // probability shouldnt be greater than 1
                a = 1.0;        
            }

            r = 1.0*rand()/RAND_MAX;    // generate random number from 0 to 1
            if (r <= a) {   // accept candidate with probability a 
                olds[j] = cands[j];
                counts[j] = counts[j] + 1;      // add 1 to the accpetance count
            }

            //ratemonitor every 100 iterations
            if ((i+1)%100 == 0) {
                arate[j] = 1.0*counts[j]/100.0; // find out how many points were accepted in 100 run period
                counts[j] = 0;  // reset the counter
                
                if (arate[j] < agoal[j]) {
                    stds[j] = 0.9*stds[j];                 // if acceptance rate is too low then decrease the range
                } else if (arate[j] > agoal[j]) {
                    stds[j] = 1.1*stds[j];                 // if acceptance rate is too high then increase the range
                }
                
                cout << olds[j] << ": " << stds[j] << " | ";    // print the burn in values
            }
        }
    
        if ((i+1)%100 == 0) {
            cout << endl;
        }
    }

    // ACTUAL RUN 
    for (int i=0; i<npoints; ++i) {
        for (int j=0; j<7; ++j) {
            for (int k=0; k<7; ++k) {
                cands[k] = olds[k];
            }
            cands[j] = stats.rand_uniform(olds[j],stds[j]); // generate candidate point centered at old point within +- stdav
            a = RMFtarget(cands,olds,cov,INIT,true);    // get probability of acceptance from the likelihood
            if (a>1.0) {    // probability shouldnt be greater than 1
                a = 1.0;        
            }

            r = 1.0*rand()/RAND_MAX;    // generate random number from 0 to 1
            if (r <= a) {   // accept candidate with probability a 
                olds[j] = cands[j];
                counts[j] = counts[j] + 1;      // add 1 to the accpetance count
            }
        }
        
        // print out the bulk properties from the sampling
        out << setprecision(10) << olds[0]*INIT[0][3] << "  " << olds[1]*INIT[1][3] << "  " << olds[2]*INIT[2][3] << "  " << olds[3]*INIT[3][3] << "  " << olds[4]*INIT[4][3] 
        << "  " << olds[5]*INIT[5][3] << "  " << olds[6]*INIT[6][3] << endl;
    }
    
    // print out the acceptance rates at the end (should be close the agoal)
    cout << "acc rates: " << 1.0*counts[0]/npoints << "  " << 1.0*counts[1]/npoints << "  " << 1.0*counts[2]/npoints << "  " << 1.0*counts[3]/npoints
    << "  " << 1.0*counts[4]/npoints << "  " << 1.0*counts[5]/npoints << "  " << 1.0*counts[6]/npoints << endl;
    
    // cleanup the covariance matrix
    dm2.cleanup(cov,7);
    return 0;
}

// Calculate astrophysical priors given a set of parameters
// nMR is number of mass radius points from NICER
// nL is number of Ligo mass and TD points
// output is of the form (BA, kf, mstar, K, J, L, h, NICER RADII, ..., NICER PRESSURES, ..., LIGO RADIUS, LIGO TD, MAXIMUM MASS, PRESSURE OF MMAX)
double MCMC :: Observables(string paramset, double** crust, int nrowscrust, double** NICER, int nMR, double** LIGO, int nL, string outname, int nruns) {
    double BA, p0, kf, J, mstar, K, L, h,cv;   // bulk parameters
    double params[7];    // arrays for the calculated parameters
    int npoints = 100;      // npoints for the eos
    int nrows1, n;      // nrows after adding crust and spline
    cv = conv2.energyCONV(0,1);   // MeV/fm3 to unitless
    bool unstable = false;  
    
    // Observables
    double ILQR[4];
    int ncols = 7;          // ncols on the EOS (nb, en, pr, dpde, mub, dPdrho, Keff)
    int encol = 1; int prcol = 2; int dpdecol = 3; int nbcol = 0; int Keffcol = 6;

    ofstream out(outname);     // output file
    double** bulks;                     // array to store bulk priors
    int nrows = dm2.rowcount(paramset);    // count nrows for bulk priors
    dm2.importdata(paramset,bulks);        // import bulk priors

    double lkl;
    // loop through all bulk priors and calculate astrophysical properties
    #pragma omp parallel num_threads(12)
    #pragma omp for ordered schedule(static,1) private(BA,kf,J,mstar,K,L,h,p0,params,ILQR,n,nrows1,lkl)
    for (int i=0; i<nruns; ++i) {
        BA = bulks[i][0];
        kf = bulks[i][1]; 
        J = bulks[i][4]; 
        mstar = bulks[i][2]*mNuc; 
        K = bulks[i][3]; 
        L = bulks[i][5]; 
        h = bulks[i][6];
        p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
        
        double** core1; double** splinecore1; double** eosc; double** PRM;

        // prRM array (last row is max mass config)
        PRM = new double*[nMR+nL+1];
        for (int j=0; j<(nMR+nL+1); j++) {
            PRM[j] = new double[3];
        }
        
        // FILL WITH MASS DATA FROM NICER AND LIGO
        for (int k=0; k<nMR; ++k) {
            PRM[k][2] = NICER[k][0];
        }
        for (int k=0; k<nL; ++k) {
            PRM[nMR+k][2] = LIGO[k][0];   
        }

        nuc1.get_parameters(BA,p0,J,mstar,K,L,h,params,false,true);     // solve for parameters given bulk properties
        nuc1.get_EOS(params,core1,npoints,false,unstable);     // get the EOS (nb,en,pr,dpde,mub,dPdrho,Keff)
        n = nuc1.splineeos(core1, npoints, 20, splinecore1,false);    // interpolate the EOS with 20 additional points per gap (for speed) so the eos now has ~ 20*100 points
        nrows1 = nuc1.ThermalCrust(crust,splinecore1,eosc,n,nrowscrust,false,nbcol,prcol,Keffcol);      // get the crust transition density and add crust
        nm1.pretovconv(eosc, encol, prcol, cv, cv, nrows1);     // convert the EOS to dimensionless units
        lkl = nm1.MCMCtovpoint(1e-4,eosc,nrows1,ncols-1,encol,prcol,dpdecol,100,8.0*cv,nMR,nL,NICER,LIGO,PRM);        // get radii for masses, and central pressures to calc ILQ
        nm1.RloveMCMC(eosc,dpdecol,encol,prcol,1e-5,ILQR,PRM[nMR][0],nrows1,ncols-1);                // get TD for 1.4 
        
        dm2.cleanup(core1,npoints); dm2.cleanup(splinecore1,n); dm2.cleanup(eosc,nrows1);   // cleanup arrays
        #pragma omp ordered
        cout << i << "  " << params[0] << "  " << params[1] << "  " << params[2] << "  " << params[3] << "  " << params[4] << "  " << params[5] << "  " << params[6] << endl;
        out << setprecision(10) << BA << "  " << kf << "  " << mstar/mNuc << "  " << K << "  " << J << "  " << L << "  " << h << "  ";  // output bulk priors
        
        // output observables >>  (NICER RADII, ..., NICER PRESSURES, ..., LIGO RADIUS, LIGO TD, MAXIMUM MASS, PRESSURE OF MMAX)
        for (int j=0; j<nMR; ++j) {
            out << conv2.rnonetokm(PRM[j][1]) << "  ";
        }
        for (int j=0; j<nMR; ++j) {
            out << PRM[j][0]/cv << "  ";
        }
        out << conv2.rnonetokm(ILQR[3]) << "  " << ILQR[1] << "  " << PRM[nMR+nL][2] << "  " << PRM[nMR+nL][0]/cv << endl;

        dm2.cleanup(PRM,nMR+nL+1);  // cleanup observables array
    }

    dm2.cleanup(bulks,nrows);
    return 0;
}

// target liklihood function for the calibration of models to neutron star observables and XEFT 
// takes in the proposed bulk properties and the old bulk properties to compute a chi square value
double targetAstro(double* &props, double* &olds, double** cov, double** DATA, double** crust, int nrowscrust, bool flag, double** NICER, int nMR, double** XEFTdata, int nrowsXEFT, double** LIGO, int nL) {
    double chisqp = 0; double chisq0 = 0;   // chisq of prior
    double vec1[7]; double vec2[7];         // temporary arrays for matrix mulitplication
    double BA, p0, kf, J, mstar, K, L, h, res;   // bulk parameters
    double R14TDtheory0, R14TDtheoryp;            // theoretical values
    double oldparams[7]; double newparams[7];       // arrays for the calculated parameters
    double** eos1; double** eos2; double** core1; double** core2; double** splinecore1; double** splinecore2; double** PRM1; double** PRM2;       // core eos, interpolated core, and total eos
    int npoints = 100;      // npoints for the eos
    int nrows1; int nrows2; int n = 0;      // nrows after adding crust and spline
    double cv = conv2.energyCONV(0,1);   // MeV/fm3 to unitless
    bool unstable = false;

    // PRM matricies that store observables
    PRM1 = new double*[nMR+nL+1];
    for (int j=0; j<(nMR+nL+1); j++) {
        PRM1[j] = new double[3];
    }

    PRM2 = new double*[nMR+nL+1];
    for (int j=0; j<(nMR+nL+1); j++) {
        PRM2[j] = new double[3];
    }

    // Fill masses with NICER DATA AND LIGO DATA
    for (int k=0; k<nMR; ++k) {
        PRM1[k][2] = NICER[k][0];   
        PRM2[k][2] = NICER[k][0];
    }
    for (int k=0; k<nL; ++k) {
        PRM1[nMR+k][2] = LIGO[k][0];   
        PRM2[nMR+k][2] = LIGO[k][0];
    }
    
    ofstream out("DEBUG.txt");
    int nbcol = 0; int encol = 1; int prcol = 2; int dpdecol = 3; int keffcol = 6; int ncols = 7;

    // Observables
    double R14TDactual = LIGO[0][1];
    double R14TDerrU = LIGO[0][2]; double R14TDerrL = LIGO[0][3];   // upper and lower errors since error is asymmetrical
    double ILQR1[4]; double ILQR2[4];   // vectors to store tidal deformability data
    double lkl0, lklp, max0, maxp;

    //Make sure proposed changes are physical
    BA = props[0]*DATA[0][3]; kf = props[1]*DATA[1][3]; mstar = props[2]*DATA[2][3]*mNuc; 
    K = props[3]*DATA[3][3]; J = props[4]*DATA[4][3]; L = props[5]*DATA[5][3]; h = props[6]*DATA[6][3];
    p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
    nuc1.get_parameters(BA,p0,J,mstar,K,L,h,newparams,false,flag);     // solve for parameters given bulk properties
    // Check if lambda is negative (eff rho mass is imaginary)
    if (newparams[6] < 0 || newparams[5] < 0) {
        return -1.0;
    }
    nuc1.get_EOS(newparams,core1,npoints,false,unstable);
    dm2.cleanup(core1,npoints);

    // check to make sure the EOS is stable
    if (unstable == true) {
        cout << "EOS is unstable for: " << L << endl;
        dm2.cleanup(core2,npoints);
        return -1.0;
    }

    // initialize the temp arrays and get the prior distribution X^2 = (x-mu)^T Cov (x-mu)
    vec1[0] = 0; vec1[1] = 0; vec1[2] = 0; vec1[3] = 0; vec1[4] = 0; vec1[5] = 0; vec1[6] = 0;
    vec2[0] = 0; vec2[1] = 0; vec2[2] = 0; vec2[3] = 0; vec2[4] = 0; vec2[5] = 0; vec2[6] = 0;
    for (int i=0; i<7; ++i) {
        for (int j=0; j<7; ++j) {
            vec1[i] = vec1[i] + cov[i][j]*(olds[j]-DATA[j][2]);
            vec2[i] = vec2[i] + cov[i][j]*(props[j]-DATA[j][2]);
        }
        chisq0 = chisq0 + vec1[i]*(olds[i]-DATA[i][2]);
        chisqp = chisqp + vec2[i]*(props[i]-DATA[i][2]);
    }
    double prior = exp(chisq0/2.0 - chisqp/2.0);
    
    // output for DEBUG
    out << setprecision(10) << olds[0]*DATA[0][3] << "  " << olds[1]*DATA[1][3] << "  " << olds[4]*DATA[4][3] << "  " << olds[2]*DATA[2][3]*mNuc << "  " << olds[3]*DATA[3][3] << "  " << olds[5]*DATA[5][3] << "  " << olds[6]*DATA[6][3] << endl;
    out << setprecision(10) << props[0]*DATA[0][3] << "  " << props[1]*DATA[1][3] << "  " << props[4]*DATA[4][3] << "  " << props[2]*DATA[2][3]*mNuc << "  " << props[3]*DATA[3][3] << "  " << props[5]*DATA[5][3] << "  " << props[6]*DATA[6][3] << endl;
    
    // BEGIN PARALLEL
    #pragma omp parallel sections private(BA,kf,J,mstar,K,L,h,p0,n) 
    {
    #pragma omp section
    {
    BA = olds[0]*DATA[0][3]; kf = olds[1]*DATA[1][3]; mstar = olds[2]*DATA[2][3]*mNuc;      // get bulk properties by rescaling 
    K = olds[3]*DATA[3][3]; J = olds[4]*DATA[4][3]; L = olds[5]*DATA[5][3]; h = olds[6]*DATA[6][3];
    p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0); // get density
    cout << BA << "  " << kf << "  " << mstar << "  " << K << "  " << J << "  " << L << "  " << h << endl;
    nuc1.get_parameters(BA,p0,J,mstar,K,L,h,oldparams,false,flag);     // solve for parameters given bulk properties
    nuc1.get_EOS(oldparams,core1,npoints,false,unstable);     // get the EOS (nb,en,pr,dpde,mub,dPdrho,effK)
    n = nuc1.splineeos(core1, npoints, 20, splinecore1,false);    // interpolate the EOS with 10 additional points per gap (for speed)
    nrows1 = nuc1.ThermalCrust(crust,splinecore1,eos1,n,nrowscrust,false,nbcol,prcol,keffcol);      // get the crust transition density adn add crust
    nm1.pretovconv(eos1, encol, prcol, cv, cv, nrows1);     // convert the EOS to dimensionless units
    lkl0 = nm1.MCMCtovpoint(1e-4,eos1,nrows1,ncols,encol,prcol,dpdecol,100,5.0*cv,nMR,nL,NICER,LIGO,PRM1);        // get 1.4 solar mass star radius
    nm1.RloveMCMC(eos1,dpdecol,encol,prcol,1e-5,ILQR1,PRM1[nMR][0],nrows1,ncols-1);                // get TD for 1.4
    max0 = PRM1[nMR+nL][2];
    R14TDtheory0 = ILQR1[1];

    // store Astro Data
    for (int k=0; k<nMR; ++k) {
        olds[7+k] = conv2.rnonetokm(PRM1[k][1]);     // store radii
        olds[7+nMR+k] = PRM1[k][0]; // store icp
    }
    olds[7+nMR*2] = conv2.rnonetokm(ILQR1[3]);         // store the 1.4 radius
    olds[7+nMR*2+1] = R14TDtheory0;             // store 1.4 TD
    olds[7+nMR*2+2] = max0;      // Mmax
    olds[7+nMR*2+3] = PRM1[nMR+nL][0];      // icp max
    dm2.cleanup(splinecore1,n);
    dm2.cleanup(eos1,nrows1);
    dm2.cleanup(core1,npoints);
    
    }
    #pragma omp section
    {
    // Proposed change
    BA = props[0]*DATA[0][3]; kf = props[1]*DATA[1][3]; mstar = props[2]*DATA[2][3]*mNuc; 
    K = props[3]*DATA[3][3]; J = props[4]*DATA[4][3]; L = props[5]*DATA[5][3]; h = props[6]*DATA[6][3];
    p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
    cout << BA << "  " << kf << "  " << mstar << "  " << K << "  " << J << "  " << L << "  " << h << endl;
    nuc1.get_parameters(BA,p0,J,mstar,K,L,h,newparams,false,flag);     // solve for parameters given bulk properties
    nuc1.get_EOS(newparams,core2,npoints,false,unstable);     // get the EOS (nb,en,pr,dpde,mub,dPdrho,Keff)
    n = nuc1.splineeos(core2, npoints, 20, splinecore2,false);        // interpolate the EOS with 10 additional points per gap (for speed)
    nrows2 = nuc1.ThermalCrust(crust,splinecore2,eos2,n,nrowscrust,false,nbcol,prcol,keffcol);
    nm1.pretovconv(eos2, encol, prcol, cv, cv, nrows2);         // convert the EOS to dimensionless units
    lklp = nm1.MCMCtovpoint(1e-4,eos2,nrows2,ncols,encol,prcol,dpdecol,100,5.0*cv,nMR,nL,NICER,LIGO,PRM2);    // get 1.4 solar mass star radius (SEG FAULT??)
    nm1.RloveMCMC(eos2,dpdecol,encol,prcol,1e-5,ILQR2,PRM2[nMR][0],nrows2,ncols-1);        // get TD for 1.4
    maxp = PRM2[nMR+nL][2]; // maximum mass
    R14TDtheoryp = ILQR2[1];    // tidal defermability

    // store astro data
    for (int j=0; j<nMR; ++j) {
        props[7+j] = conv2.rnonetokm(PRM2[j][1]);
        props[7+nMR+j] = PRM2[j][0];
    }
    props[7+nMR*2] = conv2.rnonetokm(ILQR2[3]);         // store the 1.4 radius
    props[7+nMR*2+1] = R14TDtheoryp;         // store the 1.4 TD
    props[7+nMR*2+2] = maxp;      // Mmax
    props[7+nMR*2+3] = PRM2[nMR+nL][0];      // icp max
    dm2.cleanup(splinecore2,n);
    dm2.cleanup(eos2,nrows2);
    dm2.cleanup(core2,npoints);
    
    }
    }   // END PARALLEL

    // cleanup the obervables array
    dm2.cleanup(PRM1,nMR+nL+1); dm2.cleanup(PRM2,nMR+nL+1);

    // specify what error to use for the tidal deformability depending on which region you are in
    double R14TDerr;
    if (R14TDtheory0 < R14TDactual) {
        R14TDerr = R14TDerrL;
    } else {
        R14TDerr = R14TDerrU;
    }

    lkl0 = lkl0*exp(-0.5*pow((R14TDactual - R14TDtheory0)/R14TDerr,2.0));
    lklp = lklp*exp(-0.5*pow((R14TDactual - R14TDtheoryp)/R14TDerr,2.0));

    double XFTchisqp = stats.CFTlkl(XEFTdata,nrowsXEFT,newparams);
    double XFTchisq0 = stats.CFTlkl(XEFTdata,nrowsXEFT,oldparams);
    double lklX = exp(XFTchisq0/2.0 - XFTchisqp/2.0);
    //res = prior*lklp/lkl0;       // Just Astro
    res = prior*lklp/lkl0*lklX;     // Astro+XEFT
    return res;
}

// calibrate model to astrophysics and XEFT
// nMR is number of mass radius points from NICER
// nL is number of Ligo mass and TD points
// DATA matrix specifies the  current means of the model and appropriate scaling factors for the bulk properties along with good starting guesses and widths for MCMC
// output is of the form (BA, kf, mstar, K, J, L, h, NICER RADII, ..., NICER PRESSURES, ..., LIGO RADIUS, LIGO TD, MAXIMUM MASS, PRESSURE OF MMAX)
int MCMC :: Astro(int runs, string covdata, double** DATA, double** crust, int nrowscrust, int burnin, double** NICER, int nMR, double** XEFTdata, int nrowsXEFT, double** LIGO, int nL, string outname) {
    srand(time(0));
    int counts[7]; double stds[7]; double arate[7]; double agoal[7]; double paramtest[7];
    double r, a, BA, kf, J, mstar, K, L, h, p0;
    ofstream out(outname);
    ofstream lout("lastpoint.txt");
    int nprops = 7+nMR*2+nL*2+2; // store 7 bulks + (icp,radius) for each obervable + 1 TD + Mmax and icp of Mmax
    double *olds = new double[nprops];      // (bulks, radius1, radius2, ... , icp1, icp2, ... , TD1, TD2,  ... , Mmax, icpmax)
    double *cands = new double[nprops];     // (bulks, radius1, radius2, ... , icp1, icp2, ... , TD1, TD2,  ... , Mmax, icpmax)
    double cv = conv2.energyCONV(0,1);   // MeV/fm3 to unitless

    double** cov;
    dm2.importdata(covdata,cov);    // import covariance matrix

    // acceptance rate goals, starting point for bulk properties and their widths, acceptance counter, averages
    for (int i=0; i<7; ++i) {
        agoal[i] = 0.5;
        counts[i] = 0;
        olds[i] = DATA[i][0];
        stds[i] = DATA[i][1];
    }

    // rescale the bulk properties
    BA = olds[0]*DATA[0][3]; kf = olds[1]*DATA[1][3]; mstar = olds[2]*DATA[2][3]*mNuc; 
    K = olds[3]*DATA[3][3]; J = olds[4]*DATA[4][3]; L = olds[5]*DATA[5][3]; h = olds[6]*DATA[6][3];
    
    // make sure theres no issue with initial guess for params
    p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
    nuc1.get_parameters(BA,p0,J,mstar,K,L,h,paramtest,false,true);
    if (paramtest[6] < 0 || paramtest[5] < 0) {
        cout << "problem" << endl;
        dm2.cleanup(cov,7);
        dm2.cleanup(DATA,7);
        delete olds;
        delete cands;
        return 0;
    } 

    // Burn in phase
    for (int i=0; i<burnin; ++i) {
        cout << i << "  " << endl;
        // go through each bulk parameter and suggest a random change
        for (int j=0; j<7; ++j) {
            // candidate points = old points
            for (int k=0; k<nprops; ++k) {
                cands[k] = olds[k];
            }

            cands[j] = stats.rand_uniform(olds[j],stds[j]); // generate candidate point centered at old point within +- stdav
            a = targetAstro(cands,olds,cov,DATA,crust,nrowscrust,true,NICER,nMR,XEFTdata,nrowsXEFT,LIGO,nL);   // probability of accepance
            if (a>1.0) {    // probability shouldnt be greater than 1
                a = 1.0;        
            }

            r = 1.0*rand()/RAND_MAX;    // generate random number from 0 to 1
            if (r <= a) {   // accept candidate with probability a 
                olds[j] = cands[j];
                counts[j] = counts[j] + 1;      // add 1 to the accpetance count
                
                // carry over astro properties
                for (int k=7; k<nprops; ++k) {
                    olds[k] = cands[k];
                }
            }

            // monitor the rate every 50 points
            if ((i+1)%50 == 0) {
                
                arate[j] = 1.0*counts[j]/50.0;
                counts[j] = 0;
                
                if (arate[j] < agoal[j]) {
                    stds[j] = 0.9*stds[j];                 // if acceptance rate is too low then decrease the range
                } else if (arate[j] > agoal[j]) {
                    stds[j] = 1.1*stds[j];                 // if acceptance rate is too high then increase the range
                }
                
                // save arate, stds, olds every 50 iterations (in case of code breaking or computer freezing can be deleted with better computer)
                for (int k=0; k<7; ++k) {
                    lout << "arate: " << arate[k] << "  ";
                }
                lout << endl;
                for (int k=0; k<7; ++k) {
                    lout << "stds: " << stds[k] << "  ";
                }
                lout << endl;
                for (int k=0; k<7; ++k) {
                    lout << "olds: " << olds[k] << "  ";
                }
                lout << endl;
            }
        }
        // print the rate, width, and bulk values every run
        cout << "arate: ";
        for (int k=0; k<7; ++k) {
            cout << arate[k] << "  ";
        }
        cout << endl;
        cout << "stds: ";
        for (int k=0; k<7; ++k) {
            cout << stds[k] << "  ";
        }
        cout << endl;
        cout << "olds: ";
        for (int k=0; k<7; ++k) {
            cout << olds[k] << "  ";
        }
        cout << endl;
        
        // save MCMC runs
        out << setprecision(10);
        for (int k=0; k<7; ++k) {                   // print params
            out << olds[k]*DATA[k][3] << "  ";
        }
        for (int k=7; k<(7+nMR); ++k) {             // print nicer radii
            out << olds[k] << "  ";
        }
        for (int k=(7+nMR); k<(7+nMR*2); ++k) {     // print nicer icp
            out << olds[k]/cv << "  ";
        }
        for (int k=(7+nMR*2); k<(7+nMR*2+nL); ++k) {        // print ligo radii and TD
            out << olds[k] << "  " << olds[k+1] << "  ";
        }
        out << olds[nprops-2] << "  " << olds[nprops-1]/cv << endl;
    }

    // print out the bulk properties values and widths at end of burn in
    cout << "--------------------------------------------------------------------" << endl;
    out << "---------------------------------------------------------------------" << endl;
    out << olds[0] << "  " << olds[1] << "  " << olds[2] << "  " << olds[3] << "  " << olds[4] << "  " << olds[5] << "  " << olds[6] << endl;
    out << stds[0] << "  " << stds[1] << "  " << stds[2] << "  " << stds[3] << "  " << stds[4] << "  " << stds[5] << "  " << stds[6] << endl;
    cout << "--------------------------------------------------------------------" << endl;
    out << "---------------------------------------------------------------------" << endl;
    
    // MCMC RUN
    for (int i=0; i<runs; ++i) {
        // go through each bulk parameter and suggest a random change
        for (int j=0; j<7; ++j) {
            // candidate points = old points
            for (int k=0; k<nprops; ++k) {
                cands[k] = olds[k];
            }

            cands[j] = stats.rand_uniform(olds[j],stds[j]); // generate candidate point centered at old point within +- stdav
            a = targetAstro(cands,olds,cov,DATA, crust, nrowscrust,true,NICER,nMR,XEFTdata,nrowsXEFT,LIGO,nL);         // probability of acceptance
            if (a>1.0) {    // probability shouldnt be greater than 1
                a = 1.0;        
            }

            r = 1.0*rand()/RAND_MAX;    // generate random number from 0 to 1
            if (r <= a) {   // accept candidate with probability a 
                olds[j] = cands[j];
                counts[j] = counts[j] + 1;      // add 1 to the accpetance count
                
                for (int k=7; k<nprops; ++k) {
                    olds[k] = cands[k];
                }
            }
        }
        
        cout << "olds: ";
        for (int k=0; k<7; ++k) {
            cout << olds[k]*DATA[k][3] << "  ";
        }
        cout << endl;

        // save MCMC runs
        out << setprecision(10);
        for (int k=0; k<7; ++k) {
            out << olds[k]*DATA[k][3] << "  ";
        }
        for (int k=7; k<(7+nMR); ++k) {
            out << olds[k] << "  ";
        }
        for (int k=(7+nMR); k<(7+nMR*2); ++k) {
            out << olds[k]/cv << "  ";
        }
        for (int k=(7+nMR*2); k<(7+nMR*2+nL); ++k) {        // print ligo radii and TD
            out << olds[k] << "  " << olds[k+1] << "  ";
        }
        out << olds[nprops-2] << "  " << olds[nprops-1]/cv << endl;
    }
    
    // print out final rate
    cout << "arate: ";
    for (int k=0; k<7; ++k) {
        cout << 1.0*counts[k]/runs << "  ";
    }
    cout << endl;
    
    dm2.cleanup(cov,7);
    delete olds;
    delete cands;
    
    return 0;
}

// calculate the likelihood for XEFT given XEFT data and a parameter set
double statistics :: CFTlkl(double** XEFTdata, int npoints, double params[7]) {

    double** EOSNM;
    double en, dens, enX, err;
    int nbcol = 0;
    int encol = 3;  // energy per particle and not energy density en/dens - mNUc
    double chisq = 0;
    int n = nuc1.get_PNM_SNM_EOS(params,EOSNM,500,false,1.0);   // get the PNM EOS t=1
    for (int i=0; i<npoints; ++i) {
        dens = XEFTdata[i][0];
        enX = XEFTdata[i][1];
        err = XEFTdata[i][2];
        en = dm2.interpolate(n,5,EOSNM,dens,nbcol,encol,true);

        chisq = chisq + pow((enX-en)/err,2.0);
    }
    dm2.cleanup(EOSNM,n);

    return chisq;
}

// NEEDS A LOT OF COMMENTS
// obtain correlations between radii of stars and their pressure at several densities
double statistics :: correlations(string MCMCDATA, double** crust, int nc) {
    double** data;
    dm2.importdata(MCMCDATA,data);
    double BA, kf, mstar, K, J, L, h, p0;
    double params[7];    // array for the calculated parameters
    bool unstable = false;
    int npoints = 100;      // npoints for the eos
    int n, nrows, nNM, np;
    int encol = 1; int prcol = 2; int dpdecol = 3; int nbcol = 0; int Keffcol = 5;
    double cv = conv2.energyCONV(0,1);   // MeV/fm3 to unitless
    int ncols = 6;          // ncols on the EOS (nb, en, pr, dpde, dPdrho, Keff)
    int nr = dm2.rowcount(MCMCDATA);

    nr = 10000;
    int sample1 = 54; // (6.0-0.7)/0.1 + 1
    int sample2 = 54; // (6.0-0.7)/0.1 + 1
    int nstars = 10; // (1.9-1.0)/0.1 + 1
    int ccols = 1+nstars+sample1+sample2;

    double** NEWDATen;
    double** NEWDATdens;
    NEWDATen = new double*[nr];
    for (int i=0; i<nr; i++) {
        NEWDATen[i] = new double[ccols];
    }

    NEWDATdens = new double*[nr];
    for (int i=0; i<nr; i++) {
        NEWDATdens[i] = new double[ccols];
    }

    #pragma omp parallel num_threads(12)
    #pragma omp for ordered schedule(static,1) private(BA,kf,mstar,K,J,L,h,p0,params,n,nrows,np,nNM)
    for (int i=0; i<nr; ++i) {
        double** core; double** splinecore; double** eosc; double** prm2; double** EOSNM;
        BA = data[i][0]; kf = data[i][1]; mstar = data[i][2]*mNuc; K = data[i][3]; J = data[i][4]; L = data[i][5]; h = data[i][6];
        p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
        
        nuc1.get_parameters(BA,p0,J,mstar,K,L,h,params,false,true);     // solve for parameters given bulk properties
        nuc1.get_EOS(params,core,npoints,false,unstable);     // get the EOS (nb,en,pr,dpde,mub,dPdrho,Keff)
        n = nuc1.splineeos(core, npoints, 20, splinecore,false);    // interpolate the EOS with 20 additional points per gap (for speed) so the eos now has ~ 20*100 points
        nrows = nuc1.ThermalCrust(crust,splinecore,eosc,n,nc,false,nbcol,prcol,Keffcol);      // get the crust transition density adn add crust
        nm1.pretovconv(eosc, encol, prcol, cv, cv, nrows);     // convert the EOS to dimensionless units
        np = nm1.MCMCcorr(1e-4,eosc,nrows,ncols-1,encol,prcol,dpdecol,150,1.0*cv,prm2);
        nNM = nuc1.get_PNM_SNM_EOS(params,EOSNM,1000,false,1.0);

        
        NEWDATen[i][0] = data[i][12];
        NEWDATdens[i][0] = data[i][12];
        //cout << NEWDAT[i][0] << endl;
        double M = 1.0;
        double R;
        for (int j=1; j<(nstars+1); ++j) {
            R = dm2.interpolate(np,3,prm2,M,2,1,true);
            NEWDATen[i][j] = R;
            NEWDATdens[i][j] = R;
            //cout << M << "  " << NEWDAT[i][j] << endl;
            M = M+0.1;
        }

        double sat = 0.15;
        double frac = 0.7;
        double edens = frac*sat*mNuc;
        double dens = frac*sat;
        double pren, prdens;
        for (int j=(1+nstars); j<(1+nstars+sample1); ++j) {
            pren = dm2.interpolate(nrows,ncols-1,eosc,edens*cv,1,2,true);
            prdens = dm2.interpolate(nrows,ncols-1,eosc,dens,0,2,true);
            NEWDATen[i][j] = pren*1.0/cv;
            NEWDATdens[i][j] = prdens*1.0/cv;
            //cout << frac << "  " << NEWDAT[i][j] << endl;
            frac = frac+0.1;
            edens = frac*sat*mNuc;
            dens = frac*sat;
        }

        frac = 0.7;
        edens = frac*sat*mNuc;
        dens = frac*sat;
        for(int j=1+nstars+sample1; j<ccols; ++j) {
            //cout << frac << " check" << endl;
            pren = dm2.interpolate(nNM,5,EOSNM,edens,1,2,true);
            prdens = dm2.interpolate(nNM,5,EOSNM,dens,0,2,true);
            NEWDATen[i][j] = pren;
            NEWDATdens[i][j] = prdens;
            //cout << j << "  " << frac << "  " << NEWDAT[i][j] << endl;
            frac = frac+0.1;
            edens = frac*sat*mNuc;
            dens = frac*sat;                // adjust for density or energy density
        }
        
        dm2.cleanup(core,npoints);
        dm2.cleanup(splinecore,n);
        dm2.cleanup(eosc,nrows);
        dm2.cleanup(prm2,np);
        dm2.cleanup(EOSNM,nNM);

        #pragma omp ordered
        cout << i << "  " << BA << "  " << kf << "  " << mstar << "  " << K << "  " << J << "  " << L << "  " << h << endl;
    }
    
    // GET CORRELATION MATRICIES
    double CORRMATRIXBEMen[1+nstars][sample1];
    double CORRMATRIXPNMen[1+nstars][sample2];
    double CORRMATRIXBEMdens[1+nstars][sample1];
    double CORRMATRIXPNMdens[1+nstars][sample2];
    double corr;

    // Loop through each astrophysical property and density
    for (int i=0; i<(1+nstars); ++i) {
        // FOR BEM
        for (int j=0; j<sample1; ++j) {
            // FOR ENERGY DENSITY
            double rij = 0; double meani = 0; double meanj = 0; double si = 0; double sj = 0;

            // compute rij and means
            for (int k=0; k<nr; ++k) {
                rij = rij + NEWDATen[k][i]*NEWDATen[k][j+1+nstars];   
                meani = meani + NEWDATen[k][i];                    
                meanj = meanj + NEWDATen[k][j+1+nstars];         
            }
            meani = meani/nr;
            meanj = meanj/nr; 

            // Compute standard deviations
            for (int k=0; k<nr; ++k) {
                si = si + pow(NEWDATen[k][i] - meani,2.0);
                sj = sj + pow(NEWDATen[k][j+1+nstars] - meanj,2.0);
            }
            si = sqrt(si/(nr-1));
            sj = sqrt(sj/(nr-1));

            // calculate correlation
            corr = (rij - nr*meani*meanj)/((nr-1)*si*sj);
            CORRMATRIXBEMen[i][j] = corr;
            
            // DO SAME FOR DENSITY
            rij = 0; meani = 0; meanj = 0; si = 0; sj = 0;

            // compute rij and means
            for (int k=0; k<nr; ++k) {
                rij = rij + NEWDATdens[k][i]*NEWDATdens[k][j+1+nstars];
                meani = meani + NEWDATdens[k][i];
                meanj = meanj + NEWDATdens[k][j+1+nstars];
            }
            meani = meani/nr;
            meanj = meanj/nr;

            // Compute standard deviations
            for (int k=0; k<nr; ++k) {
                si = si + pow(NEWDATdens[k][i] - meani,2.0);
                sj = sj + pow(NEWDATdens[k][j+1+nstars] - meanj,2.0);
            }
            si = sqrt(si/(nr-1));
            sj = sqrt(sj/(nr-1));

            // calculate correlation
            corr = (rij - nr*meani*meanj)/((nr-1)*si*sj);
            CORRMATRIXBEMdens[i][j] = corr;

        }
        // FOR PNM
        for (int j=sample1; j<(sample1+sample2); ++j) {
            // FOR ENERGY DENSITY
            double rij = 0; double meani = 0; double meanj = 0; double si = 0; double sj = 0;

            // compute rij and means
            for (int k=0; k<nr; ++k) {
                rij = rij + NEWDATen[k][i]*NEWDATen[k][j+1+nstars];
                meani = meani + NEWDATen[k][i];
                meanj = meanj + NEWDATen[k][j+1+nstars];
            }
            meani = meani/nr;
            meanj = meanj/nr;

            // Compute standard deviations
            for (int k=0; k<nr; ++k) {
                si = si + pow(NEWDATen[k][i] - meani,2.0);
                sj = sj + pow(NEWDATen[k][j+1+nstars] - meanj,2.0);
            }
            si = sqrt(si/(nr-1));
            sj = sqrt(sj/(nr-1));

            // calculate correlation
            corr = (rij - nr*meani*meanj)/((nr-1)*si*sj);
            CORRMATRIXPNMen[i][j-sample1] = corr;

            // FOR DENSITY
            rij = 0; meani = 0; meanj = 0; si = 0; sj = 0;

            // compute rij and means
            for (int k=0; k<nr; ++k) {
                rij = rij + NEWDATdens[k][i]*NEWDATdens[k][j+1+nstars];
                meani = meani + NEWDATdens[k][i];
                meanj = meanj + NEWDATdens[k][j+1+nstars];
            }
            meani = meani/nr;
            meanj = meanj/nr;

            // Compute standard deviations
            for (int k=0; k<nr; ++k) {
                si = si + pow(NEWDATdens[k][i] - meani,2.0);
                sj = sj + pow(NEWDATdens[k][j+1+nstars] - meanj,2.0);
            }
            si = sqrt(si/(nr-1));
            sj = sqrt(sj/(nr-1));

            // calculate correlation
            corr = (rij - nr*meani*meanj)/((nr-1)*si*sj);
            CORRMATRIXPNMdens[i][j-sample1] = corr;
        }
    } 

    // FIND MAX CORRELATIONS
    double frac;
    ofstream out("CORR_BEM_en.txt");
    ofstream fout("CORR_BEM_dens.txt");
    frac = 0.7;
    for (int i=0; i<sample1; ++i) {
        out << scientific << setprecision(4) << frac << "  ";
        fout << scientific << setprecision(4) << frac << "  ";
        frac = frac+0.1;
        for (int j=0; j<nstars; ++j) {
            out << scientific << setprecision(4) << CORRMATRIXBEMen[j][i] << "  ";
            fout << scientific << setprecision(4) << CORRMATRIXBEMdens[j][i] << "  ";
        }
        out << scientific << setprecision(4) << CORRMATRIXBEMen[nstars][i] << endl;
        fout << scientific << setprecision(4) << CORRMATRIXBEMdens[nstars][i] << endl;
    }

    ofstream lout("CORR_PNM_en.txt");
    ofstream mout("CORR_PNM_dens.txt");
    frac = 0.7;
    for (int i=0; i<sample2; ++i) {
        lout << scientific << setprecision(4) << frac << "  ";
        mout << scientific << setprecision(4) << frac << "  ";
        frac = frac+0.1;
        for (int j=0; j<nstars; ++j) {
            lout << scientific << setprecision(4) << CORRMATRIXPNMen[j][i] << "  ";
            mout << scientific << setprecision(4) << CORRMATRIXPNMdens[j][i] << "  ";
        }
        lout << scientific << setprecision(4) << CORRMATRIXPNMen[nstars][i] << endl;
        mout << scientific << setprecision(4) << CORRMATRIXPNMdens[nstars][i] << endl;
    }
    
    dm2.cleanup(NEWDATen,nr);
    dm2.cleanup(NEWDATdens,nr);
    nr = dm2.rowcount(MCMCDATA);
    dm2.cleanup(data,nr);
    
    return 0;

}

// get the mean of an array for a specific column
double statistics :: mean(double** array, int nrows, int col) {
    double res = 0;
    for (int i=0; i<nrows; ++i) {
        res = res + array[i][col];
    }
    res = res/nrows;
    return res;
}

// get the std dev of an array for a specific column
double statistics :: stddev(double** array, int nrows, int col) {
    double res = 0;
    for (int i=0; i<nrows; ++i) {
        res = res + pow(array[i][col]-mean(array,nrows,col),2.0);
    }
    res = sqrt(res/(nrows-1));
    return res;
}

// Get all the predicted mass radius points for an MCMC run
// Will take each parameter set and calculate MR points which can then be histogramed in a 2D space
double statistics :: MRband(string MCMCDATA, double** crust, int nc) {
    double** data;
    dm2.importdata(MCMCDATA,data);
    double MR[2];
    double BA, kf, mstar, K, J, L, h, p0, dens, cs2, icp;
    double params[7];    // array for the calculated parameters
    bool unstable = false;
    int npoints = 100;      // npoints for the eos
    int nmrs = 300; // npoints for mr
    int n, nrows, x, ks, index, space, limit;
    int encol = 1; int prcol = 2; int dpdecol = 3; int nbcol = 0; int Keffcol = 5;
    double cv = conv2.energyCONV(0,1);   // MeV/fm3 to unitless
    double pr0 = 2.0*cv;
    
    int ncols = 7;          // ncols on the EOS (nb, en, pr, dpde, mub, dPdrho, Keff)
    int nr = dm2.rowcount(MCMCDATA);
    ofstream out("MR_BAND.txt");

    #pragma omp parallel num_threads(12)
    #pragma omp for ordered schedule(static,1) private(BA,kf,mstar,K,J,L,h,p0,params,n,nrows,ks,x,dens,cs2,icp,index,space,limit,MR)
    for (int i=0; i<nr; ++i) {
        double** core; double** splinecore; double** eosc; double** dens_cs2_r_m;
        BA = data[i][0]; kf = data[i][1]; mstar = data[i][2]*mNuc; K = data[i][3]; J = data[i][4]; L = data[i][5]; h = data[i][6];
        p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
        nuc1.get_parameters(BA,p0,J,mstar,K,L,h,params,false,true);     // solve for parameters given bulk properties
        nuc1.get_EOS(params,core,npoints,false,unstable);     // get the EOS (nb,en,pr,dpde,dPdrho,Keff)
        n = nuc1.splineeos(core, npoints, 20, splinecore,false);    // interpolate the EOS with 20 additional points per gap (for speed) so the eos now has ~ 20*100 points
        nrows = nuc1.ThermalCrust(crust,splinecore,eosc,n,nc,false,nbcol,prcol,Keffcol);      // get the crust transition density and add crust
        nm1.pretovconv(eosc, encol, prcol, cv, cv, nrows);     // convert the EOS to dimensionless units
        
        //##############################
        index = dm2.findvalue(eosc,nrows,ncols,pr0,prcol,1.0);        // find initial pressure row and start here
        space = ceil(1.0*(nrows-(index+1.0))/nmrs);              // get spacing given npoints and starting pressure
        if (space < 2) {
            cout << "space to be sampled is to small for the number of samples requested: " << nrows-(index+1.0) << "/" << npoints << endl;
            exit(0);
        }
        limit = floor((nrows-(index+1.0))/space);
        
        // create array to tabulate each dens, dpde, radius, mass
        dens_cs2_r_m = new double*[limit];
        for (int j=0; j<limit; j++) {
            dens_cs2_r_m[j] = new double[4];
        }
        dm2.zero(dens_cs2_r_m,limit,4);
        
        for (int j=0; j<limit; ++j) {
            x = index + space*j;
            icp = eosc[x][prcol];
            dens = eosc[x][nbcol];                // central pressure to integrate
            cs2 = eosc[x][dpdecol];
            nm1.MCMCtovsolve(icp,1e-4,eosc,nrows,ncols,encol,prcol,dpdecol,MR);
            dens_cs2_r_m[j][0] = dens;                    // store the dens
            dens_cs2_r_m[j][1] = cs2;
            dens_cs2_r_m[j][2] = MR[1];                  // store the radii
            dens_cs2_r_m[j][3] = MR[0];                  // store the masses

            // only save points until the maximum mass is reached
            if (j>0) {
                if (dens_cs2_r_m[j][3] < dens_cs2_r_m[j-1][3]) {
                    break;
                }
            }
            
        }

        ks=0;
        // output the iteration number, the point number for the individual parameter set run, dens, radius, mass, cs^2
        #pragma omp ordered
        while (dens_cs2_r_m[ks][2] > 0 || ks==limit) {
            out << setprecision(8) << i << "  " << ks << "  " << limit << "  " << dens_cs2_r_m[ks][0] << "  " << conv2.rnonetokm(dens_cs2_r_m[ks][2]) << "  " << dens_cs2_r_m[ks][3] << "  " << dens_cs2_r_m[ks][1] << endl;
            ks = ks+1;
        }
        cout << i << endl;

        // cleanup the arrays
        dm2.cleanup(eosc,nrows);
        dm2.cleanup(dens_cs2_r_m,limit);
        dm2.cleanup(core,npoints);
        dm2.cleanup(splinecore,n);
    } 
    
    return 0;
}

// generate theoretical errors for the energy per particle at a given density
double statistics :: EOSband(string MCMCDATA) {
    double** data;
    dm2.importdata(MCMCDATA,data);  // import the MCMC parameter sample
    double BA, kf, mstar, K, J, L, h, p0, dens;
    double params[7];    // array for the calculated parameters
    bool unstable = false;
    int npoints = 300;      // npoints for the eos
    int prcol = 2; int encol = 1;
    double cv = conv2.energyCONV(0,1);   // MeV/fm3 to unitless
    
    int ncols = 5;          // ncols on the EOS (nb, en, pr, E/N, mstar)
    int nr = dm2.rowcount(MCMCDATA);
    ofstream out("EOS_BAND.txt");
    double eos_band[npoints][3]; // dens, mean EN, std EN;

    for (int j=0; j<npoints; ++j) {
        eos_band[j][0] = 0.0;
        eos_band[j][1] = 0.0;
        eos_band[j][2] = 0.0;
    }

    #pragma omp parallel num_threads(12)
    #pragma omp for ordered schedule(static,1) private(BA,kf,mstar,K,J,L,h,p0,params,dens)
    for (int i=0; i<nr; ++i) {
        double** eos;
        BA = data[i][0]; kf = data[i][1]; mstar = data[i][2]*mNuc; K = data[i][3]; J = data[i][4]; L = data[i][5]; h = data[i][6];
        p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
        nuc1.get_parameters(BA,p0,J,mstar,K,L,h,params,false,true);     // solve for parameters given bulk properties
        //nuc1.get_PNM_SNM_EOS(params,eos,npoints,false,1.0); // get PNM EOS
        nuc1.get_EOS(params,eos,npoints,false,false); // get NSM EOS
        for (int j=0; j<npoints; ++j) {
            eos_band[j][1] = eos_band[j][1] + eos[j][prcol];    // sum all the pressures for a given density
        }

        dm2.cleanup(eos,npoints);
        cout << i << endl;
    }

    // calculate mean
    for (int j=0; j<npoints; ++j) {
        eos_band[j][1] = eos_band[j][1]/nr;
    }

    // calculate std
    #pragma omp parallel num_threads(12)
    #pragma omp for ordered schedule(static,1) private(BA,kf,mstar,K,J,L,h,p0,params,dens)
    for (int i=0; i<nr; ++i) {
        double** eos;
        BA = data[i][0]; kf = data[i][1]; mstar = data[i][2]*mNuc; K = data[i][3]; J = data[i][4]; L = data[i][5]; h = data[i][6];
        p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
        nuc1.get_parameters(BA,p0,J,mstar,K,L,h,params,false,true);     // solve for parameters given bulk properties
        
        //nuc1.get_PNM_SNM_EOS(params,eos,npoints,false,1.0);
        nuc1.get_EOS(params,eos,npoints,false,false);
        for (int j=0; j<npoints; ++j) {
            eos_band[j][2] = eos_band[j][2] + pow(eos[j][prcol]-eos_band[j][1],2.0);
            eos_band[j][0] = eos[j][encol]; // fill energy density
        }

        dm2.cleanup(eos,npoints);
        cout << i << endl;
    }

    // calculate std
    for (int j=0; j<npoints; ++j) {
        eos_band[j][2] = sqrt(eos_band[j][2]/(nr-1));
    }

    for (int i=0; i<npoints; ++i) {
        out << setprecision(8) << eos_band[i][0] << "  " << eos_band[i][1] << "  " << eos_band[i][2] << endl;
    }

    return 0;
}