#include "NumMethods.hpp"
#include "quarkEOS.hpp"
#include "nucEOS.hpp"
#include "Conversions.hpp"
#include "MCMC.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <omp.h>

using namespace std;

int main() {
    EOS nuc;
    data2 dm;
    nummeth nm;
    Convert conv;
    MCMC mcmc;
    statistics stats;

    // Constants
    const double mP = 939.0; //938.27231;    // Mass of proton (MeV)
    const double mN = 939.0; //939.56542052;
    const double mNuc = (mP+mN)/2.0;
    const double pi = 4.0*atan(1.0);

    auto start = chrono :: high_resolution_clock::now();     
    
    /*
    //---------------- Calculate White Dwarf EOS -----------------------
    // Calculate White Dwarf EOS
    double nucnum = 2.0;            // number of nucleons per electron
    nuc.whitedwarf(nucnum);       // calculate EoS (mue,nb,en,pr) in MeV/fm3
    */

    // -----------------Calculate Self Consistent Nuclear EOS
    double BA, kf, p0, mstar, K, J, L, h;
    double params[7];
    // FSU GOLD2
    /*
    BA = -16.28017866;                                           // Binding energy (MeV)
    kf = 1.306395766;
    p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);                      // saturation density (fm-3)                  
    J = 37.62302751;                                          // Symmetry energy at sat (MeV)
    mstar = 0.5927305581*939.0;                                   // Effective Mass (MeV)
    K = 237.9617935;                                             // Compressibility (MeV)
    L = 112.8424899;                                             // Derivative of Symmetry Energy at sat (MeV)
    h = 0.004274483658;                                            // Self interaction strength for w meson       
    nuc.get_parameters(BA,p0,J,mstar,K,L,h,params,true,true);              // calculate coupling constants FSU Model         
    */
    
    // FSU GARNET
    params[0] = 110.3492/pow(496.939,2.0);
    params[1] = 187.695/pow(782.5,2.0);
    params[2] = 192.927/pow(763.0,2.0);
    params[3] = 3.26/(2.0*939.0);
    params[4] = -0.003551/6.0;
    params[5] = 0.0235/6.0;
    params[6] = 0.043377*2.0;
    
    /*
    BA = -16.180293338575;
    kf = 1.3124559217007;
    mstar = 0.58164567357185*939.0;
    K = 228.79445043241;
    J = 30.89216970096;
    L = 55.792292008302;
    h = 0.0238121311552134/6.0;
    p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
    nuc.FSUparamsolve(BA,p0,J,mstar,K,L,h,params,true,true);              // calculate coupling constants FSU Model
    */

    
    double** COREEOS;
    int npoints = 100;
    nuc.get_EOS(params,COREEOS,npoints,true,true); // (nb,en,pr,dpde,mub,dpdrho,Keff)
    //nuc.get_PNM_SNM_EOS(params,COREEOS,npoints,true,0.0);      // (nb,en,pr,E/N,mstar)
    dm.cleanup(COREEOS,npoints);        // comment out if getting crust EOS or using COREEOS to get MR or ILQ
    
    /*
    // Add Crust (Optional)
    double** crust; double** EOS;
    dm.importdata("Initializing_files/CRUSTEOS.txt",crust);
    int nrows = dm.rowcount("Initializing_files/CRUSTEOS.txt");
    int n = nuc.ThermalCrust(crust,COREEOS,EOS,npoints,nrows,true,0,2,6);
    dm.cleanup(crust,nrows);
    dm.cleanup(COREEOS,npoints);
    dm.cleanup(EOS,n);          // comment out if using the EOS to calc MR
    */
    
    // --------------- Calculate Bulk Properties 
    //string fileSNM = "FSUEOS_SNM.txt";
    //nuc.FSUbulkprop(fileSNM,0,1,2,3,4,params[0],params[1],params[2],params[3],params[4],params[5],params[6]);

    /*
    // ---------------- Calculate Mass Radius Profile
    string mrfilename = "FSUGARNET_XEFT";               // output MR filename title (adds _MR)
    int encol = 1;                          // input en col
    int prcol = 2;                          // input pr col
    int dpdecol = 3;
    double MR[2];
    double cv = conv.energyCONV(0,1);   // MeV/fm3 to unitless
    // --------------------------------------------------
    
    nm.pretovconv(EOS, encol, prcol, cv, cv, n);                                 //convert to unitless en, pr
    double pr0 = 2.0*cv; // start at 2 MeV/fm3
    //nm.MCMCtovsolve(5.0*cv,1e-4,EOS,n,5,encol,prcol,dpdecol,MR);
    //cout << MR[0] << "  " << conv.rnonetokm(MR[1]) << endl;
    nm.multitov(1e-4, EOS, n, 5, encol, prcol, dpdecol, 1000, mrfilename,pr0);    // calculate mass radius profile, output is (en,pr,r,m) (mev/fm3)
    dm.cleanup(EOS,n);        // comment out if calc ILQ
    */

    /*
    //------------------------ To specify what mass of star for the Tidal Deformability
    double** MRarr;
    string RMfile = "FSUGOLD_RM.txt";
    double starmass = 1.40;
    dm.importdata(RMfile,MRarr);
    int nrowsRM = dm.rowcount(RMfile);
    double cpr = dm.interpolate(nrowsRM,4,MRarr,starmass,3,1,true);
    dm.cleanup(MRarr,nrowsRM);
    cout << cpr << " MeV/fm3" << endl;
    
    
    // ---------------------- Calculate Relativistic and Newtonian ILQ
    encol = 1;                          // input en col
    prcol = 2;                          // input pr col
    dpdecol = 3;
    int np = 100;                      // npoints for ILQ
    cv = conv.energyCONV(0,1);   // MeV/fm3 to unitless
    double ILQR[4];                   // (I,L,Q, R)
    //nm.pretovconv(EOS, encol, prcol, cv, cv, n);                     // convert to unitless en, pr if didnt convert already in MR
    // ----------------------------------------------------
    
    // Single ILQ
    double icp = cpr*cv;                                    // central density to calculate 
    cout << "Relativistic: " << endl;
    nm.RloveMCMC(EOS,dpdecol,encol,prcol,1e-5,ILQR,icp,n,5);                      // Calculate I, L, Q
    cout << "Newtonian: " << endl;
    //nm.Nlove("SingleICP.txt",2,3,0,4,1e-3,ILQ);
    dm.cleanup(EOS,n);
    */    
   /*
    double pr0 = 14.0*cv; // start at 1 MeV/fm3
    double prf = 104*cv; 
    nm.RILOVEQ(eosDAT2,nrows,ncols+1,npoints,ncols,encol,prcol,pr0,prf);
    dm.cleanup(eosDAT,nrows);
    dm.cleanup(eosDAT2,nrows);
   */
    /*
    // ---------------------- Markov Chain Monte Carlo Run
    // Astrophysical Observables/XEFT MCMC
    double** crust; double** init; double** NICER; double** XEFT; double** LIGO;
    dm.importdata("Initializing_files/CRUSTEOS.txt",crust);
    dm.importdata("Initializing_files/NICERDAT.txt",NICER);
    dm.importdata("Initializing_files/LIGODAT.txt",LIGO);
    dm.importdata("Initializing_files/startfileAstro_XEFT_FSUGOLD.txt",init);
    dm.importdata("Initializing_files/XEFTMCMC.txt",XEFT);
    int nc = dm.rowcount("Initializing_files/CRUSTEOS.txt");
    int ns = dm.rowcount("Initializing_files/NICERDAT.txt");
    int nx = dm.rowcount("Initializing_files/XEFTMCMC.txt");
    int nL = dm.rowcount("Initializing_files/LIGODAT.txt");
    mcmc.Astro(5000,"Initializing_files/invcovmatrix_FSUGOLD.txt",init, crust, nc,500,NICER,ns,XEFT,nx,LIGO,nL,"postXEFT_FSUGOLD_test.txt");
    dm.cleanup(crust,nc);
    dm.cleanup(init,7);
    dm.cleanup(NICER,ns);
    dm.cleanup(XEFT,nx);
    dm.cleanup(LIGO,nL);
    */
    /*
    // Get sample of the prior
    double** init;
    dm.importdata("startfileAstro_FSUGOLD.txt",init);
    mcmc.RMF(100000,init,"invcovmatrix_FSUGOLD.txt");
    dm.cleanup(init,7);
    */
    /*
    // Calculate Astrophysical priors
    double** crust; double** NICER; double** LIGO;
    dm.importdata("CRUSTEOS.txt",crust);
    dm.importdata("NICERDAT.txt",NICER);
    dm.importdata("LIGODAT.txt",LIGO);
    int nc = dm.rowcount("CRUSTEOS.txt");
    int ns = dm.rowcount("NICERDAT.txt");
    int nL = dm.rowcount("LIGODAT.txt");
    mcmc.Observables("prior.txt",crust,nc,NICER,ns,LIGO,nL,"prior_FSUGOLD_p2.txt",50000);
    dm.cleanup(crust,nc);
    dm.cleanup(NICER,ns);
    dm.cleanup(LIGO,nL);
    */
    
    
    // Get Correlations
    //double** crust;
    //dm.importdata("CRUSTEOS.txt",crust);
    //int nc = dm.rowcount("CRUSTEOS.txt");
    //stats.correlations("postXEFT_FSUGARNET.txt",crust,nc);
    //stats.MRband("postXEFT_FSUGOLD.txt",crust,nc);
    //stats.EOSband("Python_and_outputs/postXEFT_FSUGARNET.txt");
    //dm.cleanup(crust,nc);
    

    // Get sample for posterior
    /*
    string file = "prior_FSUGARNET.txt";
    double** array;
    dm.importdata(file,array);
    ofstream out("prior_FSUGARNET_conv.txt");

    int j;
    for (int i=1; i<1001; ++ i) {
        j = i*20 - 1;
        out << setprecision(10) << array[j][0] << "  " << array[j][1] << "  " << array[j][2] << "  " << array[j][3] << "  " << array[j][4] << "  " << array[j][5] << "  " << array[j][6]*6.0 << endl;
    }   
    dm.cleanup(array,20000);
    */
    /*
    // output couplings to compute finite nuclei ( gsoms2, gwomw2, gpomp2, b, c, h, lambda, ms)
    double FSUGOLD2_ms = 497.47934496;
    double FSUGOLD2_gs2 = 108.09284824407266936231;
    double FSUGARNET_ms = 496.939;
    double FSUGARNET_gs2 = 110.3492;
    string MCMCfile = "Python_and_outputs/postXEFT_FSUGARNET.txt";
    double** bulkarray;
    dm.importdata(MCMCfile,bulkarray);
    int nrows = dm.rowcount(MCMCfile);
    ofstream out("FSUGARNET_XEFT_params.txt");
    for (int i=0; i<nrows; ++i) {
        BA = bulkarray[i][0]; kf = bulkarray[i][1];
        mstar = bulkarray[i][2]*mNuc; K = bulkarray[i][3];
        J = bulkarray[i][4]; L = bulkarray[i][5]; h = bulkarray[i][6];
        p0 = 2.0/(3.0*pow(pi,2.0))*pow(kf,3.0);
        cout << i << endl;
        nuc.get_parameters(BA,p0,J,mstar,K,L,h,params,false,true);
        // print out couplings and the rescaled sigma mass 
        out << setprecision(10) << params[0] << "  " << params[1] << "  " << params[2] << "  " << params[3] << "  " << params[4] << "  " << params[5] << "  " << params[6] 
            << "  " << sqrt((params[0]*pow(FSUGARNET_ms,2.0)/FSUGARNET_gs2)*pow(FSUGARNET_ms,2.0)) << endl;
    }
    dm.cleanup(bulkarray,nrows);
    */
   
    auto stop = chrono :: high_resolution_clock::now();
    auto duration = chrono :: duration_cast<chrono :: milliseconds>(stop - start);
    cout << setprecision(5) << duration.count()/1000.0 << "s" << endl;
    return 0;
    
}