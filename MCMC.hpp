class statistics {
    public:
        double rand_normal(double mean, double stddev);
        double rand_uniform(double cntr, double w);
        double CFTlkl(double** XEFTdata, int npoints, double params[7]);
        double mean(double** array, int nrows, int col);
        double stddev(double** array, int nrows, int col);
        double correlations(string MCMCDATA, double** crust, int nc);
        double MRband(string MCMCDATA, double** crust, int nc);
        double EOSband(string MCMCDATA);
};

class MCMC {
    public:
        int LDM(int npoints, string nucdata, int BAcol, int Acol, int Zcol);
        int RMF(int npoints, double** INIT, string covdata);
        double Observables(string paramset, double** crust, int nrowscrust, double** NICER, int nMR, double** LIGO, int nL, string outname, int nruns);
        int Astro(int runs, string covdata, double** DATA, double** crust, int nrowscrust, int burnin, double** NICER, int nMR, double** XEFTdata, int nrowsXEFT, double** LIGO, int nL, string outname);
};