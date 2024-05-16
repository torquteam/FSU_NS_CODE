class EOS {
public:
    void whitedwarf(double nucNum);
    // FSU Model
    double get_en(double kf, double t, double gss, double gsoms2, double gww, double gwomw2, double gpp, double gpomp2, double b, double c, double h, double lambda);
    double get_pr(double kf, double t, double gss, double gsoms2, double gww, double gwomw2, double gpp, double gpomp2, double b, double c, double h, double lambda);
    double get_K(double kf, double gss, double gwomw2, double gsoms2, double b, double c, double h, double gww);
    double get_J(double kf, double gss, double gpomp2, double gww, double lambda);
    double get_L(double kf, double gss, double gsoms2, double gww, double gwomw2, double gpomp2, double h, double lambda, double b, double c);
    int get_parameters(double BA, double p0, double J, double mstar, double K, double L, double h, double params[7], bool print, bool flag);
    int get_EOS(double params[7], double** &eos, int npoints, bool print, bool unstable);
    void get_bulkproperties(string eosSNM, int nbcol, int encol, int prcol, int BAcol, int mstrcol, double gsoms2, double gwomw2, double gpomp2, double b, double c, double h, double lambda);
    
    double get_t_betaeq(double kf, double gsoms2, double gwomw2, double gpomp2, double b, double c, double lambda, double h);
    int splineeos(double** eosarr, int nrows, int npoints, double** &neweos, bool print);
    double FSUgssnewton(double dens, double gsoms2, double b, double c, double t, double eps);
    double chneutralityeq(double kf, double gsoms2, double gwomw2, double gpomp2, double b, double c, double lambda, double h, double t);
    double FSUgwwnewton(double gwomw2, double gpomp2, double lambda, double h, double dens, double t);
    double effK2(double kf, double gss, double gsoms2, double gww, double gwomw2, double gpomp2, double gpp, double h, double lambda, double b, double c, double t);
    int ThermalCrust(double** crust, double** core, double** &neweos, int nrowscore, int nrowscrust, bool print, int nbcol, int prcol, int Keffcol);
    int get_PNM_SNM_EOS(double params[7], double** &eos, int npoints, bool print, double t);
};
