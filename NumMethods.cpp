//
//  NumMethods.cpp
//  NS Code
//
//  Created by Marc Salinas on 1/29/19.
//  Copyright © 2019 Marc Salinas. All rights reserved.
//

#include "NumMethods.hpp"
#include "Conversions.hpp"
#include "nucEOS.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <omp.h>
data2 dm;
Convert conv;
EOS nuc;
using namespace std;

const double pi = 4*atan(1);
const double mE = 0.511;        // Mass of electron (MeV)
const double mU = 2.3;          // Mass of up quark (MeV)
const double mD = 4.8;          // Mass of down quark (MeV)
const double mS = 95.0;         // Mass of strange quark (MeV)
const double mP = 938.27231;    // Mass of proton (MeV)
const double mMU = 105.6583745; // Mass of muon (MeV)

// Data Manipulation class
//*********************************************************************************
//*********************************************************************************

// find the number of rows in data
int data2 :: rowcount(string txtfile) {
    int count = 0;                  // Set initial count to 0
    string line;                    // Initialize variable to store line's string temporarily
    ifstream file(txtfile);         // Read the txt file
    if (!file) {                    // Throw error if file isn't found
        cout << "RowCount: Error opening file from path: " << txtfile << endl;
        file.close();
        exit(0);
    }
    
    // Use getline function from ifstream to count non empty lines
    while (getline(file, line)) {
        count++;
    }
    file.close();
    return count; // Return number of lines
}

//----------------------------------------------------------------------------------------

// identify the type of file
string data2 :: filetype(string txtfile) {
    
    ifstream file;
    string line;
    size_t pos;
    file.open(txtfile);
    
    if(!file){                      // Throw error if file isn't found
        cout << "Filetype: Error opening file from path: " << txtfile << endl;
        file.close();
        exit(0);
    }
    
    getline(file,line);             // get line from file
    pos=line.find(",");             // search
    if(pos!=string::npos) {         // string::npos is returned if string is not found
        file.close();
        return "CSV";
    } else {
        file.close();
        return "TXT";
    }
}

//----------------------------------------------------------------------------------------

// find number of columns in data
int data2 :: colcount(string txtfile) {
    string line;
    ifstream file(txtfile);
    if (!file) {                    // Throw error if file isn't found
        cout << "Colcount: Error opening file from path: " << txtfile << endl;
        file.close();
        exit(0);
    } else {
        if (filetype(txtfile) == "TXT") {
            int numcols=0;
            getline(file,line);
            stringstream iss(line);
            do {
                std::string sub;
                iss >> sub;
                if (sub.length()) {
                    ++numcols;
                }
            }
        while(iss);
        return numcols;
        } else {
            int numcols = 0;
            while(getline(file,line)) {
                stringstream linestream(line);
                string value;
                while(getline(linestream,value,',')) {
                    ++numcols;
                }
                return numcols;
            }
        }
    }
    return 0;
}

//----------------------------------------------------------------------------------------

// import data into an array
void data2 :: importdata(string txtfile, double ** &array) {
    int numrows = rowcount(txtfile);            // Get the number of rows
    int numcols = colcount(txtfile);            // Get the number of cols
    string type = filetype(txtfile);            // Get the file type
    ifstream in(txtfile);
    
    // Initialize array size
    array = new double*[numrows];
    for (int i = 0; i<numrows; i++) {
        array[i] = new double[numcols];
    }
    
    if (type == "TXT") {                        // Import data for txt
        for (int j=0; j<numrows; ++j) {
            for (int i=0; i<numcols; ++i) {
                in >> array[j][i];
            }
        }
        in.close();
    }
    if (type == "CSV") {                        // Import data for csv
        ifstream file;
        string line;
        int i=0;
        int j=0;
        file.open(txtfile);
        while(getline(file,line)) {
            stringstream linestream(line);
            string value;
            while(getline(linestream,value,',')) {
                //cout << i << "  " << j << endl;
                array[i][j] = stod(value);
                j=j+1;
            }
            i=i+1;
            j=0;
        }
        in.close();
    }
};

//----------------------------------------------------------------------------------------

// Print out an array
void data2 :: print(double** &array, int numrows, int numcols, bool file, string filename) {
    if (file == true) {
        ofstream out(filename);
        for (int i=0; i<numrows; ++i) {
            for (int j=0; j<(numcols-1); ++j) {
                out << scientific << setprecision(15) << array[i][j] << "  ";
            }
        out << scientific << setprecision(15) << array[i][numcols-1] << endl;
        }
    } else {
        for (int i=0; i<numrows; ++i) {
            for (int j=0; j<numcols; ++j) {
                cout << array[i][j] << "  ";
            }
        cout << endl;
        }
    }
};
//----------------------------------------------------------------------------------------

// Interpolate a value in an array
double data2 :: interpolate(int numrows, int numcols, double** array, double point, int pointcol, int ycol, bool sortflag) {
    
    // Error checking
    if (ycol >= numcols || ycol < 0) {                                              // check to see column exists
        cout << "Interpolation error: column is out of bounds" << endl;
        print(array,numrows,numcols,true,"debug.txt");
        exit(0);
    }

    if (numrows <= 1) {                                                           // check to see if interpolation is possible
        cout << "Interpolation error: not enough points to interpolate" << endl;
        exit(0);
    }

    // Sort data if unsorted
    if (sortflag == false) {
        sortasc(array, pointcol, numrows, numcols);
    }

    if (point < array[0][pointcol] || point > array[numrows-1][pointcol]) {       // check to see if point is in bounds
        cout << "Interpolation error: point is not in bounds: " << array[0][pointcol] << ", " << point << ", " << array[numrows-1][pointcol] << endl;
        exit(0);
    }

    // Set up bisection for the upper and lower limits of the point
    int mid = floor(numrows/2);
    double sol, slope, yint;
    int midpoint = 0;
    int upperbound = -1;
    int lowerbound = -1;
    int n;
    
    // base case for an array of 2 or 3 rows
    if (mid == 1) {
        slope = (array[numrows-1][ycol] - array[0][ycol])/(array[numrows-1][pointcol] - array[0][pointcol]);
        yint = array[numrows-1][ycol] - slope*array[numrows-1][pointcol];
        sol = slope*point + yint;
    }
    
    // Bisection to find a bounds
    n = mid;
    while(mid != 0) {
        midpoint = n - 1;   // rownumber to row location in array
        mid = floor(mid/2);
        //cout << mid << " " << array[midpoint][pointcol] << endl; //(For Debug only)
        if (point >= array[midpoint][pointcol]) {
            n = n + mid;
        } else {
            n = n - mid;
        }
    }
    
    // identify the bound as upper or lower
    if (array[midpoint][pointcol] >= point || midpoint == (numrows - 1)) {
        upperbound = midpoint;
        lowerbound = midpoint -1;
        while (array[lowerbound][pointcol] > point) {
            if (lowerbound - 1 == 0) {
                lowerbound = 0;
                upperbound = 1;
                break;
            } else {
                lowerbound = lowerbound - 1;
                upperbound = upperbound - 1;
            }
        }
    } else if (array[midpoint][pointcol] < point || midpoint == 0) {
        lowerbound = midpoint;
        upperbound = midpoint + 1;
        while (array[upperbound][pointcol] < point) {
            if (upperbound + 1 == numrows - 1) {
                lowerbound = numrows - 2;
                upperbound = numrows - 1;
                break;
            } else {
                lowerbound = lowerbound + 1;
                upperbound = upperbound + 1;
            }
        }
    }
    
    
    
    if (upperbound == -1 || lowerbound == -1 || lowerbound > upperbound) {
        cout << "Interpolation error unknown" << endl;
        exit(0);
    }
    
    // Interpolate
    slope = (array[upperbound][ycol] - array[lowerbound][ycol])/(array[upperbound][pointcol] - array[lowerbound][pointcol]);
    yint = array[upperbound][ycol] - slope*array[upperbound][pointcol];
    sol = slope*point + yint;

    return sol;
}

void data2 :: cubicspline(double**arr, double**&spline, int nrows, int xcol, int ycol) {
    int n = nrows-2;
    double a[n-1];
    double b[n];
    double c[n-1];
    double d[n];
    double cp[n-1];
    double dp[n];

    spline = new double*[nrows];
    for (int i = 0; i<nrows; i++) {
        spline[i] = new double[3];
    }

    for (int j=0; j<(n-1); ++j) {
        a[j] = 1.0/6.0*(arr[j+2][xcol] - arr[j+1][xcol]);
        b[j] = 1.0/3.0*(arr[j+2][xcol] - arr[j][xcol]);
        c[j] = 1.0/6.0*(arr[j+2][xcol] - arr[j+1][xcol]);
        d[j] = (arr[j+2][ycol] - arr[j+1][ycol])/(arr[j+2][xcol] - arr[j+1][xcol]) - (arr[j+1][ycol] - arr[j][ycol])/(arr[j+1][xcol] - arr[j][xcol]);
    }
    b[n-1] = 1.0/3.0*(arr[n+1][xcol] - arr[n-1][xcol]);
    d[n-1] = (arr[n+1][ycol] - arr[n][ycol])/(arr[n+1][xcol] - arr[n][xcol]) - (arr[n][ycol] - arr[n-1][ycol])/(arr[n][xcol] - arr[n-1][xcol]);

    
    cp[0] = c[0]/b[0];
    dp[0] = d[0]/b[0];
    for (int j=1; j<(n-1); ++j) {
        cp[j] = c[j]/(b[j] - a[j-1]*cp[j-1]);
        dp[j] = (d[j] - a[j-1]*dp[j-1])/(b[j]-a[j-1]*cp[j-1]);
    }
    dp[n-1] = (d[n-1] - a[n-2]*dp[n-2])/(b[n-1]-a[n-2]*cp[n-2]);

    spline[nrows-1][2] = 0.0;
    spline[n][2] = dp[n-1];

    for (int i=n-1; i>0; --i) {
        spline[i][2] = dp[i-1] - cp[i-1]*spline[i+1][2];
    }

    spline[0][2] = 0.0;

    for (int i=0; i<nrows; ++i) {
        spline[i][0] = arr[i][xcol];
        spline[i][1] = arr[i][ycol];
    }
}

double data2 :: splinecalc(double** spline, int x1, int x2, double x) {
    double a,b,c,d,y;
    
    a = (spline[x2][0] - x)/(spline[x2][0] - spline[x1][0]);
    b = 1.0-a;
    c = 1.0/6.0*(pow(a,3.0)-a)*pow(spline[x2][0]-spline[x1][0],2.0);
    d = 1.0/6.0*(pow(b,3.0)-b)*pow(spline[x2][0]-spline[x1][0],2.0);
    y = a*spline[x1][1] + b*spline[x2][1] + c*spline[x1][2] + d*spline[x2][2];
    //cout << spline[x1][0] << "  " << spline[x2][0] << "  " << a << "  " << b << "  " << c << "  " << d << "  " << y << endl;

    return y;
}

//----------------------------------------------------------------------------------------

// Bubble sort ascending order
bool data2 :: sortasc(double ** array, int col, int numrows, int numcols) {
    
    if (col >= numcols || col < 0) {
        cout << "Sorting error: column is out of bounds" << endl;           // Throw error if col is out of bounds
        exit(0);
    }
    
    // set up a temp array to store data
    double temp[numcols];

    int count;
    int numruns = 0;
    bool sorted = false;                                                    // Assume unsorted
    while (sorted == false) {                                               // Run while unsorted
        count = 0;                                                          // Count how many changes were made
        for (int i=0; i<(numrows-1); ++i) {
            for (int j=0; j<numcols; ++j) {                                 // store current row in temp array
                temp[j] = array[i+1][j];
            }
            if (array[i][col] > array[i+1][col]) {                          // compare two columns to see if they are in desc order
                //cout << i << endl;
                for (int j=0; j<numcols; ++j) {                             // Swap them if they are in desc order
                    array[i+1][j] = array[i][j];
                    array[i][j] = temp[j];
                }
                count = count + 1;                                          // add to count
            }
        }
        if (count == 0) {                                                   // set sorted flag to true if no more swaps are necessary
            sorted = true;
        } else {
            numruns = numruns + 1;                                          // number of swaps overall counter
        }
    }
    //cout << "Number of sort runs in col "<< col << ": " << numruns << endl;
    if (numruns == 0) {
        return true;                                                        // return true if array was sorted without any runs
    } else {
        return false;                                                       // return false if sorting was necessary
    }
}

//----------------------------------------------------------------------------------------

// cleanup pointer array
void data2 :: cleanup(double** &array, int nrows) {
    for (int i = 0; i<nrows; i++) {
        delete[] array[i];
    }
    delete[] array;
}

//----------------------------------------------------------------------------------------

// Find maximum of data array
double data2 :: findmax(double** array, int col, int numrows, int numcols) {

    if (col >= numcols || col < 0) {                                        // Throw error if col is out of bounds
        cout << "Find Max error: column is out of bounds" << endl;
        exit(0);
    }

    double max = array[0][col];                                             // find the max of an array by brute force
    for (int i=0; i<numrows; ++i) {
        if (array[i][col] > max) {
            max = array[i][col];
        }
    }
    return max;
}

//----------------------------------------------------------------------------------------

// Find minimum of data array
double data2 :: findmin(double** array, int col, int numrows, int numcols) {

    if (col >= numcols || col < 0) {                                    // Throw error if col is out of bounds
        cout << "Find min error: column is out of bounds" << endl;
        exit(0);
    }

    double min = array[0][col];                                         // Find the minimum by brute force
    for (int i=0; i<numrows; ++i) {
        if (array[i][col] < min) {
            min = array[i][col];
        }
    }
    return min;
}

// Find a given point in an array
int data2 :: findvalue(double **array, int nrows, int ncols, double point, int pcol, double tol) {

    if (pcol >= ncols || pcol < 0) {                                  // Throw error if column is out of bounds
        cout << "Find value error: column is out of bounds" << endl;
        exit(0);
    }

    double error = abs(point-array[0][pcol]);                           // set initial error
    double newerror;                                                 // error tolerance of point
    int irow = 0;
    for (int i=0; i<nrows; ++i) {
        newerror = abs(point-array[i][pcol]);                           // calculate how new error
        if (newerror < error) {                                         // smallest error is point within tolerance
            error = newerror;                                           // set the old error as new error
            irow = i;
        }
    }

    if (abs(point - array[irow][pcol])/point > tol) {
        cout << "Value " << point << " not within tolerance: " << endl; // throw error if the value is not close enough within tolerance
        cout << "Closest value is: " << array[irow][pcol] << endl;      // show closest value found regardless
        dm.print(array,nrows,ncols,true,"findvalue.txt");
        exit(0);
    } else {
        return irow;
    }
}

int data2 :: findvalueMCMC(double **array, int nrows, int ncols, double point, int pcol, double tol, int fullrow) {

    if (pcol >= ncols || pcol < 0) {                                  // Throw error if column is out of bounds
        cout << "Find value error: column is out of bounds" << endl;
        exit(0);
    }

    double error = abs(point-array[0][pcol]);                           // set initial error
    double newerror;                                                 // error tolerance of point
    int irow = 0;
    for (int i=0; i<nrows; ++i) {
        newerror = abs(point-array[i][pcol]);                           // calculate how new error
        if (newerror < error) {                                         // smallest error is point within tolerance
            error = newerror;                                           // set the old error as new error
            irow = i;
        }
    }

    if (abs(point - array[irow][pcol])/point > tol) {
        cout << "Value " << point << " not within tolerance: " << endl; // throw error if the value is not close enough within tolerance
        cout << "Closest value is: " << array[irow][pcol] << endl;      // show closest value found regardless
        dm.print(array,nrows,ncols,true,"findvalue.txt");
        dm.print(array,fullrow,ncols,true,"find.txt");
        exit(0);
    } else {
        return irow;
    }
}

void data2 :: zero(double** array, int nrows, int ncols) {
    for (int i=0; i<nrows; ++i) {
        for (int j=0; j<ncols; ++j) {
            array[i][j] = 0.0;
        }
    }
}



//*********************************************************************************
//*********************************************************************************

// Numerical Methods class
//*********************************************************************************
//*********************************************************************************

void nummeth :: pretovconv(double** &eos, int encol, int prcol, double enconv, double prconv, int nrows) {
    for (int i=0; i<nrows; ++i) {
        eos[i][encol] = eos[i][encol]*enconv;
        eos[i][prcol] = eos[i][prcol]*prconv;
    }
}

// dm/dr TOV equation (unitless)
double nummeth :: dmdr(double r, double p, double m, double en) {
    double eq;
    eq = r*r*en; // Dimensionless
    return eq;
}

// dp/dr TOV equation (unitless)
double nummeth :: dpdr(double r, double p, double m, double en) {
    double eq;
    eq = -0.5*(en+p)*(m+pow(r,3)*p)/(r*(r-m)); //Dimensionless
    return eq;
}

// dm/dr TOV equation (unitless)
double nummeth :: Ndmdr(double r, double en) {
    double eq;
    eq = pow(r,2)*en; // Dimensionless
    return eq;
}

// dp/dr TOV equation (unitless)
double nummeth :: Ndpdr(double r, double m, double en) {
    double eq;
    eq = -(en*m)/(2.0*pow(r,2)); //Dimensionless
    return eq;
}

// Newtonian love differential equation (unitless)
double nummeth :: Ndnldr(double r, double dpde, double en, int l, double nl) {
    double eq;
    eq = (l*(l+1.0))/r - 0.5*en*r/dpde - (nl+1.0)*nl/r;
    return eq;
}

double nummeth :: Rdnldr(double r, double m, double dpde, double en, double p, int l, double nl) {
    double eq, Fr, Qr;
    Fr = (r + 0.5*pow(r,3)*(p-en) )/(r-m);
    Qr = r/(r-m) * ( 0.5*(5.0*en + 9.0*p + (en+p)/dpde) - l*(l+1)/pow(r,2) ) - pow( (m+pow(r,3)*p)/(pow(r,2)-m*r) , 2);
    eq = -pow(nl,2)/r - Fr*nl/r - r*Qr;
    return eq;
}

double nummeth :: dIdr(double r, double m, double en, double p, double v, double wbar, double w, double I) {
    double eq, j, jp;
    jp = -(en+p)*exp(-v)*r;
    j = exp(-v)*(1.0-m/r);
    eq = 2.0/3.0*(pow(r,4.0)/j)*(en+p)*exp(-v)*wbar/w - 0.5*jp*I/j;
    return eq;
}

double nummeth :: dvdr(double r, double m, double p) {
    double eq;
    eq = (m+pow(r,3)*p)/(r*(r-m));
    return eq;
}

double nummeth :: dYdr(double r, double Y, double en, double p, double v, double m, double wbar) {
    double eq, jp, j;

    jp = -(en+p)*exp(-v)*r;
    j = exp(-v)*(1.0-m/r);
    eq = -2.0*jp*wbar/(r*j) -(0.5*jp/j + 4.0/r)*Y;
    return eq;
}

double nummeth :: dwdr(double r, double Y) {
    double eq;
    eq = Y;
    return eq;
}

// Single star icp tov solve
void nummeth :: tovsolve(double icp, double h, double** eos, int nrows, int ncols, int encol, int prcol, int dpdecol, double MR[2], bool print) {
    
    if (prcol >= ncols || prcol < 0 || encol >= ncols || encol < 0) {   // throw error for out of bounds columns
        cout << "TOV error: column is out of bounds" << endl;
        exit(0);
    }

    double r = 1e-20;       // set initial radius close to zero
    double p = icp;         // set central pressure
    double m = 0;           // inital mass to zero
    double en,dpde;
    double k1,k2,k3,k4;
    double l1,l2,l3,l4;
    double pmin = eos[0][prcol];

    if (print == true) {
        ofstream out("SingleICP.txt");
        while (p>pmin) {
            en = dm.interpolate(nrows, ncols, eos, p, prcol, encol, true);  // Interpolate the nuclear/quark eos
            dpde = dm.interpolate(nrows, ncols, eos, p, prcol, dpdecol, true);
            out << scientific << setprecision(15) << en << "  " << p << "  " << dpde << "  " << r << "  " << m << endl;
            // Integration routine (Possible to optimize stepsize)
            k1 = dmdr(r,p,m,en);
            l1 = dpdr(r,p,m,en);
            k2 = dmdr(r+h/2.0,p+l1*h/2.0,m+k1*h/2.0,en);
            l2 = dpdr(r+h/2.0,p+l1*h/2.0,m+k1*h/2.0,en);
            k3 = dmdr(r+h/2.0,p+l2*h/2.0,m+k2*h/2.0,en);
            l3 = dpdr(r+h/2.0,p+l2*h/2.0,m+k2*h/2.0,en);
            k4 = dmdr(r+h,p+l3*h,m+k3*h,en);
            l4 = dpdr(r+h,p+l3*h,m+k3*h,en);
            m = m + 1.0/6.0*h*(k1 + 2.0*k2 + 2.0*k3 + k4);
            p = p + 1.0/6.0*h*(l1 + 2.0*l2 + 2.0*l3 + l4);
            r = r + h;
        }
    } else {
        while (p>pmin) {
            en = dm.interpolate(nrows, ncols, eos, p, prcol, encol, true);  // Interpolate the nuclear/quark eos
            dpde = dm.interpolate(nrows, ncols, eos, p, prcol, dpdecol, true);
            // Integration routine (Possible to optimize stepsize)
            
            k1 = dmdr(r,p,m,en);
            l1 = dpdr(r,p,m,en);
            k2 = dmdr(r+h/2.0,p+l1*h/2.0,m+k1*h/2.0,en);
            l2 = dpdr(r+h/2.0,p+l1*h/2.0,m+k1*h/2.0,en);
            k3 = dmdr(r+h/2.0,p+l2*h/2.0,m+k2*h/2.0,en);
            l3 = dpdr(r+h/2.0,p+l2*h/2.0,m+k2*h/2.0,en);
            k4 = dmdr(r+h,p+l3*h,m+k3*h,en);
            l4 = dpdr(r+h,p+l3*h,m+k3*h,en);
            
            // for white dwarf
            /*
            k1 = Ndmdr(r,en);
            l1 = Ndpdr(r,m,en);
            k2 = Ndmdr(r+h/2.0,en);
            l2 = Ndpdr(r+h/2.0,m+k1*h/2.0,en);
            k3 = Ndmdr(r+h/2.0,en);
            l3 = Ndpdr(r+h/2.0,m+k2*h/2.0,en);
            k4 = Ndmdr(r+h,en);
            l4 = Ndpdr(r+h,m+k3*h,en);
            */
            m = m + 1.0/6.0*h*(k1 + 2.0*k2 + 2.0*k3 + k4);
            p = p + 1.0/6.0*h*(l1 + 2.0*l2 + 2.0*l3 + l4);
            r = r + h;
        }
    }
    
    MR[0] = m;   // save the final star mass
    MR[1] = r;   // save final star radius
    //cout << "icp: " << icp << "  M: " << m << "  R: " << conv.rnonetokm(r) << endl;
}

// Multiple central pressure star configurations
void nummeth :: multitov(double h, double** eos, int nrows, int ncols, int encol, int prcol, int dpdecol, int npoints, string filesave, double pr0) {
    double en, icp;
    int space;
    
    if (prcol >= ncols || prcol < 0 || encol >= ncols || encol < 0) { // throw error for out of bounds columns
        cout << "TOV error: column is out of bounds" << endl;
        exit(0);
    }

    if (npoints > nrows) {
        cout << " too many points: calculating each pr point" << endl;
        npoints = nrows;
    }
    
    double MR[2];                                            // initalize the mass radius array
    ofstream out(filesave + "_RM.txt");                         // File name save
    int index = dm.findvalue(eos,nrows,ncols,pr0,prcol,0.1);

    space = ceil(1.0*(nrows-(index+1.0))/npoints);

    while (index < nrows) {
        icp = eos[index][prcol];
        tovsolve(icp, h, eos, nrows, ncols, encol, prcol, dpdecol, MR, false);                 // solve tov equations
        en = dm.interpolate(nrows, ncols, eos, icp, prcol, encol, true);
        out << en*conv.energyCONV(1,0) << "  " << icp*conv.energyCONV(1,0) << "  " << conv.rnonetokm(MR[1]) << " " << MR[0] << endl;
        //cout << to_string(index) + "/" + to_string(nrows) << "  " << "mue: " << dm.interpolate(nrows, ncols, eos, icp, prcol, 0, true) << " icp: " << icp*conv.energyCONV(1,0) << " ice: " << en*conv.energyCONV(1,0) << " RM: " << conv.rnonetokm(MR[0][1]) << " " << MR[0][0] << endl;
        index = index + space;  //spacing of each mass radius point in pressure space
        //cout << "LOVE: " << Nlove("SingleICP.txt", 3, 2, 0, MR[0][1], 1e-3) << endl;
    }
}

void nummeth :: speedofsound(string EOS, int mubcol, int muecol, int nnbcol, int nnqcol, int nbcol, int xpcol, int Escol, int dEscol, int nencol, int encol, int nprcol, int prcol, int Xcol) {
    double** arr;
    double X, ce2, cs2, nnb, nen, npr, dEs, xp, Es, dxdn, mub, mue;
    int i;
    dm.importdata(EOS, arr);
    int nrows = dm.rowcount(EOS);
    ofstream out("EOS_SOS.txt");
    double cmev1 = conv.energyCONV(6, 0);
    
    // arr columns:  0 ,  1 , 2 , 3 , 4 ,  5  , 6 , 7 ,  8 ,  9 , 10 , 11 , 12 , 13
    // arr columns: mub, mue, nh, nq, nt, beta, xp, Es, dEs, nen, ten, npr, tpr, X
    nnb = arr[0][nnbcol];
    nen = arr[0][nencol];
    npr = arr[0][nprcol];
    dEs = arr[0][dEscol];
    Es = arr[0][Escol];
    xp = arr[0][xpcol];
    mub = arr[0][mubcol];
    mue = arr[0][muecol];
    X = arr[0][Xcol];
    
    ce2 = (arr[1][nprcol] - arr[0][nprcol])/(arr[1][nencol] - arr[0][nencol]);
    dxdn = (arr[1][xpcol] - arr[0][xpcol])/(arr[1][nnbcol] - arr[0][nnbcol]);
    cs2 = ce2 - nnb/(nen + npr) * ( pow(nnb,2.0)*(-4.0*dEs*(1.0-2.0*xp) + 4.0*Es*(1.0-2.0*xp)/(3.0*nnb)) ) * dxdn;
    out << mub << "  " << mue << "  " << nnb*cmev1 << "  " << nen*cmev1 << "  " << npr*cmev1 << "  " << X << "  " << ce2 << "  " << cs2 << endl;
    
    nnb = arr[1][nnbcol];
    nen = arr[1][nencol];
    npr = arr[1][nprcol];
    dEs = arr[1][dEscol];
    Es = arr[1][Escol];
    xp = arr[1][xpcol];
    mub = arr[1][mubcol];
    mue = arr[1][muecol];
    X = arr[1][Xcol];
    
    ce2 = (arr[2][nprcol] - arr[0][nprcol])/(arr[2][nencol] - arr[0][nencol]);
    dxdn = (arr[2][xpcol] - arr[0][xpcol])/(arr[2][nnbcol] - arr[0][nnbcol]);
    cs2 = ce2 - nnb/(nen + npr) * ( pow(nnb,2.0)*(-4.0*dEs*(1.0-2.0*xp) + 4.0*Es*(1.0-2.0*xp)/(3.0*nnb)) ) * dxdn;
    out << mub << "  " << mue << "  " << nnb*cmev1 << "  " << nen*cmev1 << "  " << npr*cmev1 << "  " << X << "  " << ce2 << "  " << cs2 << endl;
    
    // NUCLEAR PHASE
    i = 2;
    X = arr[i][Xcol];
    while (X < 0) {
        nnb = arr[i][nnbcol];
        nen = arr[i][nencol];
        npr = arr[i][nprcol];
        dEs = arr[i][dEscol];
        Es = arr[i][Escol];
        xp = arr[i][xpcol];
        mub = arr[i][mubcol];
        mue = arr[i][muecol];
        
        ce2 = ((arr[i+1][nprcol] - arr[i-1][nprcol])/(arr[i+1][nencol] - arr[i-1][nencol]) + (arr[i+2][nprcol] - arr[i-2][nprcol])/(arr[i+2][nencol] - arr[i-2][nencol]))*0.5;
        dxdn = ((arr[i+1][xpcol] - arr[i-1][xpcol])/(arr[i+1][nnbcol] - arr[i-1][nnbcol]) + (arr[i+2][xpcol] - arr[i-2][xpcol])/(arr[i+2][nnbcol] - arr[i-2][nnbcol]))*0.5;
        cs2 = ce2 - nnb/(nen + npr) * ( pow(nnb,2.0)*(-4.0*dEs*(1.0-2.0*xp) + 4.0*Es*(1.0-2.0*xp)/(3.0*nnb)) ) * dxdn;
        out << mub << "  " << mue << "  " << nnb*cmev1 << "  " << nen*cmev1 << "  " << npr*cmev1 << "  " << X << "  " << ce2 << "  " << cs2 << endl;
        i = i+1;
        X = arr[i][Xcol];
    }
    
    dm.cleanup(arr, nrows); // delete original array
}

void nummeth :: twopointD(double **oriarr, double ** &newarr, int ncols, int nrows, int ycol, int xcol) {
    double dydx1, dydx2;
    
    // create new array with an additonal column for the derivative
    newarr = new double*[nrows];
    for (int i = 0; i<nrows; i++) {
        newarr[i] = new double[ncols+1];
    }
    
    // fill the new array with the same data
    for (int j=0; j<ncols; ++j) {
        for (int i=0; i<nrows; ++i) {
            newarr[i][j] = oriarr[i][j];
        }
    }
    
    newarr[0][ncols] = (newarr[1][ycol] - newarr[0][ycol])/(newarr[1][xcol] - newarr[0][xcol]);
    newarr[1][ncols] = (newarr[2][ycol] - newarr[0][ycol])/(newarr[2][xcol] - newarr[0][xcol]);
    
    for (int i=2; i<(nrows-2); ++i) {
        dydx1 = (newarr[i+1][ycol] - newarr[i-1][ycol])/(newarr[i+1][xcol] - newarr[i-1][xcol]);
        dydx2 = (newarr[i+2][ycol] - newarr[i-2][ycol])/(newarr[i+2][xcol] - newarr[i-2][xcol]);
        newarr[i][ncols] = 0.5*(dydx1 + dydx2);
    }
    
    newarr[nrows-2][ncols] = (newarr[nrows-1][ycol] - newarr[nrows-3][ycol])/(newarr[nrows-1][xcol] - newarr[nrows-3][xcol]);
    newarr[nrows-1][ncols] = (newarr[nrows-1][ycol] - newarr[nrows-2][ycol])/(newarr[nrows-1][xcol] - newarr[nrows-2][xcol]);
}

// calculate Newtonian love number, inertia, and quadrupole moment
void nummeth :: Nlove(string txtfile, int dpdecol, int rcol, double encol, double mcol, double h, double ILQ[1][3]) {
    double** eos;
    dm.importdata(txtfile, eos);
    int nrows = dm.rowcount(txtfile);
    int ncols = dm.colcount(txtfile);

    double r0 = 5e-2;   // Set the integration starting point (close to zero)
    double en, cs2, m, pst, m2, X, PHI, e0, e2, w, r, a, b, c, d, X0;
    double Q = 0;
    double I0 = 0;
    double I2 = 0;
    int l = 2;          // Set the order to 2
    double nl = l;      // boundary condition for clairaut-radau equation
    en = dm.interpolate(nrows, ncols, eos, r0, rcol, encol, true);      // central energy density
    cs2  = dm.interpolate(nrows, ncols, eos, r0, rcol, dpdecol, true);  // central speed of sound
    double R = dm.findmax(eos,rcol,nrows,ncols);                        // max radius of unperterbed star
    m = dm.interpolate(nrows, ncols, eos, R, rcol, mcol, true);         // max mass of unperterbed star
    w = sqrt(m/pow(R,3))*0.77;      // rotational speed of star (0.77 to get to kepler frequency)
    pst = 1/6.0*pow(w*r0,2);        // initial condition for p*
    m2 = 1.0/30*en*pow(w,2)*pow(r0,5)/cs2;      // initial condition for m2
    double PHI0 = 0;        //Boundary condition for phi2'                               

    r = r0;             //integrate from center
    X = 0;              //try phi2(r0)=0
    PHI = PHI0;         //initialize phi2'
    while (r<R) {
        en = dm.interpolate(nrows, ncols, eos, r, rcol, encol, true);  // Interpolate the nuclear/quark eos
        cs2  = dm.interpolate(nrows, ncols, eos, r, rcol, dpdecol, true);   //speed of sound
        X = X + h*PHI;
        PHI = PHI + h*( 6/pow(r,2)*X - 0.5*en*X/cs2 -2.0*PHI/r - 1.0/6*en*pow(r,2)*pow(w,2)/cs2 );
        r = r+h;
    }
    a = X;      //define a(R)   X = a(R) + X(0)b(R)
    c = PHI;    //define c(R)   PHI = c(R) + X(0)d(R)

    r = r0;             //integrate again from center 
    X = 1;              //try phi2(r0)=1
    PHI = PHI0;         //phi2'=0
    while (r<R) {
        en = dm.interpolate(nrows, ncols, eos, r, rcol, encol, true);  // Interpolate the nuclear/quark eos
        cs2  = dm.interpolate(nrows, ncols, eos, r, rcol, dpdecol, true);   //speed of sound
        X = X + h*PHI;
        PHI = PHI + h*( 6/pow(r,2)*X - 0.5*en*X/cs2 -2.0*PHI/r - 1.0/6*en*pow(r,2)*pow(w,2)/cs2 );
        r = r+h;
    }  
    b = X - a;  //define b(R)
    d = PHI - c;    //define d(R)
    X0 = -(c+3.0*a/r)/(d+3.0*b/r);  // find the actual X(r0) from the boundary conditions

    double X1,X2,PHI1,PHI2,phi;
    r = r0;
    X1 = 0;         // solve the diff eqs again to get the actual solution X(r) = a(r) + X(0)b(r)
    X2 = 1;         
    PHI1 = PHI0;    // solve the diff eqs again to get the actual solution PHI(r) = c(r) + X(0)d(r)
    PHI2 = PHI0;
    while (r<R) {
        en = dm.interpolate(nrows, ncols, eos, r, rcol, encol, true);  // Interpolate the nuclear/quark eos
        cs2  = dm.interpolate(nrows, ncols, eos, r, rcol, dpdecol, true);   //speed of sound
        m = dm.interpolate(nrows, ncols, eos, r, rcol, mcol, true);         //unperterbed mass
        nl = nl + h*Ndnldr(r,cs2,en,l,nl);              // clairaut-radau equation
        pst = pst + h*(1/3.0*pow(w,2)*r - 0.5*m2/pow(r,2));     //p* equation
        m2 = m2+ h*pow(r,2)*en*pst/cs2;         // second order correction to mass

        X1 = X1 + h*PHI1;
        X2 = X2 + h*PHI2;
        PHI1 = PHI1 + h*( 6/pow(r,2)*X1 - 0.5*en*X1/cs2 -2.0*PHI1/r - 1.0/6*en*pow(r,2)*pow(w,2)/cs2 );
        PHI2 = PHI2 + h*( 6/pow(r,2)*X2 - 0.5*en*X2/cs2 -2.0*PHI2/r - 1.0/6*en*pow(r,2)*pow(w,2)/cs2 );
        a = X1;
        b = X2 - a;
        c = PHI1;
        d = PHI2 - c;
        phi = a + b*X0; //phi2 particular solution

        //U2 = 1.0/6*pow(r*w,2);      // potential due to rotation
        e0 = 2.0*pow(r,2)*pst/m;    // first order correction to r
        e2 = -pow(r,2)/m*(phi+1/3.0*pow(w,2)*pow(r,2));     // second order correction to r
        I0 = I0 + h*2.0/3.0*(en*pow(r,4));          // first order inertia
        //I2 = I2 - h*2.0/3.0*(pow(r,4)*(e0-0.2*e2)*dedr);        // second order correction inertia
        I2 = I2 - h*1.0/3.0*(pow(r,2)*(0.2*e2-e0)*cs2*m*en);        // second order correction inertia
        Q = phi*pow(R,3);       // Quadrupole moment
        r = r+h;
    }
    double k2 = 0.5*(l-nl)/(l+1+nl);        // love number from clairaut radau eq
    //cout << "I love Q: " << (I0+I2)/(m*pow(r,2)) << "  " << k2 << "  " << Q*1.989e33*pow(2.95e5,2)/pow(10,49) << endl;
    //cout << "I love Q: " << I0/pow(m,3) << "  " << k2 << "  " << Q*m/pow(w*I0,2) << endl;
    //cout << "Check: (Boundary Condition, Love NUmber) (" << phi + 3.0*2.0*k2*U2/r << " , " << 3.0/2.0*phi/pow(w*r,2) << ")" << endl;         // check if matching condition is satisfied and love number matches
    //cout << "w: " << w << " wk: " << sqrt((m+m2)/pow(R+e0-e2/2.0,3)) << endl;   //compare frequency to theoretical kepler frequency
    ILQ[0][0] = 4.0*I0/pow(m,3);
    ILQ[0][1] = k2*pow(2.0*r/m,5)/3.0;
    ILQ[0][2] = 2*Q*m/pow(w*I0,2);
    cout << "k2: "  << k2 << "  I: " << I0 << endl;
    dm.cleanup(eos, nrows);
}

// relativistic love number, inertia, and quadrupole moment
void nummeth :: Rlove(string txtfile, int dpdecol, int rcol, double encol, double prcol, double mcol, double h, double ILQR[4]) {
    double** eos;
    dm.importdata(txtfile, eos);
    int nrows = dm.rowcount(txtfile);
    int ncols = dm.colcount(txtfile);

    double r0 = 1e-5;   // Set the integration starting point (close to zero)
    double en, p, cs2, m, r, v, wbar,Y,vaR;
    double r1, r2, r3, r4, l1, l2, l3, l4, p1, p2, p3, p4, m1, m2, m3, m4;
    double I = 0;
    int l = 2;          // Set the order to 2
    double nl = l*1.0;      // boundary condition for clairaut-radau equation
    en = dm.interpolate(nrows, ncols, eos, r0, rcol, encol, true);      // central energy density
    cs2  = dm.interpolate(nrows, ncols, eos, r0, rcol, dpdecol, true);  // central speed of sound
    double R = dm.findmax(eos,rcol,nrows,ncols);                        // max radius of unperterbed star
    double M = dm.interpolate(nrows,ncols,eos,R,rcol,mcol,true);                 

    r = r0;             //integrate from center
    v = 0;
    m = 0;
    while (r<R) {
        m = dm.interpolate(nrows, ncols, eos, r, rcol, mcol, true);
        p = dm.interpolate(nrows, ncols, eos, r, rcol, prcol, true);
        m1 = dvdr(r,m,p);
        m2 = dvdr(r+h/2.0,m,p);
        m3 = dvdr(r+h/2.0,m,p);
        m4 = dvdr(r+h,m,p);
        v = v + 1.0/6.0*h*(m1 + 2.0*m2 + 2.0*m3 + m4);
        r = r+h;
    }
    vaR = v;

    cout << "r: " << r << "  m: " << m << endl;
    cout << "vaR: " << vaR << endl;

    r = r0;
    v = log(1.0-M/R) - vaR;
    wbar = r0;
    Y = r0;
    while (r<R) {
        en = dm.interpolate(nrows, ncols, eos, r, rcol, encol, true);  // Interpolate the nuclear/quark eos
        p = dm.interpolate(nrows, ncols, eos, r, rcol, prcol, true);
        cs2  = dm.interpolate(nrows, ncols, eos, r, rcol, dpdecol, true);   //speed of sound
        m = dm.interpolate(nrows, ncols, eos, r, rcol, mcol, true);         //unperterbed mass

        r1 = Rdnldr(r,m,cs2,en,p,l,nl);
        r2 = Rdnldr(r+h/2.0,m,cs2,en,p,l,nl+r1*h/2.0);
        r3 = Rdnldr(r+h/2.0,m,cs2,en,p,l,nl+r2*h/2.0);
        r4 = Rdnldr(r+h,m,cs2,en,p,l,nl+r3*h);
        l1 = dwdr(r,Y);
        p1 = dYdr(r,Y,en,p,v,m,wbar);
        l2 = dwdr(r+h/2.0,Y+p1*h/2.0);
        p2 = dYdr(r+h/2.0,Y+p1*h/2.0,en,p,v,m,wbar+l1*h/2.0);
        l3 = dwdr(r+h/2.0,Y+p2*h/2.0);
        p3 = dYdr(r+h/2.0,Y+p2*h/2.0,en,p,v,m,wbar+l2*h/2.0);
        l4 = dwdr(r+h,Y+p3*h);
        p4 = dYdr(r+h,Y,en,p+p3*h,v,m,wbar+l3*h);
        nl = nl + 1.0/6.0*h*(r1 + 2.0*r2 + 2.0*r3 + r4); // clairaut-radau equation
        wbar = wbar + 1.0/6.0*h*(l1 + 2.0*l2 + 2.0*l3 + l4); 
        Y = Y + 1.0/6.0*h*(p1 + 2.0*p2 + 2.0*p3 + p4);
        
        I = Y*pow(r,4)/(3.0*wbar + Y*r);
        m1 = dvdr(r,m,p);
        m2 = dvdr(r+h/2.0,m,p);
        m3 = dvdr(r+h/2.0,m,p);
        m4 = dvdr(r+h,m,p);
        v = v + 1.0/6.0*h*(m1 + 2.0*m2 + 2.0*m3 + m4);
        r = r+h;
    }
    cout << "check: " << exp(v) << "  " << (1.0-M/R) << endl;

    double C = M/(2.0*R);   // compactness
    cout << "Compactness: " << C << endl;
    double d1 = 2.0*C*(6.0-3.0*nl+15.0*C*nl-24.0*C);
    double d2 = 4.0*pow(C,3)*(13.0-11.0*nl+3.0*C*nl-2.0*C+2.0*pow(C,2)*(nl+1.0));
    double d3 = 3.0*pow(1.0-2.0*C,2)*(2.0-nl+2.0*C*nl-2.0*C)*log(1.0-2.0*C);
    double k2 = 8.0*pow(C,5)/5.0*pow(1.0-2.0*C,2)*(2.0-2.0*C+2.0*nl*C-nl)/(d1+d2+d3);
    ILQR[0] = 4.0*I/pow(m,3);         // DIMENSIONLESS I
    //ILQ[0][1] = k2*pow(2.0*r/m,5)/3.0;  // DIMENSIONLESS LOVE
    ILQR[1] = 2.0/3.0*k2*pow(C,-5);       // Tidal Deformability
    ILQR[2] = 0;
    dm.cleanup(eos, nrows);
    cout << "nl: " << nl << endl;
    cout << "k2: " << k2 << "  I: " << I << endl;
    cout << "R: " << conv.rnonetokm(R) << endl;
    cout << "M: " << M << endl;
    cout << "TDeform: " << 2.0/3.0*k2*pow(C,-5) << endl;
}

void nummeth :: RILOVEQ(double** eosdat2, int nrows, int ncols, int npoints, int dpdecol, int encol, int prcol, double pr0, double prf) {
    double MR[2];
    double icp;
    int index = dm.findvalue(eosdat2,nrows,ncols,pr0,prcol,0.1);
    int final = dm.findvalue(eosdat2,nrows,ncols,prf,prcol,0.1);
    double ILQR[4];
    ofstream out("RILoveQ.txt");
    
    if (npoints > (final-index)) {
        cout << "too many points!" << endl;
        exit(0);
    }
    int space = floor((final-index)/npoints);

    while (index < final) {
        icp = eosdat2[index][prcol];
        tovsolve(icp, 1e-6, eosdat2,nrows,ncols,encol,prcol,dpdecol,MR,true);   // solve tov, output is (en,pr,cs2,r,m)
        Rlove("SingleICP.txt", 2, 3, 0, 1, 4, 1e-3, ILQR);
        out << eosdat2[index][0] << "  " << icp*conv.energyCONV(1,0) << "  " << MR[0] << "  " << ILQR[0] << "  " <<ILQR[1] << "  " << ILQR[2] << endl;
        index = index + space;
    }
}

// Get specific radius for 
double nummeth :: MCMCtovpoint(double h, double** eos, int nrows, int ncols, int encol, int prcol, int dpdecol, int npoints, double pr0, int nMR, int nL, double** NICER, double** LIGO, double** &PRM) {
    double icp;
    int x;
    double MR[2];   // initalize the mass radius array
    int row[nMR+nL];

    if (npoints > nrows) {
        cout << " too many points: calculating each pr point" << endl;
        npoints = nrows;
    }
    
    // get the spacing for a given npoints and starting pressure (at most npoints)
    int index = dm.findvalue(eos,nrows,ncols,pr0,prcol,1.0);        // find initial pressure row and start here
    int space = ceil(1.0*(nrows-(index+1.0))/npoints);              // get spacing given npoints and starting pressure
    if (space < 2) {
        cout << "space to be sampled is to small for the number of samples requested: " << nrows-(index+1.0) << "/" << npoints << endl;
        exit(0);
    }
    
    int limit = floor((nrows-(index+1.0))/space);
     // create array to tabulate each icp, radius, mass
    double** prm;
    prm = new double*[limit];
    for (int i=0; i<limit; i++) {
        prm[i] = new double[3];
    }
    dm.zero(prm,limit,3);
    
    //#pragma omp parallel num_threads(6)
    //{
    //#pragma omp for private(x,icp,MR)
    for (int j=0; j<limit; ++j) {
        x = index + space*j;
        icp = eos[x][prcol];                // central pressure to integrate
        MCMCtovsolve(icp,h,eos,nrows,ncols,encol,prcol,dpdecol,MR);
        prm[j][0] = icp;                    // store the icp
        prm[j][1] = MR[1];  // store the radii
        prm[j][2] = MR[0];                  // store the masses
        
        // remove if running parallel (stops at maxmass)
        if (j>0) {
            if (prm[j][2] < prm[j-1][2]) {
                break;
            }
        }
    }
    //}

    double maxmass = dm.findmax(prm,2,limit,3);
    int maxrow = dm.findvalue(prm,limit,3,maxmass,2,0.001);
    PRM[nMR+nL][0] = prm[maxrow-1][0];     // max mass icp
    PRM[nMR+nL][1] = prm[maxrow-1][1];     // max mass radius
    PRM[nMR+nL][2] = maxmass;         // max mass

    // find given masses
    for (int i=0; i<(nMR+nL); ++i) { 
        if (PRM[i][2] >= prm[0][2] && PRM[i][2] <= prm[maxrow][2]) {  // check if mass is bounded
            row[i] = dm.findvalueMCMC(prm,maxrow+1,3,PRM[i][2],2,0.1,limit);  // if it is then find the mass searching from row 0 to maxmassrow
            if (prm[row[i]][2] < PRM[i][2]) {                     // make sure to save the row right above the found mass (for interpolation)
                row[i] = row[i] + 1;                            // make sure to save the row right above the found mass (for interpolation)
            } 
        } else if (PRM[i][2] < prm[0][2]){              // check if mass is too low
            cout << "mass too low" << endl;
            exit(0);
        } else if (PRM[i][2] > prm[maxrow][2])          // check if mass is too big
            row[i] = -2;
    }

    double** spline1;
    double** spline2;
    dm.cubicspline(prm,spline1,maxrow+1,2,0);
    dm.cubicspline(prm,spline2,maxrow+1,2,1);
    for (int i=0; i<(nMR+nL); ++i) {
        if (row[i] == -2) {     // if the mass was too large then mark the icp 0, <radius> calc later
            PRM[i][0] = 0.0;
        } else {                // if the mass was found then interpolate to get the icp,r
            PRM[i][0] = dm.splinecalc(spline1,row[i]-1,row[i],PRM[i][2]);
            PRM[i][1] = dm.splinecalc(spline2,row[i]-1,row[i],PRM[i][2]);
        }
    }

    // check to make sure spline doesn't mess up
    for (int i=0; i<(nMR+nL); ++i) {
        if (row[i] > 0) {
            if (prm[row[i]][1]>prm[row[i]-1][1]) { //if increasing
                if (PRM[i][1]<prm[row[i]-1][1] || PRM[i][1]>prm[row[i]][1]) {
                    cout << prm[row[i]-1][1] << "  " << PRM[i][1] << "  " << prm[row[i]][1] << endl;
                    cout << "correction to spline made A: " << PRM[i][1];
                    PRM[i][0] = dm.interpolate(maxrow+1,3,prm,PRM[i][2],2,0,true);
                    PRM[i][1] = dm.interpolate(maxrow+1,3,prm,PRM[i][2],2,1,true);
                    cout << " to " << PRM[i][1] << endl;
                }
            } else {
                if (PRM[i][1]>prm[row[i]-1][1] || PRM[i][1]<prm[row[i]][1]) {
                    cout << prm[row[i]][1] << "  " << PRM[i][1] << "  " << prm[row[i]-1][1] << endl;
                    cout << "correction to spline made B: " <<  PRM[i][1];
                    PRM[i][0] = dm.interpolate(maxrow+1,3,prm,PRM[i][2],2,0,true);
                    PRM[i][1] = dm.interpolate(maxrow+1,3,prm,PRM[i][2],2,1,true);
                    cout << " to " << PRM[i][1] << endl;
                }
            }
        }
    }
    
    // ##################### LINE INTEGRAL FORMULATION FOR LIKELIHOOD ##########################
    double t = 0;
    double dt = 0.01;
    double lklhood = 1.0;
    double R, M, R0, R1, M0, M1, ds, lint, rexp, dM, norm;
    for (int j=0; j<nMR; ++j) {
        lint = 0;
        rexp = 0;
        for (int i=0; i<(maxrow); ++i) {
            R0 = prm[i][1]; R1 = prm[i+1][1];
            M0 = prm[i][2]; M1 = prm[i+1][2];
            ds = sqrt(pow(R1-R0,2.0)+pow(M1-M0,2.0))*dt;
            dM = (M1-M0)*dt;
        
            t = 0;
            R = (1.0-t)*R0 + t*R1;
            M = (1.0-t)*M0 + t*M1;

            while(t<1.0) {
                lint = lint + ds*exp(-0.5*pow((conv.rnonetokm(R)-NICER[j][2])/NICER[j][3],2.0))*exp(-0.5*pow((M-NICER[j][0])/NICER[j][1],2.0));
                rexp = rexp + dM*R*exp(-0.5*pow((M-NICER[j][0])/NICER[j][1],2.0));
                t = t + dt;
                R = (1.0-t)*R0 + t*R1;
                M = (1.0-t)*M0 + t*M1;
            }
        }
        norm = 2.0*pow(erf((NICER[j][0]-prm[0][2])/(sqrt(2.0)*NICER[j][1]))-erf((NICER[j][0]-prm[maxrow][2])/(sqrt(2.0)*NICER[j][1])),-1.0);
        PRM[j][1] = rexp/(sqrt(2.0*pi)*NICER[j][1])*norm;
        lklhood = lklhood*lint;
    }

    dm.cleanup(prm,limit);          
    dm.cleanup(spline1,maxrow+1);
    dm.cleanup(spline2,maxrow+1);

    return lklhood;
}


// Single star icp tov solve
void nummeth :: MCMCtovsolve(double icp, double h, double** eos, int nrows, int ncols, int encol, int prcol, int dpdecol, double MR[2]) {
    double r = 1e-20;       // set initial radius close to zero
    double p = icp;         // set central pressure
    double m = 0;           // inital mass to zero
    double en;
    double k1,k2,k3,k4;
    double l1,l2,l3,l4;
    double pmin = eos[0][prcol];
    while (p>pmin) {
        en = dm.interpolate(nrows, ncols, eos, p, prcol, encol, true);  // Interpolate the nuclear/quark eos
        // Integration routine (Possible to optimize stepsize)
        k1 = dmdr(r,p,m,en);
        l1 = dpdr(r,p,m,en);
        k2 = dmdr(r+h/2.0,p+l1*h/2.0,m+k1*h/2.0,en);
        l2 = dpdr(r+h/2.0,p+l1*h/2.0,m+k1*h/2.0,en);
        k3 = dmdr(r+h/2.0,p+l2*h/2.0,m+k2*h/2.0,en);
        l3 = dpdr(r+h/2.0,p+l2*h/2.0,m+k2*h/2.0,en);
        k4 = dmdr(r+h,p+l3*h,m+k3*h,en);
        l4 = dpdr(r+h,p+l3*h,m+k3*h,en);
        m = m + 1.0/6.0*h*(k1 + 2.0*k2 + 2.0*k3 + k4);
        p = p + 1.0/6.0*h*(l1 + 2.0*l2 + 2.0*l3 + l4);
        r = r + h;
    }
    
    MR[0] = m;   // save the final star mass
    MR[1] = r;   // save final star radius
    //cout << "icp: " << icp << "  M: " << m << "  R: " << conv.rnonetokm(r) << endl;
}


// relativistic love number, inertia, and quadrupole moment
void nummeth :: RloveMCMC(double** EOS, int dpdecol, int encol, int prcol, double h, double ILQR[4], double icp, int nrows, int ncols) {
    double r0 = 1e-5;   // Set the integration starting point (close to zero)
    double en, cs2, r, v, wbar,Y,vaR;
    double I = 0;
    int l = 2;          // Set the order to 2
    double nl = l*1.0;      // boundary condition for clairaut-radau equation       
    double M, R;
    double r1,r2,r3,r4;
    double p1,p2,p3,p4;
    double m1,m2,m3,m4;
    double b1,b2,b3,b4;
    double a1,a2,a3,a4;
    double l1,l2,l3,l4;
    double pmin = EOS[0][prcol];
    double p = icp;         // set central pressure
    double m = 0;           // inital mass to zero
    
    r = r0;             //integrate from center
    v = 0;
    while (p>pmin) {
        en = dm.interpolate(nrows, ncols, EOS, p, prcol, encol, true);  // Interpolate the nuclear/quark eos
        a1 = dmdr(r,p,m,en);
        l1 = dpdr(r,p,m,en);
        m1 = dvdr(r,m,p);
        a2 = dmdr(r+h/2.0,p+l1*h/2.0,m+a1*h/2.0,en);
        l2 = dpdr(r+h/2.0,p+l1*h/2.0,m+a1*h/2.0,en);
        m2 = dvdr(r+h/2.0,m+a1*h/2.0,p+l1*h/2.0);
        a3 = dmdr(r+h/2.0,p+l2*h/2.0,m+a2*h/2.0,en);
        l3 = dpdr(r+h/2.0,p+l2*h/2.0,m+a2*h/2.0,en);
        m3 = dvdr(r+h/2.0,m+a2*h/2.0,p+l2*h/2.0);
        a4 = dmdr(r+h,p+l3*h,m+a3*h,en);
        l4 = dpdr(r+h,p+l3*h,m+a3*h,en);
        m4 = dvdr(r+h,m+a3*h,p+l3*h);
        m = m + 1.0/6.0*h*(a1 + 2.0*a2 + 2.0*a3 + a4);
        p = p + 1.0/6.0*h*(l1 + 2.0*l2 + 2.0*l3 + l4);
        v = v + 1.0/6.0*h*(m1 + 2.0*m2 + 2.0*m3 + m4);
        r = r+h;
    }
    M = m;
    R = r;
    vaR = v;

    r = r0;
    v = log(1.0-M/R) - vaR;
    wbar = r0;
    Y = r0;
    p = icp;
    m = 0;
    while (p>pmin) {
        en = dm.interpolate(nrows, ncols, EOS, p, prcol, encol, true);  // Interpolate the nuclear/quark eos
        cs2  = dm.interpolate(nrows, ncols, EOS, p, prcol, dpdecol, true);   //speed of sound
        
        a1 = dmdr(r,p,m,en);
        b1 = dpdr(r,p,m,en);
        r1 = Rdnldr(r,m,cs2,en,p,l,nl);
        l1 = dwdr(r,Y);
        p1 = dYdr(r,Y,en,p,v,m,wbar);
        m1 = dvdr(r,m,p);

        a2 = dmdr(r+h/2.0, p+b1*h/2.0, m+a1*h/2.0, en);
        b2 = dpdr(r+h/2.0, p+b1*h/2.0, m+a1*h/2.0, en);
        r2 = Rdnldr(r+h/2.0, m+a1*h/2.0, cs2, en, p+b1*h/2.0, l, nl+r1*h/2.0);
        l2 = dwdr(r+h/2.0, Y+p1*h/2.0);
        p2 = dYdr(r+h/2.0, Y+p1*h/2.0,en, p+b1*h/2.0, v+m1*h/2.0, m+a1*h/2.0, wbar+l1*h/2.0);
        m2 = dvdr(r+h/2.0, m+a1*h/2.0, p+b1*h/2.0);

        a3 = dmdr(r+h/2.0, p+b2*h/2.0, m+a2*h/2.0,en);
        b3 = dpdr(r+h/2.0, p+b2*h/2.0, m+a2*h/2.0,en);
        r3 = Rdnldr(r+h/2.0, m+a2*h/2.0, cs2, en, p+b2*h/2.0, l, nl+r2*h/2.0);
        l3 = dwdr(r+h/2.0, Y+p2*h/2.0);
        p3 = dYdr(r+h/2.0, Y+p2*h/2.0,en, p+b2*h/2.0, v+m2*h/2.0, m+a2*h/2.0, wbar+l2*h/2.0);
        m3 = dvdr(r+h/2.0, m+a2*h/2.0, p+b2*h/2.0);

        a4 = dmdr(r+h, p+b3*h, m+a3*h,en);
        b4 = dpdr(r+h, p+b3*h, m+a3*h,en);
        r4 = Rdnldr(r+h, m+a3*h, cs2,en, p+b3*h, l, nl+r3*h);
        l4 = dwdr(r+h, Y+p3*h);
        p4 = dYdr(r+h, Y+p3*h,en, p+b3*h, v+m3*h, m+a3*h, wbar+l3*h);
        m4 = dvdr(r+h, m+a3*h, p+b3*h);

        m = m + 1.0/6.0*h*(a1 + 2.0*a2 + 2.0*a3 + a4); 
        p = p + 1.0/6.0*h*(b1 + 2.0*b2 + 2.0*b3 + b4); 
        nl = nl + 1.0/6.0*h*(r1 + 2.0*r2 + 2.0*r3 + r4); // clairaut-radau equation
        wbar = wbar + 1.0/6.0*h*(l1 + 2.0*l2 + 2.0*l3 + l4); 
        Y = Y + 1.0/6.0*h*(p1 + 2.0*p2 + 2.0*p3 + p4);
        v = v + 1.0/6.0*h*(m1 + 2.0*m2 + 2.0*m3 + m4);
        r = r+h;
    }
    I = Y*pow(r,4)/(3.0*wbar + Y*r);

    double C = M/(2.0*R);   // compactness
    double d1 = 2.0*C*(6.0-3.0*nl+15.0*C*nl-24.0*C);
    double d2 = 4.0*pow(C,3)*(13.0-11.0*nl+3.0*C*nl-2.0*C+2.0*pow(C,2)*(nl+1.0));
    double d3 = 3.0*pow(1.0-2.0*C,2)*(2.0-nl+2.0*C*nl-2.0*C)*log(1.0-2.0*C);
    double k2 = 8.0*pow(C,5)/5.0*pow(1.0-2.0*C,2)*(2.0-2.0*C+2.0*nl*C-nl)/(d1+d2+d3);
    ILQR[0] = 4.0*I/pow(m,3);         // DIMENSIONLESS I
    ILQR[1] = 2.0/3.0*k2*pow(C,-5);       // Tidal Deformability
    ILQR[2] = 0;
    ILQR[3] = R;
    //cout << "nl: " << nl << endl;
    cout << "mass: " << M << endl;
    //cout << "Compactness: " << C << endl;
    cout << "k2: " << k2 << "  I: " << I << endl;
    cout << "TDeform: " << 2.0/3.0*k2*pow(C,-5) << endl;
    cout << "Radius: " << R*2.95324203 << endl;
}

// Find intersection point of 2 arrays
double nummeth :: findintersect(double **arr1, int nrows1, int ncols1, int xcol1, int ycol1, double **arr2, int nrows2, int ncols2, int xcol2, int ycol2) {
    double cand;
    double m1, b1, m2, b2;

    if (xcol1 >= ncols1 || xcol1 < 0) {                                     // Throw error if columns are out of bounds
        cout << "Find Intersect error: column is out of bounds" << endl;
        exit(0);
    }

    if (xcol2 >= ncols2 || xcol2 < 0) {                                     // Throw error if columns are out of bounds
        cout << "Find Intersect error: column is out of bounds" << endl;
        exit(0);
    }

    // Running linear spline interpolations
    for (int i=0; i<(nrows1-1); ++i) {
        for (int j=0; j<(nrows2-1); ++j) {
            m1 = (arr1[i+1][ycol1] - arr1[i][ycol1])/(arr1[i+1][xcol1] - arr1[i][xcol1]);
            m2 = (arr2[j+1][ycol2] - arr2[j][ycol2])/(arr2[j+1][xcol2] - arr2[j][xcol2]);
            b1 = arr1[i+1][ycol1] - m1*arr1[i+1][xcol1];
            b2 = arr2[j+1][ycol2] - m2*arr2[j+1][xcol2];
            cand = (b2 - b1)/(m1 - m2); // candidate intersection point of two lines
            // check if the candidate point is in the spline range
            if (cand >= arr1[i][xcol1] && cand <= arr1[i+1][xcol1] && cand >= arr2[j][xcol2] && cand <= arr2[j+1][xcol2]) {
                return cand;
            }
        }
    }
    cout << "No intersection point found " << endl;
    exit(0);
}

// Get specific radius for 
int nummeth :: MCMCcorr(double h, double** eos, int nrows, int ncols, int encol, int prcol, int dpdecol, int npoints, double pr0, double** &prm2) {
    double icp;
    int x;
    double MR[2];   // initalize the mass radius array

    if (npoints > nrows) {
        cout << " too many points: calculating each pr point" << endl;
        npoints = nrows;
    }
    
    // get the spacing for a given npoints and starting pressure (at most npoints)
    int index = dm.findvalue(eos,nrows,ncols,pr0,prcol,1.0);        // find initial pressure row and start here
    int space = ceil(1.0*(nrows-(index+1.0))/npoints);              // get spacing given npoints and starting pressure
    if (space < 2) {
        cout << "space to be sampled is to small for the number of samples requested: " << nrows-(index+1.0) << "/" << npoints << endl;
        exit(0);
    }
    
    int limit = floor((nrows-(index+1.0))/space);
     // create array to tabulate each icp, radius, mass
    double** prm;
    prm = new double*[limit];
    for (int i=0; i<limit; i++) {
        prm[i] = new double[3];
    }
    dm.zero(prm,limit,3);
    
    //#pragma omp parallel num_threads(6)
    //{
    //#pragma omp for private(x,icp,MR)
    int count = 0;
    for (int j=0; j<limit; ++j) {
        x = index + space*j;
        icp = eos[x][prcol];                // central pressure to integrate
        MCMCtovsolve(icp,h,eos,nrows,ncols,encol,prcol,dpdecol,MR);
        prm[j][0] = icp;                    // store the icp
        prm[j][1] = conv.rnonetokm(MR[1]);  // store the radii
        prm[j][2] = MR[0];                  // store the masses
        count = count+1;
        

        // remove if running parallel (stops at maxmass)
        if (j>0) {
            if (prm[j][2] < prm[j-1][2]) {
                break;
            }
        }
    }

    prm2 = new double*[count];
    for (int i=0; i<count; i++) {
        prm2[i] = new double[3];
    }

    for (int i=0; i<count; ++i) {
        for (int j=0; j<3; ++j) {
            prm2[i][j] = prm[i][j];
        }
    }

    dm.cleanup(prm,limit);

    //}
    return count;
}
