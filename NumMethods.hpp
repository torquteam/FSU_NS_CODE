//
//  NumMethods.hpp
//  NS Code
//
//  Created by Marc Salinas on 1/29/19.
//  Copyright © 2019 Marc Salinas. All rights reserved.
//

#ifndef NumMethods_hpp
#define NumMethods_hpp

#include <stdio.h>
#include <string>
using namespace std;

class data2 {
public:
    int rowcount(string txtfile);
    string filetype(string txtfile);
    int colcount(string txtfile);
    void importdata(string txtfile, double **&array);
    void print(double **&array, int numrows, int numcols, bool file, string filename);
    double interpolate(int numrows, int numcols, double** array, double point, int pointcol, int ycol, bool sortflag);
    void cubicspline(double**arr, double**&spline, int nrows, int xcol, int ycol);
    double splinecalc(double** spline, int x1, int x2, double x);
    bool sortasc(double ** array, int col, int numrows, int numcols);
    void cleanup(double** &array, int nrows);
    double findmin(double** array, int col, int numrows, int numcols);
    double findmax(double** array, int col, int numrows, int numcols);
    int findvalue(double **array, int nrows, int ncols, double point, int pcol, double tol);
    int findvalueMCMC(double **array, int nrows, int ncols, double point, int pcol, double tol, int fullrow);
    void zero(double** array, int nrows, int ncols);
};

class nummeth {
public:
    void pretovconv(double** &eos, int encol, int prcol, double enconv, double prconv, int nrows);
    void tovsolve(double icp, double h, double** eos, int nrows, int ncols, int encol, int prcol, int dpdecol, double MR[2], bool print);
    void multitov(double h, double** eos, int nrows, int ncols, int encol, int prcol, int dpdecol, int npoints, string filesave, double pr0);
    void speedofsound(string EOS, int mubcol, int muecol, int nnbcol, int nnqcol, int nbcol, int xpcol, int Escol, int dEscol, int nencol, int encol, int nprcol, int prcol, int Xcol);
    void twopointD(double **oriarr, double ** &newarr, int ncols, int nrows, int ycol, int xcol);
    void Nlove(string txtfile, int dpdecol, int rcol, double encol, double mcol, double h, double ILQ[1][3]);
    void Rlove(string txtfile, int dpdecol, int rcol, double encol, double prcol, double mcol, double h, double ILQR[4]);
    void RILOVEQ(double** eosdat2, int nrows, int ncols, int npoints, int dpdecol, int encol, int prcol, double pr0, double prf);
    double MCMCtovpoint(double h, double** eos, int nrows, int ncols, int encol, int prcol, int dpdecol, int npoints, double pr0, int nMR, int nL, double** NICER, double** LIGO, double** &PRM);
    void MCMCtovsolve(double icp, double h, double** eos, int nrows, int ncols, int encol, int prcol, int dpdecol, double MR[2]);
    void RloveMCMC(double** EOS, int dpdecol, int encol, int prcol, double h, double ILQR[4], double icp, int nrows, int ncols);
    double findintersect(double **arr1, int nrows1, int ncols1, int xcol1, int ycol1, double **arr2, int nrows2, int ncols2, int xcol2, int ycol2);
    int MCMCcorr(double h, double** eos, int nrows, int ncols, int encol, int prcol, int dpdecol, int npoints, double pr0, double** &prm);
private:
    double dmdr(double r, double p, double m, double en);
    double dpdr(double r, double p, double m, double en);
    double Ndmdr(double r, double en);
    double Ndpdr(double r, double m, double en);
    double Ndnldr(double r, double dpde, double en, int l, double nl);
    double Rdnldr(double r, double m, double dpde, double en, double p, int l, double nl);
    double dIdr(double r, double m, double en, double p, double v, double wbar, double w, double I);
    double dvdr(double r, double m, double p);
    double dYdr(double r, double Y, double en, double p, double v, double m, double wbar);
    double dwdr(double r, double Y);
};

#endif /* NumMethods_h */
