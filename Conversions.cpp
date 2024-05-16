//
//  Conversions.cpp
//  NS Code
//
//  Created by Marc Salinas on 1/28/19.
//  Copyright © 2019 Marc Salinas. All rights reserved.
//

#include "Conversions.hpp"
#include <iostream>
#include <cmath>

using namespace std;

/*
      (to)|  Mevfm3 | none | gcc | J/m3 | 1/km2 | dyne/cm2 | MeV4 |
 (from)   |                                                       |
 ---------|-------------------------------------------------------|
 MeVfm3   |   0,0     0,1    0,2   0,3     0,4      0,5       0,6 |
 none     |   1,0     1,1    1,2   1,3     1,4      1,5       1,6 |
 gcc      |   2,0     2,1    2,2   2,3     2,4      2,5       2,6 |
 J/m3     |   3,0     3,1    3,2   3,3     3,4      3,5       3,6 |
 1/km2    |   4,0     4,1    4,2   4,3     4,4      4,5       4,6 |
 dyne/cm2 |   5,0     5,1    5,2   5,3     5,4      5,5       5,6 |
 MeV4     |   6,0     6,1    6,2   6,3     6,4      6,5       6,6 |
 -----------------------------------------------------------------|
 
*/

static const int numtype = 7;                   // Number of units available
static double enCONV[numtype][numtype];

void Convert :: convert() {
    
    // Set up identity conversions
    for (int i=0; i<numtype; ++i) {
        enCONV[i][i] = 1.0;
    }
    
    enCONV[0][1] = 0.000290173;                 // MeVfm3 to None
    enCONV[0][2] = 1.7827 * pow(10,12);         // MeVfm3 to gcc
    enCONV[0][5] = 1.6021766*pow(10,33);        // MeVfm3 to dyne/cm2
    enCONV[0][6] = pow(197.32698,3);               // MeVfm3 to Mev4
    enCONV[1][2] = 6.143577506*pow(10,15);      // unitless to gcc
    enCONV[1][3] = 5.52145404 * pow(10,35);     // unitless to J/m3
    enCONV[1][4] = 0.00456073;                  // unitless to 1/km2
    enCONV[2][3] = 8.98735962*pow(10,19);       // gcc to J/m3
    enCONV[3][4] = 8.26001326 * pow(10,-39);    // J/m3 to 1/km2
    
    // Transpose the reciprocal
    for (int i=1; i<numtype; ++i) {
        for (int j=0; j<i; ++j) {
            enCONV[i][j] = 1.0/enCONV[j][i];
        }
    }
    
    // Cyclic permutations to fill missing conversions
    enCONV[0][3] = enCONV[0][1]*enCONV[1][3];
    enCONV[0][4] = enCONV[0][1]*enCONV[1][4];
    enCONV[1][5] = enCONV[1][0]*enCONV[0][5];
    enCONV[1][6] = enCONV[1][0]*enCONV[0][6];
    enCONV[2][4] = enCONV[2][1]*enCONV[1][4];
    enCONV[2][5] = enCONV[2][0]*enCONV[0][5];
    enCONV[2][6] = enCONV[2][0]*enCONV[0][6];
    enCONV[3][5] = enCONV[3][1]*enCONV[1][0]*enCONV[0][5];
    enCONV[3][6] = enCONV[3][1]*enCONV[1][0]*enCONV[0][6];
    enCONV[4][5] = enCONV[4][1]*enCONV[1][0]*enCONV[0][6];
    enCONV[4][6] = enCONV[4][1]*enCONV[1][0]*enCONV[0][6];
    enCONV[5][6] = enCONV[5][0]*enCONV[0][6];
    
    // Transpose the reciprocal
    for (int i=1; i<numtype; ++i) {
        for (int j=0; j<i; ++j) {
            enCONV[i][j] = 1/enCONV[j][i];
        }
    }
}

double Convert :: energyCONV(int i, int j) {
    convert();
    return enCONV[i][j];
}

// ** unitless r to km **
double Convert :: rnonetokm(double none) {
    double conv = none * 2.95324203;
    return conv;
}


