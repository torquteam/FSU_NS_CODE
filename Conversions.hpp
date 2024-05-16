//
//  Conversions.h
//  NS Code
//
//  Created by Marc Salinas on 1/28/19.
//  Copyright © 2019 Marc Salinas. All rights reserved.
//

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

class Convert {
public:
    void convert();
    void printit();
    double energyCONV(int i, int j);
    double rnonetokm(double none);
};
