#include <iostream>
#include <cmath>
#include <algorithm> // min max
#include <iomanip> // precision

#include "cgwr.hpp"

int main(int argc, char **argv) {

    double long  j1, j2, j3, j;
    double long  l1, l2, l3;
    double long  m1, m2 , m3, m;

// CG
/*******************************************************/
    j1=0.5; j2=0.5;
    m1=0.5; m2=-0.5;
    j=1; m=0;

    double long  CG[]={j1,j2,m1,m2,j,m};
    _CGWR op_CG(CG, _CGWR::C_CG);
    std::cout << std::setprecision(20) << "CG: " << op_CG.CG() << "\n";
/*******************************************************/


// Wigner_3j
/*******************************************************/
    j1=0.5; j2=1; j3=1.5;
    m1=-0.5; m2=1; m3=-0.5;

    double long  W3j[]={j1,j2,j3,m1,m2,m3};
    _CGWR op_W3j(W3j, _CGWR::C_W3j);
    std::cout << std::setprecision(20) << "Wigner_3j: " << op_W3j.W3j() << "\n";
/*******************************************************/


// Wigner_6j
/*******************************************************/
//    j1=2.5; j2=1; j3=1.5;
//    l1=2; l2=0.5; l3=3;

    j1=1; j2=2; j3=2;
    l1=2; l2=3; l3=2;

    double long  W6j[]={j1,j2,j3,l1,l2,l3};
    _CGWR op_W6j(W6j, _CGWR::C_W6j);
    std::cout << std::setprecision(20) << "Wigner_6j: " << op_W6j.W6j() << "\n";
/*******************************************************/


    return EXIT_SUCCESS;
}
