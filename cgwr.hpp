#ifndef CGWR_HPP
#define CGWR_HPP

#include <iostream>
#include <cmath>
#include <algorithm>

//#include <boost/rational.hpp>
//#include <boost/math/special_functions/factorials.hpp>
//#include <boost/math/special_functions/powm1.hpp>

/*
- CG:
< j1 j2 m1 m2 | J M >

- Wigner 3j or 6j:
(j1,j2,j3) or {j1,j2,j3}
(l1,l2,l3)    {m1,m2,m3}

- Wigner 9j:
{j1 j2 j3}
{j4 j5 j6}
{j7 j8 j9}

-Racah:
W(j1 j2 l2 l1 | j3 l3)

CG translated as array=[j1,j2,m1,m2,J,M] with _CGWR::C_CG token
or
Wigner 3j translated as array=[j1,j2,j3,m1,m2,m3] with _CGWR::C_W3j token
or
Wigner 6j translated as array=[j1,j2,j3,l1,l2,l3] with _CGWR::C_W6j token
or
and so on with usual notations

*/

class _CGWR { // Clebsch-Gordan Wigner Racah
public:
    // Token
    static const unsigned int C_CG;
    static const unsigned int C_W3j;
    static const unsigned int C_W6j;
    static const unsigned int C_W9j;
    static const unsigned int C_RACAH;

    // Init the class with the coeff and the token
    _CGWR(double long [], unsigned int);

    _CGWR& operator= ( const _CGWR& other ); //TODO

    double long CG(void);
    double long W3j(void);
    double long W6j(void);
    double long W9j(void); //TODO
    double long Racah(void);

private:
    double long j1, j2, j3, j;
    double long l1, l2, l3, l;
    double long m1, m2 , m3, m;

    unsigned int selected_Symb;

    // mid computations
    double long fact(double long ); // factorial
    double long Delta(double long , double long , double long ); // triangle relation
    double long delta_k(double long , double long ); // kronecker
    bool is_integer(double long ); // return true if interger
    bool is_halfinteger(double long ); // return true if half interger

};

#endif // CGWR_HPP

// The Theory of Atomic Structure and Spectra - Robert D. Cowan (1981)
