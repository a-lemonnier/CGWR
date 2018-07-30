#include "cgwr.hpp"

// Token definitions
/***************************************/
const unsigned int _CGWR::C_CG    =   1;
const unsigned int _CGWR::C_W3j   =   2;
const unsigned int _CGWR::C_W6j   =   3;
const unsigned int _CGWR::C_W9j   =   4;
const unsigned int _CGWR::C_RACAH =   5;
/***********/


// Copy arguments
/***************************************/
_CGWR::_CGWR(double long M[], unsigned int Symb) {
   
    if (Symb==C_CG) {
        j1=M[0]; j2=M[1];
        m1=M[2]; m2=M[3];
        j=M[4];  m=M[5];
    }
    if (Symb==C_W3j) {
        j1=M[0]; j2=M[1]; j3=M[2];
        m1=M[3]; m2=M[4]; m3=M[5];
    }
    if (Symb==C_W6j) {
        j1=M[0]; j2=M[1]; j3=M[2];
        l1=M[3]; l2=M[4]; l3=M[5];
    }
    if (Symb==C_W9j) {
        //TODO
    }
    if (Symb==C_RACAH) {
        j1=M[0]; j2=M[1]; j3=M[4];
        l1=M[2]; l2=M[3]; l3=M[5];
    }
    if (Symb!=C_CG &&
        Symb!=C_W3j&&
        Symb!=C_W6j&&
        Symb!=C_W9j&&
        Symb!=C_W9j&&
        Symb!=C_RACAH) {
            std::cerr << "/!\\ invalid symbol.\n";
            exit(EXIT_FAILURE);
    }
    selected_Symb=Symb;
}
/***********/


_CGWR& _CGWR::operator= ( const _CGWR &A ) {
    if (this != &A) {
     this->j=A.j;   
     this->j1=A.j1;
     this->j2=A.j2;
     this->j3=A.j3;
     this->l=A.l;   
     this->l1=A.l1;
     this->l2=A.l2;
     this->l3=A.l3;
     this->m=A.m;   
     this->m1=A.m1;
     this->m2=A.m2;
     this->m3=A.m3; 
     this->selected_Symb=A.selected_Symb;
    }
    return *this;
}


// Compute CG
/***************************************/
double long _CGWR::CG() {
    if (selected_Symb!=C_CG) {
        std::cerr << "_CGWR::CG(): /!\\ invalid symbol.\n";
        exit(EXIT_FAILURE);
    }

// Selection rules
/***************************************/
    bool test = abs(m1)<=j1;
    test &= abs(m2)<=j2;
    test &= abs(j1-j2)<=j;
    test &= j<=j1+j2;
    test &= abs(m)<=j;
    test &= m==m1+m2;
/***********/

    double long CG1=0, CG2=0;

    if (test) {
        CG1=sqrt(
            (2*j+1)
            * fact(j2+j-j1)
            * fact(j+j1-j2)
            * fact(j1+j2-j)
            * fact(j1+m1) * fact(j1-m1)
            * fact(j2+m2) * fact(j2-m2)
            * fact(j+m) * fact(j-m)
            / fact(j1+j2+j+1)
        );

        double long zmin, zmax;
        zmin=std::floor(std::max(
                std::max(static_cast<double long>(0.0),m2+j1-j),
                j2-j-m1));
        zmax=std::floor(
                std::min(j1+j2-j,
                std::min(j2+m2,j1-m1)));

        for(double long z=zmin; z<zmax+1; z++)
            CG2+=pow(-1,z) / fact(z)
                / fact(j1+j2-j-z)
                / fact(j+z-j1-m2)
                / fact(j2+m2-z)
                / fact(z+j-j2+m1)
                / fact(j1-m1-z);
    }
    return CG1*CG2;
}
/***********/


// Compute Wigner 3j
/***************************************/
double long _CGWR::W3j(void) {

    if (selected_Symb!=C_W3j) {
        std::cerr << "_CGWR::W3j(): /!\\ invalid symbol.\n";
        exit(EXIT_FAILURE);
    }

    double long res=0.0;
    double long sum=0.0;
    double long kmin, kmax;

    bool test = j1+j2>=j3;
    test &= j2+j3>=j1;
    test &= j3+j1>=j2;
    test &= m1+m2+m3==0;
    test &= abs(m1)>=0 && j1>=abs(m1);
    test &= abs(m2)>=0 && j2>=abs(m2);
    test &= abs(m3)>=0 && j3>=abs(m3);
    test &= is_integer(j1+j2+j3, j1-j2-m3);
    test &= is_integer(j1,j2,j3,m1,m2,m3);
    test |= is_halfinteger(j1,j2,j3,m1,m2,m3);

    if (test) {
        kmin=std::max(
                std::max(static_cast<double long>(0.0),j2-j3-m1),
                j1-j3+m2);
        kmax=std::min(std::min(j1+j2-j3,j1-m1),j2+m2);

        res = pow(-1,j1-j2-m3);
        res *= sqrt(
                fact(j1+j2-j3)
                * fact(j1-j2+j3)
                * fact(-j1+j2+j3)
                * fact(j1-m1) * fact(j1+m1)
                * fact(j2-m2) * fact(j2+m2)
                * fact(j3-m3) * fact(j3+m3)
                / fact(j1+j2+j3+1));

        for(double long k=kmin; k<kmax+1; k++)
            sum += pow(-1,k) / fact(k)
                    / fact(j1+j2-j3-k)
                    / fact(j1-m1-k)
                    / fact(j2+m2-k)
                    / fact(j3-j2+m1+k)
                    / fact(j3-j1-m2+k);
    }
    return res*sum;
}
/***********/


// Compute Wigner 6j
/***************************************/
double long _CGWR::W6j(void) {

    if (selected_Symb!=C_W6j) {
        std::cerr << "_CGWR::W6j(): /!\\ invalid symbol.\n";
        exit(EXIT_FAILURE);
    }

    
//Selection rules
// j_i has to obey to triangular delta below
/***************************************/
    auto delta = [] (double long a, 
                     double long b, 
                     double long c) {
        return (abs(a-b)<=c)&&(c<=a+b);
    };
    
    bool test=delta(j1,j2,j3);
    test &= delta(j1,l2,l3);
    test &= delta(l1,j2,l3);
    test &= delta(l1,l2,j3);
    
    test &= is_integer(j1+j2+j3);
    
/***********/

    double long res=0, sum=0;
    
    if (test) {
        double long kmin, kmax;
        kmin=std::max(
            std::max(j1+j2+j3, j1+l2+l3),
                      std::max(l1+j2+l3, l1+l2+j3) );
        kmax=std::min(
            std::min(j1+j2+l1+l2, j2+j3+l2+l3),
                      j3+j1+l3+l1);
        
        res= Delta(j1,j2,j3)
        * Delta(j1,l2,l3)
        * Delta(l1,j2,l3)
        * Delta(l1,l2,j3);
        
        for(double long k=kmin; k<kmax+1;k++)
            sum += pow(-1,k) * fact(k+1)
            / fact(k-j1-j2-j3)
            / fact(k-j1-l2-l3)
            / fact(k-l1-j2-l3)
            / fact(k-l1-l2-j3)
            / fact(j1+j2+l1+l2-k)
            / fact(j2+j3+l2+l3-k)
            / fact(j3+j1+l3+l1-k);
    }
    return res*sum;
}
/***********/


// Compute Racah symbol
/***************************************/
double long _CGWR::Racah(void) {
    return pow(-1,j1+j2+l1+l2)*W6j();
}
/***********/


// Int Factorial
/***************************************/
inline double long _CGWR::fact(double long n) {
    if (n<0 || !is_integer(n)) {
        std::cerr << "_CGWR::fact(): /!\\ bad number: "
                  << n << ".\n";
        exit(EXIT_FAILURE);
    }
    if (n==0) return 1;
    return n*fact(n-1);
    // ::math::factorial<double long >(abs(n));
}
/***********/


// Triangle function
/***************************************/
double long _CGWR::Delta(double long a,
                         double long b,
                         double long c) {
    double long res = fact(a+b-c)
                    * fact(a-b+c)
                    * fact(-a+b+c)
                    / fact(a+b+c+1);
    if (res<0) {
        std::cerr << "_CGWR::Delta(): /!\\ complex number: sqrt(" 
                  << res << ").\n";
        exit(EXIT_FAILURE);
    }
    return sqrt(res);
}
/***********/


// Kronecker symbol
/***************************************/
double long _CGWR::delta_k(double long a,
                           double long b) {
    return a==b;
}
/***********/


// Misc. methods
/***************************************/
template<typename _T>
bool _CGWR::is_integer(_T n) {
  return std::floor(n) == n;
}

template<typename _T, typename... _args>
bool _CGWR::is_integer(_T n, _args... args) {
    return is_integer(n) && is_integer(args...);    
}

template<typename _T>
bool _CGWR::is_halfinteger(_T s) {
  return std::floor(s*2) == 2*s && std::floor(s)!=s;
}

template<typename _T, typename... _args>
bool _CGWR::is_halfinteger(_T n, _args ...args) {
    return is_halfinteger(n) && is_halfinteger(args...);    
}

/***********/
