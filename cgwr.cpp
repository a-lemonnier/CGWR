#include "cgwr.hpp"

const unsigned int _CGWR::C_CG    =   1;
const unsigned int _CGWR::C_W3j   =   2;
const unsigned int _CGWR::C_W6j   =   3;
const unsigned int _CGWR::C_W9j   =   4;
const unsigned int _CGWR::C_RACAH =   5;

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



_CGWR& _CGWR::operator= ( const _CGWR& other ) {
    //TODO
    return *this;
}


double long _CGWR::CG() {
    if (selected_Symb!=C_CG) {
        std::cerr << "_CGWR::CG(): /!\\ invalid symbol.\n";
        exit(EXIT_FAILURE);
    }

    bool test = abs(m1)<=j1;
    test &= abs(m2)<=j2;
    test &= abs(j1-j2)<=j;
    test &= j<=j1+j2;
    test &= abs(m)<=j;
    test &= m==m1+m2;

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

    test &= is_integer(j1+j2+j3);
    test &= is_integer(j1-j2-m3);

    test &= m1+m2+m3==0;

    test &= abs(m1)>=0 && j1>=abs(m1);
    test &= abs(m2)>=0 && j2>=abs(m2);
    test &= abs(m3)>=0 && j3>=abs(m3);

    test &= is_integer(j1) || is_halfinteger(j1);
    test &= is_integer(j2) || is_halfinteger(j2);
    test &= is_integer(j3) || is_halfinteger(j3);
    test &= is_integer(m1) || is_halfinteger(m1);
    test &= is_integer(m2) || is_halfinteger(m2);
    test &= is_integer(m3) || is_halfinteger(m3);

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

double long _CGWR::W6j(void) {

    if (selected_Symb!=C_W6j) {
        std::cerr << "_CGWR::W6j(): /!\\ invalid symbol.\n";
        exit(EXIT_FAILURE);
    }

    //Selection rules
    //TODO
    
    double long res, sum=0.0;
    double long kmin, kmax;
    kmin=std::max(
            std::max(j1+j2+j3, j1+l2+l3),
            std::max(l1+j2+l3, l1+l2+j3) );
    kmax=std::min(
            std::min(j1+j2+l1+l2, j2+j3+l2+l3),
            j3+j1+l3+l1);
    
    res= Delta(j1,j2,j3)
                * Delta(j1,l2,l3)
                * Delta(l1,j2,j3)
                * Delta(l1,l2,j3);
    //std::cout << "\n" <<   << "\n";

    for(double long k=kmin; k<kmax+1;k++)
        sum += pow(-1,k) * fact(k+1)
                / fact(k-j1-j2-j3)
                / fact(k-j1-l2-l3)
                / fact(k-l1-j2-l3)
                / fact(k-l1-l2-j3)
                / fact(j1+j2+l1+l2-k)
                / fact(j2+j3+l2+l3-k)
                / fact(j3+j1+l3+l1-k);

    return res*sum;
}

double long _CGWR::Racah(void) {
    return pow(-1,j1+j2+l1+l2)*W6j();
}

double long _CGWR::fact(double long n) {
    if (n<0 || !is_integer(n)) {
        std::cerr << "_CGWR::fact(): /!\\ bad number: "
                  << n << ".\n";
        exit(EXIT_FAILURE);
    }
    double long res=1;
    for(double long i=1;i<n+1;i++)
        res*=i;
    return res;
    // ::math::factorial<double long >(abs(n));
}

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

double long _CGWR::delta_k(double long a,
                            double long b) {
    double long res=0;
    if (a==b) res=1;
    return res;
}

bool _CGWR::is_integer(double long n) {
  return std::floor(n) == n;
}

bool _CGWR::is_halfinteger(double long s) {
  return std::floor(s*2) == 2*s && std::floor(s)!=s;
}
