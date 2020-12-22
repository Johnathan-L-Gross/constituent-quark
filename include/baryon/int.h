//
//  int.h
//  
//
//  Created by Johnathan Gross
//
//

#ifndef ____int__
#define ____int__

#ifdef D6_JLG
#define D3_JLG
#define _complex_JLG
#endif

#include "polynomials/polynomials.h"

inline Double phix(const Double x){return M_PI*(x+1);}

class int1{
    /*--------------------------------------------------------------------------
     Class to hold integration variables for Gaussian quadrature radial only

     returns values of the following variables:
     bool polytype: Hermite (true) or Laguerre (false) quadrature
     unsigned int N: max HO excitation
     std::size_t n: # of abscissa
     Double x(i), xw(i), dx(i): abscissa, weights, weight function
        if(polytype) XX=r, else XX=r^2
     	if(polytype) dx=XW*r^2, else dx=XW*r/2
     Double r(i) and r2(i) give r and r^2 agnositic to polytype
     Double lag(k,l,i): Laguerre polynomials
     	L_k^(l+1/2)(r^2)
     Double rad(k,l,i): radial wavefunction
     	R_kl(r^2)=N_nl*r^l*L_k^(l+1/2)(r^2) (no exp)
     Double harm(n,l): HO normalization N_nl
     
     Can only construct if given polytype, N, and reference to XX and XW.
     
     lag[x][l][k], lag(k,l,x)
     rad[x][k][l], rad(k,l,x)
     x from 0 to n-1
     l from 0 to N
     k from 0 to N/2
     -------------------------------------------------------------------------*/
public:
    //constructor
    int1(const bool,//polytype
         const unsigned int,//N
         const Dvector&,//XX
         const Dvector&);//XW

    inline bool polytype() const{return polytype_;}//get value of polytype
    inline unsigned int N() const{return N_;}//get value of N
    inline std::size_t n() const{return n_;}//get value of n

    //value access
    inline Double lag(std::size_t k,std::size_t l,std::size_t x) const{
#ifdef TESTING_STATES_JLG
        return lag_.at(x).at(l).at(k);
#else
        return lag_[x][l][k];
#endif
    }
    inline Double rad(std::size_t k,std::size_t l,std::size_t x) const{
#ifdef TESTING_STATES_JLG
        return rad_.at(x).at(k).at(l);
#else
        return rad_[x][k][l];
#endif
    }
    inline Double x(std::size_t x) const{
#ifdef TESTING_STATES_JLG
        return XX_.at(x);
#else
        return XX_[x];
#endif
    }
    inline Double r(std::size_t i) const{
        if(polytype_) return x(i);
        else return std::sqrt(x(i));
    }
    inline Double r2(std::size_t i) const{
        if(polytype_) return sqr(x(i));
        else return x(i);
    }
    inline Double xw(std::size_t x) const{
#ifdef TESTING_STATES_JLG
        return XW_.at(x);
#else
        return XW_[x];
#endif
    }
    inline Double dx(std::size_t x) const{
#ifdef TESTING_STATES_JLG
        return dx_.at(x);
#else
        return dx_[x];
#endif
    }
    inline Double harm(std::size_t n,std::size_t l) const{
#ifdef TESTING_STATES_JLG
        return harm_.at(n).at(l);
#else
        return harm_[n][l];
#endif
    }

protected:
    const unsigned int N_;

private:
    const bool polytype_;
    const std::size_t n_;
    const Dvector &XX_,&XW_;
    std::vector<Dmatrix> lag_,rad_;
    Dvector dx_;
    Dmatrix harm_;

    void setup();

};
//t is true for Hermite, false for Laguerre
#ifdef D3_JLG
class int3:public int1{
    /*--------------------------------------------------------------------------
     Class to hold integration variables for Gaussian quadrature in radial and
        altitude

     inherits from int1
     
     returns values of the following variables:
     std::size_t l: # of angular abscissa
     Double t(i), tw(i), td(i): angular abscissa, and weights
     Double Y(l,m,i), Y2(l,m,i): normalized associated Legendre polynomials
        P_l^m(ct) and P_l^m(-ct)=(-1)^(l+m)*P_l^m(ct) where ct=cos(theta)

     Can only construct if given pointers to LX and LW and either an instance
        of in1, or arguments to construct it.
     
     Y[th][ll][m], Y(ll,m,th)
     th from 0 to l-1
     ll from 0 to L
     m from 0 to ll
     -------------------------------------------------------------------------*/
public:
    //constructor from int1 and parameters
    int3(const int1&,//A
         const unsigned int,//L
         const Dvector&,//LX
         const Dvector&);//LW
    //constructor from parameters
    int3(const bool polytype,
         const unsigned int N,
         const Dvector&XX,
         const Dvector&XW,
         const unsigned int L,
         const Dvector&LX,
         const Dvector&LW):
    int3(int1(polytype,N,XX,XW),L,LX,LW){}

    inline std::size_t l() const{return l_;}//get value of l
    inline unsigned int L() const{return L_;}//get value of L

    //get value of Ylm
    inline Double Y(std::size_t ll,std::size_t m,std::size_t x) const{
#ifdef TESTING_STATES_JLG
        return Y_.at(x).at(ll).at(m);
#else
        return Y_[x][ll][m];
#endif
    }
    inline Double Y(std::size_t ll,std::size_t x) const{return Y(ll,0,x);}
    //m=0
    
    //get value of Ylm2
    inline Double Y2(std::size_t ll,std::size_t m,std::size_t x) const{
#ifdef TESTING_STATES_JLG
        return Y2_.at(x).at(ll).at(m);
#else
        return Y2_[x][ll][m];
#endif
    }
    inline Double Y2(std::size_t ll,std::size_t x) const{return Y2(ll,0,x);}
    //m=0
    
    
    inline Double t(std::size_t x) const{
#ifdef TESTING_STATES_JLG
        return LX_.at(x);
#else
        return LX_[x];
#endif
    }
    inline Double tw(std::size_t x) const{
#ifdef TESTING_STATES_JLG
        return LW_.at(x);
#else
    	return LW_[x];
#endif
    }
    inline Double dt(std::size_t x) const{return tw(x);}
    
    inline Double cos(std::size_t x)const{return t(x);}
    inline Double sin(std::size_t x)const{return std::sqrt(1-sqr(t(x)));}
    inline Double theta(std::size_t x)const{return std::acos(t(x));}

protected:
    const unsigned int L_;

private:
    const std::size_t l_;
    const Dvector &LX_,&LW_;
    std::vector<Dmatrix> Y_,Y2_;

    void setup();
};
#endif
#ifdef D6_JLG
class int6:public int3{
    /*--------------------------------------------------------------------------
     Class to hold integration variables for 3D Gaussian quadrature

     inherits from int3

     returns values and const pointers to the following variables:
	 std::size_t m: # phase abscissa
     std::size_t m2: # of positive phase abscissa
     Double p(i), pw(i): phase abscissa and weights
     Double dphi(i): phase weight function
     Double phi(i): actual values of points
     CDouble phase(M,i): exp(I M phi)/sqrt(2Pi) and complex
        conjugate, phase(-M,i)=phase2[x][M]

     Can only construct if given pointers to PX and PW and either an instance
        of int3, or arguments to construct it.

     phase[x][M], phase(M,x)
     x from 0 to m
     M from 0 to L
     -------------------------------------------------------------------------*/
public:
    //constructor from int3 and parameters
    int6(const int3&,//A
         const Dvector&,//PX
         const Dvector&);//PW
    //constructor from int1 and parameters
    int6(const int1&A,
         const unsigned int L,
         const Dvector&LX,
         const Dvector&LW,
         const Dvector&PX,
         const Dvector&PW):
    int6(int3(A,L,LX,LW),PX,PW){}
    //constructor from parameters
    int6(const bool polytype,
         const unsigned int N,
         const Dvector&XX,
         const Dvector&XW,
         const unsigned int L,
         const Dvector&LX,
         const Dvector&LW,
         const Dvector&PX,
         const Dvector&PW):
    int6(int3(int1(polytype,N,XX,XW),L,LX,LW),PX,PW){}

    inline std::size_t m() const{return m_;}//get value of m
    inline std::size_t m2() const{return m2_;}//get value of m2

    inline Double p(std::size_t x) const{
#ifdef TESTING_STATES_JLG
        return PX_.at(x);
#else
        return PX_[x];
#endif
    }//get value of PX
    inline Double pw(std::size_t x) const{
#ifdef TESTING_STATES_JLG
        return PW_.at(x);
#else
        return PW_[x];
#endif
    }//get value of PW
    inline Double phi(std::size_t x) const{
#ifdef TESTING_STATES_JLG
        return phi_.at(x);
#else
        return phi_[x];
#endif
    }//get value of phi
    inline Double dphi(std::size_t x) const{
#ifdef TESTING_STATES_JLG
        return dphi_.at(x);
#else
        return dphi_[x];
#endif
    }//get value of dphi
    inline CDouble phase(int mm,std::size_t x) const{
#ifdef TESTING_STATES_JLG
        return (mm>0?phase_.at(x).at(mm):phase2_.at(x).at(-mm));
#else
        return (mm>0?phase_[x][mm]:phase2_[x][-mm]);
#endif
    }//get value of phase


private:
    const std::size_t m_,m2_;
    const Dvector &PX_,&PW_;
    Dvector dphi_,phi_;
    CDmatrix phase_,phase2_;

    void setup();
};
#endif

#endif /* defined(____int__) */
