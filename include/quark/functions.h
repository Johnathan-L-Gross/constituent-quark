//
//  functions.h
//  
//
//  Created by Johnathan Gross
//

#ifndef functions_h
#define functions_h

#include "baryon/funch.h"
#include "polynomials/polynomials.h"
#include <cmath>
#include <stdexcept>

inline constexpr long double Mpi_3=1.047197551196597746154214461093167628065723133125035273658L;
inline constexpr long double M2pi_3=2.094395102393195492308428922186335256131446266250070547316L;
inline constexpr long double M_spi=0.5641895835477562869480794515607725858440506293289988568440857217106424684414934144867436602021073636L;
inline constexpr long double M16pi_9=5.58505360638185464615581045916356068301719004333352145951101260854722916673103821978317302283043034307937437602361L;
inline constexpr long double M16_9spi=1.00300370408490006568547458055248459705609000769599796772281906081891994389598829242087761813707975723205039627797L;
inline constexpr long double M8pi_3=8.37758040957278196923371568874534102452578506500028218926651891282084375009655732967475953424564551461906156403541L;
inline constexpr long double M4_9spi=0.25075092602122501642136864513812114926402250192399949193070476520472998597399707310521940453426993930801259906949L;

constexpr Double cos12(const Double cosx,
                       const Double cosy,
                       const Double px,
                       const Double py){
    return std::sqrt((1.-sqr(cosx))*(1.-sqr(cosy)))*std::cos(px-py)+cosx*cosy;
}
constexpr Double dot(const Double x,
                     const Double y,
                     const Double c){
    return x*y*c;
}
constexpr Double dot(const Double x,
                     const Double tx,
                     const Double px,
                     const Double y,
                     const Double ty,
                     const Double py){
    return dot(x,y,cos12(tx,ty,px,py));
}

#ifndef _reduced_JLG
#define _reduced_JLG
//reduced masses
constexpr Double mrho(const marray&m){return m[0]*m[1]/(m[0]+m[1]);}
constexpr Double mlambda(const marray&m){
    return (m[0]+m[1])*m[2]/(m[0]+m[1]+m[2]);
}
#endif

constexpr Double fermat_length(const Double r,
                               const Double l,
                               const Double c,
                               const marray&m){
    //The Fermat point is the point that connects to all three points of a
    //triangle with minimum distance. If one angle of the triangle is greater
    //than 120 degrees, then the length is the sum of the two adjacent sides
    //returns the total length
    const Double x=m[1]/(m[0]+m[1]),
    y=m[0]/(m[0]+m[1]);
    const Double t1=std::atan2(l*std::sqrt(1-c*c),x*r+l*c),
    t2=atan2(l*std::sqrt(1-c*c),y*r-l*c);
    const Double r13=std::sqrt(x*x*r*r+l*l+2*x*r*l*c),
    r23=std::sqrt(y*y*r*r+l*l-2*y*r*l*c);

    if(t1>=M2pi_3)
        return r+r13;
    else if(t2>=M2pi_3)
        return r+r23;
    else if(t1+t2<=Mpi_3)//t3>=M2pi_3
        return r13+r23;
    else
        return std::sqrt(r13*r13+r*r-2*r13*r*std::cos(t1+Mpi_3));
}

constexpr Double pairwise_length(const Double r,
                                 const Double l,
                                 const Double c,
                                 const marray&m){
    const Double x=m[1]/(m[0]+m[1]),
    y=m[0]/(m[0]+m[1]);
    const Double r13=std::sqrt(x*x*r*r+l*l+2*x*r*l*c),
    r23=std::sqrt(y*y*r*r+l*l-2*y*r*l*c);
    return r+r13+r23;
}

//smearing sigma value
constexpr Double sigmaij(const Double s0,
                         const Double s,
                         const marray&m){
    //full Capstick-Isgur sigma_ij
    if(0!=s0 and 0!=s)
        return std::sqrt(s0*s0*(1+pow(4*m[0]*m[1]/sqr(m[0]+m[1]),4))/2
                         +s*s*sqr(2*m[0]*m[1]/(m[0]+m[1])));
    else if(0!=s) return 2*s*m[0]*m[1]/(m[0]+m[1]);
    else if(0!=s0) return s0;
    else return 0;
}
constexpr Double sigmaij2(const Double s0,
                          const Double s,
                          const marray&m){
    //reduced Gross-Capstick sigma_ij
    if(0!=s0 and 0!=s)
        return std::sqrt(sqr(s0)+sqr(2*s*m[0]*m[1]/(m[0]+m[1])));
    else if(0!=s) return 2*s*m[0]*m[1]/(m[0]+m[1]);
    else if(0!=s0) return s0;
    else return 0;
}
constexpr Double sigmaijk(const Double sij,
                          const Double gk){
    return sij*gk/std::sqrt(sqr(gk)+sqr(sij));
}

constexpr Double gamma1(const Double a1,
                        const Double a0,
                        const Double g0,
                        const Double x){
    //returns the value of g1 such that
    //a1*exp(-x^2/(4g1^2))=a0*exp(-x^2/(4g0^2))
    return pow(4/sqr(x)*std::log(a1/a0)+1/sqr(g0),-.5);
}

//identity function
constexpr Double ident(const Double,const marray&){return 1;}

//funcitons that return the (non)relativistic kinetic energy
constexpr Double relkin(const Double x,
                        const marray&m){
    return sqrt(sqr(x)+sqr(m[2]));
}
constexpr Double nonrelkin(const Double x,
                           const marray&m){
    return sqr(x)/(2*m[2]);
}
constexpr Double nonrelrho(const Double x,
                           const marray&m){
    return sqr(x)/(2*mrho(m));
}
constexpr Double nonrellambda(const Double x,
                              const marray&m){
    return sqr(x)/(2*mlambda(m));

}
/*
constexpr Double nonreltotrho(const Double x,
                              const marray&m){
    return x*x*(1/m[0]+1/m[1])/2;
}
constexpr Double nonreltotlam(const Double x,
                              const marray&m){
    return x*x*(1/m[2]+1/(m[0]+m[1]))/2;
}
 */
//repeat of nonrelrho and nonrellambda

constexpr Double HOE(const int n,const int l){return 2*n+l+Double(1.5);}
//Harmonic oscillator energy level

//functions that take the HO frequency and return the harmonic oscillator potential
inline func1 HO(const Double omega){
    return [=](const Double x,const marray&m){
        return mrho(m)*sqr(omega*x)/2;
    };
}
inline func1 HOrho(const Double omega){return HO(omega);}
inline func1 HOlambda(const Double omega){
    return [=](const Double x,const marray&m){
        return mlambda(m)*sqr(x*omega)/2;
    };
}
inline func1 HOtotrho(const Double omega){
    return [=](const Double x,const marray&m){
        return omega*m[0]*(1+m[2]/(m[0]+m[2]))*sqr(x)/4;
    };
}
inline func1 HOtotlam(const Double omega){
    return [=](const Double x,const marray&m){
        return omega*m[0]*m[2]/(m[0]+m[2])*sqr(x);
    };
}

//linear functions
inline func1 lin(Double b){
    return [=](const Double x,const marray&m[[maybe_unused]]){
        return b*x;
    };
}
inline func1 smearlin(const Double b,
                      const Double s){
    return [=](const Double x,
               const marray&m[[maybe_unused]]){
        const Double sr=s*x;
        return b/s * ((sr+0.5/sr)*std::erf(sr)+std::exp(-sqr(sr))*M_spi);
    };
}
#ifdef D3_JLG
//returns total nonrelativistic pairwise confining potential
inline func3 lintot(const Double b){
    return [=](const Double x,
               const Double y,
               const Double t,
               const marray&m){
        return b*pairwise_length(x,y,t,m);
    };
}
//function that returns the nonrelativistic linear confining ponential
inline func3 conf(const Double b){
    return [=](const Double x,
               const Double y,
               const Double t,
               const marray&m){
        return b*fermat_length(x,y,t,m);
    };
}
#endif

//color Coulomb
inline func1 coulomb(const Double a){
    return [=](const Double x,
               const marray&m[[maybe_unused]]){
        return -a/(1.5*x);
    };
}
inline func1 smearcoulomb(const Double a,const Double s){
    return [=](const Double x,
               const marray&m[[maybe_unused]]){
        return -a*std::erf(s*x)/(1.5*x);
    };
}
inline func1 fullcoulomb(const std::vector<std::pair<Double,Double>>&a,
                         const Double s){
    return [s,&a](const Double x,
                  const marray&m[[maybe_unused]]){
        Double result=0;
#pragma omp parallel for reduction(+:result)
        for(std::size_t k=0;k<a.size();k++){
            const Double xk=x*sigmaijk(s,a[k].second);
            result+=a[k].first*erf(xk);
        }
        return result/(-1.5*x);
    };
}

//EM Coulomb
inline func1 EMcoulomb(const Double a,
                       const marray&ch){
    return [a,&ch](const Double x,
                   const marray&m[[maybe_unused]]){
        return a*ch[0]*ch[1]/x;
    };
}
inline func1 EMsmearcoulomb(const Double a,
                            const Double g,
                            const marray&ch){
    return [a,g,&ch](const Double x,
                   const marray&m[[maybe_unused]]){
        return a*ch[0]*ch[1]*std::erf(g*x)/x;
    };
}

//color contact spin-spin
inline func1 contact(const Double a){
    return [a](const Double x[[maybe_unused]],
               const marray&m){
        return M16pi_9*a/(m[0]*m[1]);
    };
}
inline func1 smearcontact(const Double a,
                          const Double s){
    return [a,s](const Double x,
                 const marray&m){
        return M16_9spi*a/(m[0]*m[1])*std::exp(-sqr(s*x))*s*s*s;
    };
}
inline func1 fullcontact(const std::vector<std::pair<Double,Double>>&a,
                         const Double s){
    return [s,&a](const Double x,const marray&m){
        Double result=0;
#pragma omp parallel for reduction(+:result)
        for(std::size_t k=0;k<a.size();k++){
            const Double sk=sigmaijk(s,a[k].second);
            result+=a[k].first*std::pow(sk,3)*std::exp(-sqr(sk*x));
        }
        return result*M16_9spi/(m[0]*m[1]);
    };
}

//EM contact spin-spin
inline func1 EMcontact(const Double a,
                       const marray&ch){
    return [a,&ch](const Double x[[maybe_unused]],
                   const marray&m){
        return -a*M8pi_3*ch[0]*ch[1]/(m[0]*m[1]);
    };
}
inline func1 EMsmearcontact(const Double a,
                            const Double g,
                            const marray&ch){
    return [a,g,&ch](const Double x,
                   const marray&m){
        return -M8pi_3*a*ch[0]*ch[1]/(m[0]*m[1])*std::exp(-sqr(g*x))*std::pow(g*M_spi,3);
    };
}

//tensor
inline func1 tensor(const Double a){
    return [a](const Double x,
               const marray&m){
        return a/(1.5*m[0]*m[1]*std::pow(x,3));
    };
}
inline func1 fulltensor(const std::vector<std::pair<Double,Double>>&a,
                        const Double s){
    const Double s2=s*s;
    return [s,s2,&a](const Double x,
                     const marray&m){
        Double result1=0,
 		result2=0;
#pragma omp parallel for reduction(+:result1,result2)
        for(std::size_t k=0;k<a.size();k++){
            const Double xk=x*s*a[k].second/std::sqrt(sqr(a[k].second)+s2);
            result1+=a[k].first*std::erf(xk);
            result2+=a[k].first*xk*(2*sqr(xk)+3)*std::exp(-sqr(xk));
        }
        return (result1/1.5-M4_9spi*result2)/(m[0]*m[1]*x*x*x);
    };
}

//spin-orbit
inline func1 SOcm(const Double a){
    return [a](const Double x,
               const marray&m[[maybe_unused]]){
        return a/(-1.5*std::pow(x,3));
    };
}
/*inline fucn1 SOt(const Double a){
    return [=](const double x,const marray&m[[maybe_unused]]){
        return
    };
}*/
/*inline fucn1 SOv(const Dvector&a,const Dvector&g,const Double s){
 	if(a.size()!=g.size()) throw std::invalid_argument("a and g not same length");
    return [&a,&g,=](const double x,const marray&m){
        return
    };
}*/
/*inline fucn1 SOs(const Dvector&a,const Dvector&g,const Double s){
 	if(a.size()!=g.size()) throw std::invalid_argument("a and g not same length");
    return [&a,&g,=](const double x,const marray&m){
        return
    };
}*/

//relativistic suppression terms
inline func1 beta(const Double eps){
    if(-0.5==eps) return ident;
    else return [eps](const Double p,
                      const marray&m){
        const Double pp=p*p;
        if(m[0]==m[1]) return std::pow(1+1/(1+sqr(m[0])/pp),0.5+eps);
        else return std::pow(1+pp/std::sqrt((pp+sqr(m[0]))*(pp+sqr(m[1]))),
                        0.5+eps);
    };
}
inline func1 delta(const Double eps){
    if(-0.5==eps) return ident;
    else return [eps](const Double p,
                      const marray&m){
        const Double pp=p*p;
        if(m[0]==m[1]) return std::pow(sqr(m[0])/(pp+sqr(m[0])),0.5+eps);
        else return std::pow(m[0]*m[1]/std::sqrt((pp+sqr(m[0]))*(pp+sqr(m[1]))),
                        0.5+eps);
    };
}
#endif /* functions_h */
