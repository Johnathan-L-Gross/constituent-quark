//
//  bar_calc.h
//  
//
//  Created by Johnathan Gross
//
//

#ifndef ____bar_calc__
#define ____bar_calc__

#ifdef D6_JLG
#define D3_JLG
#define _complex_JLG
#endif

#include "polynomials/polynomials.h"
#include "funch.h"
#include "int.h"
#include "sparse/sparse.h"

#ifdef _long_JLG
#warning "hamcalc only uses double precision for coupling constants"
#endif

constexpr marray mrot(const int i,
                      const marray&m){
    return {m[mod(i,3)],m[mod(1+i,3)],m[mod(2+i,3)]};
}

#ifndef _reduced_JLG
#define _reduced_JLG
//reduced masses
constexpr Double mrho(const marray&m){return m[0]*m[1]/(m[0]+m[1]);}
constexpr Double mlambda(const marray&m){
    return (m[0]+m[1])*m[2]/(m[0]+m[1]+m[2]);
}
#endif

/*ham is resultant matrix
 states is list of states
 params is collection of parameters related to integration
    polytype: Hermite (true) or Laguerre (false) quadrature
 		if(polytype) XX=r, else XX=r^2
 		if(polytype) dx=XW*XX^2, else dx=XW*XX/2
 omega is variational parameter
 mass is list of particle masses
 f is list terms to evaluate, functions of r and/or p
 k is list of variable type:
 	2s place: 0 for 1D lambda, 1 for 1D rho
 	1s place: 0 for position, 1 for momentum,
 		used in scalar, 1D spin-spin, and Dirac delta
 	ignored otherwise, still need to include
 tens is list of term type:
 	0 for scalar
 	1 for spin-spin
 	2 for spin tensor
 	3 for spin vector 2 body (do not include 1/rho in f)
 	4 for spin vector 3 body
 	5 for Dirac delta
 spin is list of spin configurations (formerly tens2), one of:
 	tens=0,2,3,4,5: ignored, still need to include a value to comply
 	tens=1: 11,22,33 for S_i^2
 		12,21,13,31,23,32 for S_i.S_j
 		only 12 and 21 allowed for 1D
 rot is list of Moshinsky rotation index
 pars is list of paramaters for use in the vector, ignored otherwise
 mm,nn terms to evaluate: both >=0 for element mmXnn, both <0 for full matrix
 */
void hamcalc(std::vector<sparse>&,//ham
             const imatrix&,//state
             const int1&,//params
             const Double,//omega
             const marray&,//mass
             const std::vector<func1>&,//f
             const ivector&,//k
             const ivector&,//tens
             const ivector&,//spin
             const ivector&,//rot
             const std::vector<marray>&,//pars
             const bool=false,//sym
             const Double=0);//error
inline void hamcalc(std::vector<sparse>&ham,
                    const imatrix&state,
                    const int1&params,
                    const Double omega,
                    const marray&mass,
                    const std::vector<func1>&f,
                    const ivector&k,
                    const ivector&tens,
                    const ivector&spin,
                    const ivector&rot,
                    const bool sym=false,
                    const Double error=0){
    //no pars
    hamcalc(ham,state,params,omega,mass,f,k,tens,spin,rot,
            std::vector<marray>(f.size(),{0,0,0}),sym,error);
}
//1D
#ifdef D3_JLG
void hamcalc(std::vector<sparse>&,//&ham
             const imatrix&,//states
             const int3&,//params
             const Double,//omega
             const marray&,//mass
             const std::vector<func3>&,//f
             const ivector&,//k
             const ivector&,//tens
             const ivector&,//spin
             const ivector&,//rot
             const bool=false,//sym
             const Double=0);//error
//3D
#endif
#ifdef D6_JLG
void hamcalc(std::vector<CDmatrix>&,//ham
             const imatrix&,//states
             const int6&,//params
             const Double,//omega
             const marray&,//mass
             const std::vector<func6>&,//f
             const ivector&,//k
             const ivector&,//tens
             const ivector&,//spin
             const ivector&,//rot
             const bool=false,//sym
             const Double=0);//error
//6D
#endif



#endif /* defined(____bar_calc__) */
