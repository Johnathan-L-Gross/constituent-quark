//
//  quark.h
//
//
//  Created by Johnathan Gross
//

#ifndef quark_h
#define quark_h

#define D3_JLG
#include "baryon/funch.h"
#include "sparse/sparse.h"
#include "baryon/int.h"
#include "baryon/wf.h"

#ifdef _long_JLG
#warning "calculation only uses double precision for coupling constants"
#endif

typedef std::function<void(sparse&,const Double)> hfunc;

hfunc strong(const std::vector<std::pair<Double,Double>>&,//a
             const Double,//b
             const Double,//sigma0
             const Double,//s
             const Double,//eCoul
             const Double,//eCont
             const Double,//eTen
             const wf2&,//wf
             const int3&,//par
             const Double=0);//error
/*
 Calculates the hamiltonian matrix for the strong terms
 kinetic energy + 3 body confinement + Coulomb + Contact + Tensor
 */
inline hfunc strong(const std::vector<std::pair<Double,Double>>&a,
                    const Double b,
                    const Double sigma0,
                    const Double s,
                    const Double eCoul,
                    const Double eCont,
                    const wf2&wf,
                    const int3&par,
                    const Double error=0){
    //overload for eCont=eTen
    return strong(a,b,sigma0,s,eCoul,eCont,eCont,wf,par,error);
}

hfunc electric(const Double,//aEM
               const Double,//gEM
               const marray&,//ch
               const wf2&,//wf
               const int1&,//par
               const Double=0);//error
/*
 Calculates the hamiltonian matrix for the electric Coulomb terms
 */
hfunc electromagnetic(const Double,//aEM
                      const Double,//gEM
                      const Double,//eEM
                      const marray&,//ch
                      const wf2&,//wf
                      const int1&,//par
                      const Double=0);//error
/*
 Cacluates the hamiltonian matrix for the electromagnetic terms
 Coulomb + magnetic
 */

inline hfunc stel(const std::vector<std::pair<Double,Double>>&a,
                  const Double aEM,
                  const Double b,
                  const Double sigma0,
                  const Double s,
                  const Double gEM,
                  const Double eCoul,
                  const Double eCont,
                  const Double eTen,
                  const marray&ch,
                  const wf2&wf,
                  const int3&par,
                  const Double error=0){
    //calculates strong and electric terms
    return [=,&a,&ch,&wf,&par](sparse&ham,const Double omega){
        sparse h1(wf.size1(),wf.size1(),error),
        h2(wf.size1(),wf.size1(),error);
        strong(a,b,sigma0,s,eCoul,eCont,eTen,wf,par,error)(h1,omega);
        electric(aEM,gEM,ch,wf,par,error)(h2,omega);
        ham=h1+h2;
        ham=(ham+ham.transpose())/2;
        ham.set_error(error,true);
    };
}

inline hfunc stelmag(const std::vector<std::pair<Double,Double>>&a,
                     const Double aEM,
                     const Double b,
                     const Double sigma0,
                     const Double s,
                     const Double gEM,
                     const Double eCoul,
                     const Double eCont,
                     const Double eTen,
                     const Double eEM,
                     const marray&ch,
                     const wf2&wf,
                     const int3&par,
                     const Double error=0){
    //calcualtes strong and electromagnetic terms
    return [=,&a,&ch,&wf,&par](sparse&ham,const Double omega){
        sparse h1(wf.size1(),wf.size1(),error),
        h2(wf.size1(),wf.size1(),error);
        strong(a,b,sigma0,s,eCoul,eCont,eTen,wf,par,error)(h1,omega);
        electromagnetic(aEM,gEM,eEM,ch,wf,par,error)(h2,omega);
        ham=h1+h2;
        ham=(ham+ham.transpose())/2;
        ham.set_error(error,true);
    };
}

#endif // quark_h
