//
//  energy.h
//
//
//  Created by Johnathan Gross on 5/23/14.
//
//

#ifndef ____energy_JLG_
#define ____energy_JLG_

#include "sparse/sparse.h"
#include <functional>

#ifdef _long_JLG
#warning "Armadillo does not support long double types, all calculations were done with double"
#endif

void state_sort(tvector&,//ilist
                const Dvector&,//en
                const bool=false);//flag
/*sorts ilist into order of increasing energy as defined by en.
 Effect: stable_sort of ilist if flag is true, sort of ilist if flag is false.
 */

Double energy(const sparse&,//ham
              const std::size_t,//N
              const std::size_t=-1);//k
/*
 Finds Nth energy of ham using k eigenvalues.
 Output: Nth energy level
 */

Double energy(const sparse&,//ham
              Dvector&,//eval
              const std::size_t,//N
              const std::size_t=-1);//k
/*overload to store eigenvalues
 Output: Nth energy level
 Effect: modifies WF
 */

void energy(const sparse&,//ham
            sparse&,//evec
            Dvector&,//eval
            const std::size_t=-1);//k
/*overload to store eigenvalues and eigenvectors using k eigenvalues
 Effect: modifies evec and eval
 */

std::array<Double,4> minen(std::function<void(sparse&,Double)>&,//Hfunc
                           Double,//par
                           Double,//step
                           const std::size_t,//N
                           const std::size_t=-1,//k
                           Double=0,//error
                           const bool=true,//rel
                           const bool=true,//pos
                           const bool=false);//disp
/*Finds parameter par which minimizes Nth energy level
 Error is the fractional error if rel is true, absolute error if rel is false
 If disp, will output intermediate values to std::cerr
 Output: array with {min par, error par, min energy, error energy}
 Effect: calls Hfunc to produce Hamiltonian matrix
 */
inline std::array<Double,4> minen(std::function<void(Dmatrix&,Double)>&Hfunc,
                                  Double par,
                                  Double step,
                                  const std::size_t N,
                                  const std::size_t k=-1,
                                  Double error=0,
                                  const bool rel=true,
                                  const bool pos=true,
                                  const bool disp=false){
    std::function<void(sparse&,Double)> Hf=[&Hfunc](sparse&D,Double x){
        Dmatrix DD;
        Hfunc(DD,x);
        D=DD;
    };
    return minen(Hf,par,step,N,k,error,rel,pos,disp);
}
/*
 overload for dense matrix
 */
#endif //ifndef ____energy_JLG_
