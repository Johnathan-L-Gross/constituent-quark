//
//  states.h
//
//
//  Created by Johnathan Gross
//
//

/*-------------------------------------------------------------
 J AND S VALUES ARE TWICE THEIR PHYSICAL VALUES
 J=1/2 is stored as J=1
 J=3/2 is stored as J=3
 This was done to store them as integers.
 
 Positive parity is stored as 0 for (-1)^0 and negative as 1 for (-1)^1
 
 states are labeled in this order:
 Pf,J,Pj,M,N,L,np,lp,nl,ll,S,Ps
 0  1 2  3 4 5 6  7  8  9  10 11
 */



#ifndef ____states_JLG_
#define ____states_JLG_

#include "sparse/sparse.h"
#include "polynomials/polynomials.h"
#include <cmath>
#ifdef TESTING_STATES_JLG
#include <iostream>
#endif

typedef std::array<Double,3> marray;

inline constexpr int SMIN_B=1;
inline constexpr int SMAX_B=3;

void statematrix(const int,//Pf
                 const int,//J
                 const int,//P
                 const int,//M
                 const int,//Nmin
                 const int,//Nmax
                 imatrix&,//statelist
                 const bool=true);//sym
/*
 Lists all 3D HO states with 1-2 flavor symmetry Pf, total angular momentum J,
    parity P, and projection M. Uses S=1/2 and S=3/2 spins. Will return states
    with energy quantum number between Nmin and Nmax, inclusive. If sym,
    will only give states that are antisymmetric in 1-2 (in color, flavor,
    spin, and space).
 If sort, will double check with std::sort.
 
 Pf,J,P,M,N,L,np,lp,nl,ll,S,Ps
 0 ,1,2,3,4,5,6 ,7 ,8 ,9,10,11
 Sorts by increasing N, increasing L, increasing lp+ll, increasing np,
    increasing lp, increasing S, antisym before sym,
 Effect: empties and replaces statelist with list of labels for each state
 */

inline void statematrixm(const int Pf,
                         const int J,
                         const int P,
                         const int Nmin,
                         const int Nmax,
                         imatrix&statelist,
                         const bool sym=true){
    statematrix(Pf,J,P,J,Nmin,Nmax,statelist,sym);
    return;
}
//pseudo-overload statematrix without M projection
inline void statematrix(const int Pf,
                        const int J,
                        const int P,
                        const int M,
                        const int Nmax,
                        imatrix&statelist,
                        const bool sym=true){
    statematrix(Pf,J,P,M,0,Nmax,statelist,sym);
    return;
}
//overload statematrix without Nmin
inline void statematrix(const int Pf,
                        const int J,
                        const int P,
                        const int Nmax,
                        imatrix&statelist,
                        const bool sym=true){
    statematrix(Pf,J,P,J,0,Nmax,statelist,sym);
    return;
}
//overload statematrix without Nmin or M projection

void combinestates(const std::vector<imatrix>&,//sin
                   imatrix&);//sout
/*
 Combines state matrices into one.
 Effect: replaces sout with a state matrix combining all of sin
 */

void changeM(const int,//M
             imatrix&);//statelist
/*
 Changes M value of states in statelist.
 Effect: changes value of 3rd element of each item in statelist
 */


inline int sym12(const ivector&s){return mod(1+s[0]+s[7]+s[11],2);}
/*
 Returns the total 1-2 symmetry of the wavefunction
 Output: 1 if antisymmetric and 0 if symmetric
 */
inline int sym12ss(const ivector&s){return mod(s[7]+s[11],2);}
/*
 Returns the spin-space 1-2 symmetry of the wavefunction
 Output: 1 if antisymmetric and 0 if symmetric
 */

//******************************************************************************
//Moshinsky

void transform(const imatrix&,//s1 outer state
               const imatrix&,//s2 inner state
               const marray&,//mass unrotated masses
               const int,//k rotation index
               sparse&);//trans transformation matrix
/*Creates transformation matrix between two states using Moshinsky brackets.
 Assumes that flavor wave function is either totally symmetric or totally
 antisymmetric, no mixed symmetry
 Effect: resizes trans and assigns values
 */
inline void transform(const imatrix&s,//state
                      const marray&mass,//unrotated masses
                      const int k,//rotation index
                      sparse&trans){//transformation matrix
	transform(s,s,mass,k,trans);
}
//Same as above, but between states in the same basis.
void transformss(const imatrix&,//s1 outer state
                 const imatrix&,//s2 inner state
                 const marray&,//mass unrotated masses
                 const int,//k rotation index
                 sparse&);//trans transformation matrix
/*Creates transformation matrix between two states using Moshinsky brackets.
 Only applies spin-space transformations
 Effect: resizes trans and assigns values
 */
inline void transformss(const imatrix&s,//state
                        const marray&mass,//unrotated masses
                        const int k,//rotation index
                        sparse&trans){//transformation matrix
    transformss(s,s,mass,k,trans);
}
//Same as above, but between states in the same basis.

#endif //____states_JLG_
