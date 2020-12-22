//	
//  Moshinsky.h
//
//
//  Created by Johnathan Gross
//
//

#ifndef _Moshinsky_JLG_h
#define _Moshinsky_JLG_h

#ifdef _long_JLG
typedef long double Double;
#else
typedef double Double;
#endif //_long_JLG

Double tm3(const int m1,const int L1,const int m2,const int L2,
           //ket |l,m1,L1,m2,L2>
           const int n1,const int l1,const int n2,const int l2,const int l,
           //bra <l,n1,l1,n2,l2|
           const Double X,const Double Y,const Double Z,const int k); //_k
/*
 Calculate 3 body Moshinsky brackets.
 k=3 corresponds to alpha'=alpha=123
 k=1 corresponds to alpha'=beta=231
 k=2 corresponds to alpha'=gamma=312
 _k<l,n1,l1,n2,l2|l,m1,L1,m2,L2>_3
 */

#ifdef TM4_JLG
//don't need tm4 for baryons, so not including it
Double tm4(const int n1,const int l1,
           const int n2,const int l2,
           const int n3,const int l3,
           //ket |>
           const int m1,const int L1,
           const int m2,const int L2,
           const int m3,const int L3,
           //bra <|
           const int l0,const int l0t,const int l,
           const Double X1,const Double X2,const Double X3,const Double X4,
           const int K);
/*
 Calculate 4 body Moshinsky brackets.
 K=k1,k2,i1,i2,i3 has allowed values
 11123, 11124, 11132, 11134, 11142, 11143, 11231, 11234, 11241, 11243, 11341, 11342,
 12123, 12132, 12142,
 21123, 21124, 21132, 21134, 21142, 21143, 21231, 21234, 21241, 21243, 21341, 21342,
 22123, 22132, 22142.
 */
#endif

Double spin3(const int S,
             const int P1,
             const int P2,
             int k)
//_k<S,P2|S,P1>_3
noexcept(
#ifdef TESTING_STATES_JLG
         false
#else
         true
#endif
         );
/*
 Calculates the transformation matrix between spin wavefunctions for 3 Spin-1/2
 particles.
 Same rule for k applies as for tm3
 */

#endif //ifndef _Moshinsky_JLG_h
