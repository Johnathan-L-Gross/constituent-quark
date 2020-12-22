//	
//  algebra.h
//  
//
//  Created by Johnathan Gross
//
//

#ifndef ____algebra__
#define ____algebra__

#ifdef W_JLG_extern
#ifndef W_JLG_map
#define W_JLG_map
#endif
#endif

#include "polynomials/polynomials.h"
#include <cmath>
#include <functional>
#include <iostream>

#ifdef W_JLG_map
#include <array>
#include <map>
typedef std::array<unsigned int,2> key2;
typedef std::array<int,5> key5;
typedef std::array<unsigned int,6> key6;
typedef std::array<unsigned int,9> key9;
typedef std::map<key2,double> map2;
typedef std::map<key5,double> map3;
typedef std::map<key6,double> map6;
typedef std::map<key9,double> map9;
#endif

#ifdef W_JLG_extern
extern map3 w3j_map_JLG;
extern map6 w6j_map_JLG;
extern map9 w9j_map_JLG;
extern map2 norm_map_JLG;
extern map2 norml_map_JLG;
#endif

#ifdef _long_JLG
#warning "GNU Scientific Library calculated with double precision, not long double"
#endif


constexpr unsigned int abs_diff(const unsigned int i,
                                const unsigned int j){
    return (i<j)?j-i:i-j;
}
constexpr bool istriangle(const unsigned int i,
                          const unsigned int j,
                          const unsigned int k){
    return i+j>=k and abs_diff(i,j)<=k and !((i+j+k)%2);
}

/*-----------------------------------------------------------------------------
 Wigner 3-j symbols
 / j1/2 j2/2 j3/2 \
 \ m1/2 m2/2 m3/3 /
 
 Clebsch-Gordan coeffidients
 C^{J/2 M/2}_{j1 m1 j2 m2}=(-1)^((j2-j1-M)/2) * sqrt(J+1) * / j1/2 j2/2  J/2 \
                                                            \ m1/2 m2/2 -M/2 /
 -----------------------------------------------------------------------------*/
#ifdef W_JLG_map
double W3J(const unsigned int,//j1
           const unsigned int,//j2
           const unsigned int,//j3
           const int,//m1
           const int,//m2
           const bool=false);//th
/*Calculateds Wigner 3-j symbol. Caution, all Js and Ms are doubled.
 Output: Wigner 3-j symbol value
 Effect: (if not already created) creates static/global map
 */
inline double W3J(const unsigned int j1,
                  const unsigned int j2,
                  const unsigned int j3,
                  const int m1,
                  const int m2,
                  const int m3,
                  const bool th=false){
    if(m1+m2+m3){
#ifdef TESTING_STATES_JLG
        if(th) throw std::invalid_argument("projections don't add up");
        else return 0.0;
#else
        return 0.0;
#endif
    }
    else return W3J(j1,j2,j3,m1,m2);
}
/*overload of W3J to add m3 (since m3 is redundant)
 Output: Wigner 3-j symbol value
 Effect: (if not already created) creates static/global map
 */
inline double W3J(const key5 k,const bool th=false){
    return W3J(k[0],k[1],k[2],k[3],k[4],th);
}
/*overload to allow using key5 type as argument. Caution, all Js and Ms are doubled.
 Output: Wigner 3-j symbol value
 Effect: (if not already created) creates static/global map
 */
#endif //W_JLG_map

double Wigner3(const unsigned int,//j1
               const unsigned int,//j2
               const unsigned int,//j3
               const int,//m1
               const int);//m2
/*Calculateds Wigner 3-j symbol. Caution, all Js and Ms are doubled.
 Similar to W3J, but doesn't create static/global map
 Output: Wigner 3-j symbol value
 Effect: none
 */
inline double Wigner3(const unsigned int j1,
                      const unsigned int j2,
                      const unsigned int j3,
                      const int m1,
                      const int m2,
                      const int m3){
    if(m1+m2+m3) return 0.0;
    else return Wigner3(j1,j2,j3,m1,m2);
}
/*overload of Wigner3 to add m3 (since m3 is redundant)
 Output: Wigner 3-j symbol value
 Effect: none
 */

#ifdef W_JLG_map
inline double CG(const unsigned int J,
                 const int M,
                 const unsigned int j1,
                 const int m1,
                 const unsigned int j2,
                 const int m2,
                 const bool th=false){
    return parity((j1-j2+M)/2)*std::sqrt(J+1)*W3J(j1,j2,J,m1,m2,-M,th);
}
/*Calculates Clebsch-Gordan coefficients. Caution, all Js and Ms are doubled.
 Output: C-G coeficienct
 Effect: calls W3J, which creates a static/global map (if not already created)
 */
#endif //W_JLG_map

inline double Clebsch(const unsigned int J,
                      const int M,
                      const unsigned int j1,
                      const int m1,
                      const unsigned int j2,
                      const int m2){
    return parity((j2-j1-M)/2)*std::sqrt(J+1)*Wigner3(j1,j2,J,m1,m2,-M);
}
/*Calculateds Clebsch-Gordan coefficients. Caution, all Js and Ms are doubled.
 Similar to CG, but doesn't create static/global map
 Output: C-G coeficienct
 Effect: none
 */

/*------------------------------------------------------------------------------
 Wigner 6-j symbols
 / j1/2 j2/2 j3/2 \
 \ j4/2 j5/2 j6/3 /

 Clebsch-Gordan coeffidients
 W(a,b,c,d,e,f)=(-1)^((a+b+c+d)/2) * / a b e \
                                     \ d c f /
------------------------------------------------------------------------------*/
#ifdef W_JLG_map
double W6J(const unsigned int,//j1
           const unsigned int,//j2
           const unsigned int,//j3
           const unsigned int,//j4
           const unsigned int,//j5
           const unsigned int,//j6
           const bool=false);//th
/*Calculateds Wigner 6-j symbol. Caution, all Js are doubled.
 Output: Wigner 6-j symbol value
 Effect: (if not already created) creates static/global map, otherwise edits it
 */

inline double W6J(const key6 k,const bool th=false){
    return W6J(k[0],k[1],k[2],k[3],k[4],k[5],th);
}
/*overload to allow using key6 type as argument
 Output: Wigner 6-j symbol value
 Effect: (if not already created) creates static/global map
 */
#endif //W_JLG_map

double Wigner6(const unsigned int,//j1
               const unsigned int,//j2
               const unsigned int,//j3
               const unsigned int,//j4
               const unsigned int,//j5
               const unsigned int);//j6
/*Calculateds Wigner 6-j symbol. Caution, all Js are doubled.
 Similar to W6J, but does not create static/global map.
 Output: Wigner 6-j symbol value
 Effect: none
 */

#ifdef W_JLG_map
inline double RW(const unsigned int j1,
                 const unsigned int j2,
                 const unsigned int j3,
                 const unsigned int j4,
                 const unsigned int j5,
                 const unsigned int j6,
                 const bool th=false){
    return parity((j1+j2+j4+j3)/2,1.)*W6J(j1,j2,j5,
                                          j4,j3,j6,th);
}
/*Calculates Racah coefficients. Caution, all Js are doubled.
 Output: Racah coeficienct
 Effect: calls W6J, which creates or edits a static/global map (if not already created)
 */
#endif //W_JLG_map

inline double Racah(const unsigned int j1,
                    const unsigned int j2,
                    const unsigned int j3,
                    const unsigned int j4,
                    const unsigned int j5,
                    const unsigned int j6){
    return parity((j1+j2+j4+j5)/2,1.)*Wigner6(j1,j2,j5,
                                              j4,j3,j6);
}
/*CCalculates Racah coefficients. Caution, all Js and Ms are doubled.
 Similar to R, but doesn't create static/global map
 Output: Racah coeficienct
 Effect: none
 */

/*-----------------------------------------------------------------------------
 Wigner 9-j symbols
 / j1/2 j2/2 j3/2 \
 | j4/2 j5/2 j6/3 |
 \ j7/2 j8/2 j9/2 /
------------------------------------------------------------------------------*/
#ifdef W_JLG_map
double W9J(const unsigned int,//j1
           const unsigned int,//j2
           const unsigned int,//j3
           const unsigned int,//j4
           const unsigned int,//j5
           const unsigned int,//j6
           const unsigned int,//j7
           const unsigned int,//j8
           const unsigned int,//9
           const bool=false);//th
/*Calculateds Wigner 9-j symbol. Caution, all Js are doubled.
 Output: Wigner 9-j symbol value
 Effect: (if not already created) creates static/global map, otherwise edits it
 */

inline double W9J(const key9 k,const bool th=false){
    return W9J(k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],k[8],th);
}
/*overload to allow using key9 type as argument
 Output: Wigner 9-j symbol value
 Effect: calls W9J, which creates or edits a static/global map (if not already created)
 */
#endif

double Wigner9(const unsigned int,//j1
               const unsigned int,//j2
               const unsigned int,//j3
               const unsigned int,//j4
               const unsigned int,//j5
               const unsigned int,//j6
               const unsigned int,//j7
               const unsigned int,//j8
               const unsigned int);//9
/*Calculateds Wigner 9-j symbol. Caution, all Js are doubled.
 Similar to W9J, but does not create static/global map.
 Output: Wigner 9-j symbol value
 Effect: none
 */

/*------------------------------------------------------------------------------------------------
 Harmonic oscillator normalization
 N=sqrt( 2 Gamma(n+1)
         -------------
         Gamma(n+l+1.5) )
-----------------------------------------------------------------------------------------------*/
#ifdef W_JLG_map
Double HN_l(const unsigned int,//n
            const unsigned int);//l
/*Gives the logarithm of the normalization coefficeint for the harmonic
 oscillator in 3D. Stores values in a static map for quick recovery.
 Output: log of normalization coefficient
 Effect: (if not already created) creates static map of values
 */
inline Double HN(const unsigned int n,
                 const unsigned int l){
    return std::exp(HN_l(n,l));
}
/*Returns the normalization coefficient for the 3D HO.
 Output: normalization coefficient
 Effect: calls harm_norm_l
 */
inline Double HN_l(key2 k){
    return HN_l(k[0],k[1]);
}
inline Double HN(key2 k){
    return HN(k[0],k[1]);
}
#endif //W_JLG_map

Double harm_norm_l(const unsigned int,//n
                   const unsigned int);//l
/*Gives the logarithm of the normalization coefficeint for the harmonic
 oscillator in 3D.
 Output: log of normalization coefficient
 */
inline Double harm_norm(const unsigned int n,
                        const unsigned int l){
    return std::exp(harm_norm_l(n,l));
}
/*Returns the normalization coefficient for the 3D HO.
 Output: normalization coefficient
 Effect: calls harm_norm_l
 */

#endif /* defined(____algebra__) */
