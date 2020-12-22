//  Creates template functions for doing calculations with polynomials.
//	Hermite, Laguerre, and Legendre
//  polynomials.h
//
//	much has been invalidated by C++17 numerics library
//
//  Created by Johnathan Gross
//
//
/*polynomails:
 hermite
 laguerre
 legendre
 */

#ifndef ____polynomials_h
#define ____polynomials_h

#include "type.h"
#include <cmath>
#include <cstdlib>

//function to return x mod y with result in [0,|y|)
template<typename T,typename U> constexpr
typename std::common_type_t<T,U> mod(const T x,const U y){
    static_assert(std::is_arithmetic<T>::value,"T not arithmetic");
    static_assert(std::is_arithmetic<U>::value,"U not arithmetic");
    if constexpr(std::is_integral<T>::value and std::is_integral<U>::value){
        if constexpr(std::is_signed<T>::value){
            const std::common_type_t<T,U> z=x%y;
            return (0>x and bool(z))?z+std::abs(y):z;
        }else return x%y;
    }else{
        const std::common_type_t<T,U> z=fmod(x,y);
        return (0>x and bool(z))?z+std::abs(y):z;
    }
}

//function to return the square of a value
template<typename T> constexpr
T sqr(const T x){
    return x*x;
}

//function to return (-1)^n
template<typename I> constexpr
typename std::make_signed_t<I> parity(const I n){
    static_assert(std::is_integral<I>::value,"not integral type");
    //static_assert(std::is_signed<I>::value,"not signed type");
    return (n%2)?-1:1;
}
template<typename I, typename T=int> constexpr
T parity(const I n,const T k){
    static_assert(std::is_integral<I>::value,"not integral type");
    static_assert(std::is_arithmetic<T>::value,"not arithmetic type");
    static_assert(std::is_signed<T>::value,"not signed type");
    return (n%2)?-k:k;
}

//******************************************************************************
//Polynomial values given
template <typename T, class U>
T polycalc(const T,//x
           const std::vector<U>&,//poly
           const std::size_t o=-1);
/*calculates value of poly at point x, up to order o if provided
 Restrictions: T and U must be arithmetic types
 Output: value of polynomial
 Effect: none
 */

//******************************************************************************
//Hermite
template<class T>
void hermite(const unsigned int,//rank
             matrix<T>&);//herm
/*Creates a matix of coefficients for the Hermite polynomials
 Restrictions: T must be arithmetic type
 Effect: changes herm to coefficients
 Warning: may not be accurate for integral types
 */

template<typename T>
void hermitex(const unsigned int,//o
              const T,//x
              std::vector<T>&);//poly
/*caluclates Hermite polynomials up to order o at x
 Restrictions: T must be floating point type
 Effect: changes poly to values of polynomials
 */

template<typename T>
T hermitex(const unsigned int,//o
           const T);//x
/*caluclates Hermite polynomial of order o at x
 Restrictions: T must be floating point type
 overload of void hermitex<T>(const int,const T)
 Effect: none
 */

//******************************************************************************
//binomial
template<typename T,typename U>
Double binomial(const T,//n
                const U);//k
/*calculates general binomial coefficient for arbitrary n and k
 Restrictions: T must be arithmetic type
 Output: binomial coefficient
 Effect: none
 */

//******************************************************************************
//Laguerre
//C++17 special functions only work for integer a
template<typename T>
void laguerre(const unsigned int,//o
              const T,//a
              matrix<T>&);//Lag
/*calculates general Laguerre polynomials up to order o for given rank a
 Restrictions: T must be floating point type
 Effect: changes Lag to coefficients
 */
template<typename T>
inline void laguerre(const unsigned int ord,matrix<T>& Lag)
{laguerre(ord,T(0),Lag);}
/*overload for a=0
 */

template<typename T>
void laguerrex(const unsigned int,//o
               const T,//a
               const T,//x
               std::vector<T>&);//poly
/*calculates the Laguerre polynomials up to order o at x for given a
 Restrictions: T must be floating point type
 Effect: changes poly to values of polynomials
 */
template<typename T>
inline void laguerrex(const unsigned int o,
                      const T x,
                      std::vector<T>& poly)
{laguerrex(o,T(0),x,poly);}
/*overload for a=0
 */

template<typename T>
T laguerrex(const unsigned int,//o
            const T,//a
            const T);//x
/*calculates the Laguerre polynomial of order o at x
 Restrictions: T must be floating point type
 Effect: none
 */
template<typename T>
inline T laguerrex(const unsigned int o,const T x)
{return laguerrex(o,T(0),x);}
/*overload for a=0
 */

//******************************************************************************
//3D Harmonic Oscillator Normalization
template<typename T>
void harmoniclx(const unsigned int,//nmax
               const unsigned int,//lmax
               matrix<T>&);//Harm
/*calculates logarithm of Harmonic Oscillator normalization factor
 Restrictions: T must be floating point type
 Effect: changes Harm
 */
template<typename T>
inline void harmonicx(const unsigned int nmax,
               const unsigned int lmax,
               matrix<T>&Harm){
    harmoniclx(nmax,lmax,Harm);
#pragma omp parallel for collapse(2)
    for(std::size_t n=0;n<=nmax;++n) for(std::size_t l=0;l<=lmax;++l)
        Harm[n][l]=std::exp(Harm[n][l]);
}
/*gives value instead of logarithm
 Effect: changes Harm
 */

inline Double harmoniclx(const unsigned int n,
                         const unsigned int l){
    return (M_LN2+std::lgamma(n+Double(1))-std::lgamma(n+l+Double(1.5)))/2;
}
/* calculates Harmonic Oscillator normalziation factor for n and l
 Effect: none
 */
inline Double harmonicx(const unsigned int n,
                 const unsigned int l){
    return std::exp(harmoniclx(n,l));
}
/*gives value instead of logarithm
 Effect: none
 */

//******************************************************************************
//Normalized Associated Legendre
//defined as P~_l^m=sqrt((2l+1)(l-m)!/2/(l+m)!)P_l^m
//Note C++17 is the phi=0, m>=0 spherical harmonics, with theta in radians
//naglegendrex(l,m,x)=sqt(2pi)*sph_legendre(l,m,arccos(x))
template<typename T>
void nalegendrex(const unsigned int,//lmax
                 const T,//x
                 matrix<T>&);//poly
/*calculates the normalized associated Legendre polynomials up to order lmax at x for all 0<=m<=l
 Restrictions: T must be floating point type
 Effect: changes poly to values of polynomials
 */

template<typename T>
void nalegendrex(const unsigned int,//l
                 const T,//x
                 std::vector<T>&);//poly
/*calculates the normalized associated Legendre polynomials at x for all 0<=m<=l
 Restrictions: T must be floating point type
 Effect: changes poly to values of polynomials
 */

template<typename T>
T nalegendrex(const int,//l
             const int,//m
             const T);//x
/*calculates the normalized associated Legendre polynomial of order l at given m at x
 Restrictions: T must be floating point type
 Effect: none
 */

/*template<typename T>
 void nalegendrex(const unsigned int,//lmax
 const unsigned int,//m
 const T,//x
 std::vector<T>&);//poly*/
/*calculates the normalized associated Legendre polynomials up to order lmax at x for given m
 Effect: changes poly to values of polynomials

 Note: creates square matrix with m>l elements = 0
 */

//******************************************************************************
//Associated Legendre
/*template<typename T>
void alegendrex(const unsigned int,//lmax
                int,//m
                const T,//x
                std::vector<T>&);//poly*/
/*calculates the associated Legendre polynomials up to order lmax at x for given m
 Effect: changes poly to values of polynomials

 Note: creates square matrix with m>l elements = 0
 */

template<typename T>
void alegendrex(const unsigned int,//l
                const T,//x
                std::vector<T>&,//poly
                const bool=true);//b
/*calculates the associated Legendre polynomials at x with l
 for all 0<=m<=l for b=true and 0>=m>=-lamx for b=false
 Restrictions: T must be floating point type
 Effect: changes poly to values of polynomials
 */

template<typename T>
void alegendrex(const unsigned int,//lmax
                 const T,//x
                 matrix<T>&,//poly
                 const bool=true);//b
/*calculates the associated Legendre polynomials up to order lmax at x for all
 0<=m<=l for b=true and 0>=m>=-l for b=false
 Restrictions: T must be floating point type
 Effect: changes poly to values of polynomials
 */

template<typename T>
T alegendrex(const unsigned int,//l
             int,//m
             const T);//x
/*calculates the associated Legendre polynomial of order l at given m at x
 Restrictions: T must be floating point type
 Effect: none
 */

//******************************************************************************
//Legendre
template<typename T>
void legendre(const unsigned int,//l
              matrix<T>&);//P
/*calculates Legendre polynomials up to rank l
 Restrictions: T must be floating point type
 Effect: changes P to coefficients
 */

template<typename T>
inline void legendrex(const unsigned int,//lmax
                      const T,//x
                      std::vector<T>&);//poly
/*special case of the associated legendre polynomials for m=0
 Restrictions: T must be floating point type
 Effect: changes poly to values of polynomials
 */

template<typename T>
inline T legendrex(const unsigned int,//l
                   const T);//x
/*calculates Legendre polynomial P_l(x)
 Restrictions: T must be floating point type
 returns value of Legendre polynomial of order o at point x.
 */

#include "polynomials.cpp"
//includes definitions of templates
#endif //ifndef ____polynomials_h
