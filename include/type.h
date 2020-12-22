//  Defines numerical types
//  type.h
//  
//
//  Created by Johnathan Gross on 5/20/15.
//
//

#ifndef _type_h
#define _type_h

#include <vector>
#ifdef _complex_JLG
#include <complex>
#endif

//allows making a matrix of values of general type without vector<vector>
template <typename T>
using matrix = std::vector<std::vector<T>>;

//typedefs to quickly create vectors and matrices (and complex) of specific types
typedef std::vector<bool> bvector;
typedef matrix<bool> bmatrix;
typedef std::vector<int> ivector;
typedef matrix<int> imatrix;
typedef std::vector<unsigned int> uivector;
typedef matrix<unsigned int> uimatrix;
typedef std::vector<long int> livector;
typedef matrix<long int> limatrix;
typedef std::vector<float> fvector;
typedef matrix<float> fmatrix;
typedef std::vector<double> dvector;
typedef matrix<double> dmatrix;
typedef std::vector<long double> ldvector;
typedef matrix<long double> ldmatrix;
typedef std::vector<std::size_t> tvector;
typedef matrix<std::size_t> tmatrix;
#ifdef _complex_JLG
typedef std::complex<float> cfloat;
typedef std::vector<cfloat> cfvector;
typedef matrix<cfloat> cfmatrix;
typedef std::complex<double> cdouble;
typedef std::vector<cdouble> cdvector;
typedef matrix<cdouble> cdmatrix;
typedef std::complex<long double> cldouble;
typedef std::vector<cldouble> cldvector;
typedef matrix<cldouble> cldmatrix;
#endif //_complex_JLG

//allows using long double instead of double without altering code
#ifdef _long_JLG
typedef long double Double;
#else
typedef double Double;
#endif
typedef std::vector<Double> Dvector;
typedef std::vector<Dvector> Dmatrix;
#ifdef _complex_JLG
typedef std::complex<Double> CDouble;
typedef std::vector<CDouble> CDvector;
typedef std::vector<CDvector> CDmatrix;
#endif //_complex_JLG

//function to return sign of x
template <typename T> constexpr
int sign(const T x, std::false_type is_signed [[maybe_unused]]) {
    //return T(0)<x;
    return T(0)!=x;
}
template <typename T> constexpr
int sign(const T x, std::true_type is_signed [[maybe_unused]]) {
    //return int(T(0)<x)-int(T(0)>x);
    return T(0)<x?1:T(0)>x?-1:0;
}
template <typename T> constexpr
int sign(const T x) {
    return sign(x,std::is_signed<T>());
}

#ifdef _complex_JLG
template<typename T> constexpr
bool Cbool(const std::complex<T> x){
    return bool(x.real()) or bool(x.imag());
}

template<typename T> inline
std::vector<T> realvec(const std::vector<std::complex<T>>&c){
    std::vector<T> r(c.size());
#pragma omp parallel for
    for(std::size_t j=0;j<c.size();j++) r[j]=c[j].real();
    return r;
}
template<typename T> inline
std::vector<T> imagvec(const std::vector<std::complex<T>>&c){
    std::vector<T> i(c.size());
#pragma omp parallel for
    for(std::size_t j=0;j<c.size();j++) i[j]=c[j].imag();
    return i;
}
template<typename T> inline
matrix<T> realmat(const matrix<std::complex<T>>&c){
    matrix<T> r(c.size(),std::vector<T>(c[0].size()));
    for(std::size_t j=0;j<c.size();j++)
        r[j].resize(c[j].size());
#pragma omp parallel for
    for(std::size_t j=0;j<c.size();j++) for(std::size_t k=0;k<c[j].size();k++)
        r[j][k]=c[j][k].real();
    return r;
}
template<typename T> inline
matrix<T> imagmat(const matrix<std::complex<T>>&c){
    matrix<T> r(c.size(),std::vector<T>(c[0].size()));
    for(std::size_t j=0;j<c.size();j++) r[j].resize(c[j].size());
#pragma omp parallel for
    for(std::size_t j=0;j<c.size();j++) for(std::size_t k=0;k<c[j].size();k++)
        r[j][k]=c[j][k].imag();
    return r;
}
#endif

#endif //_type_h
