//
//  polynomials.cpp
//
//
//  Created by Johnathan Gross
//
//

#ifndef ____polynomials_t
#define ____polynomials_t

#ifdef TESTING_STATES_JLG
#include <stdexcept>
#endif

inline constexpr long double Mln2_4rootPI=0.40696470909759526588137528361991189016367643113142736124227424162535917949848819565995750942754626061649902796339;
inline constexpr long double Ms2PI=2.50662827463100050241576528481104525300698674060993831662992357634229365460784197494659583837805726611600997266520L;
inline constexpr long double M_s2=0.70710678118654752440084436210484903928483593768847403658833986899536623923105351942519376716382078636750692311545L;
inline constexpr long double Ms3_s2=1.22474487139158904909864203735294569598297374032833506421634628362548018872865751326992971655232011740929730061330L;
inline constexpr long double Ms3_2=0.86602540378443864676372317075293618347140262690519031402790348972596650845440001854057309337862428783781307070770L;

template <typename T,class U>
T polycalc(const T x, const std::vector<U>&poly,const std::size_t o){
    static_assert(std::is_arithmetic<T>::value,"x not arithmetic type");
    static_assert(std::is_arithmetic<U>::value,"poly not arithmetic type");
    /*calculates value of polynomial or order o at x*/
    T result=0;
    std::size_t m=std::min(poly.size(),o);
    for(std::size_t i=m-1;i<m;--i) result=x*result + poly[i];
    return result;
}


//******************************************************************************
//Hermite polynomials
template<typename T>
void hermite(const unsigned int rank, matrix<T>&Herm){
    static_assert(std::is_arithmetic<T>::value,"Not arithmetic: hermite");
    //floating point preferred, but not necessary
    Herm.resize(rank+1);
    for(size_t i=0;i<=rank;++i) Herm[i].assign(i+1,0);
    //recursion initialization
    Herm[0][0]=1;
    if(0==rank) return;
    //Herm[1][0]=0;
    Herm[1][1]=2;

    //recursive definition
    for(std::size_t i=2;i<=rank;++i){
        Herm[i][0]=-2*int(i-1)*Herm[i-2][0];
        for(std::size_t j=1;j<i-1;++j)
            Herm[i][j]=2*Herm[i-1][j-1] - 2*int(i-1)*Herm[i-2][j];
        Herm[i][i]=2*Herm[i-1][i-1];
    }
    return;
}

template<typename T>
void hermitex(const unsigned int o, const T x, std::vector<T>&poly){
    static_assert(std::is_floating_point<T>::value,"Not floating point: hermitex");
    poly.assign(o+1,0);
    poly[0]=1;
    if(0==o) return;
    poly[1]=x*2;
    for(size_t i=2;i<=o;++i) poly[i]=2* (x*poly[i-1] - int(i-1)*poly[i-2]);
    return;
}

template<>
inline float hermitex(const unsigned int o, const float x){
    return std::hermitef(o,x);
}
template<>
inline long double hermitex(const unsigned int o, const long double x){
    return std::hermitel(o,x);
}
template<typename T>
T hermitex(const unsigned int o, const T x){
    static_assert(std::is_floating_point<T>::value,"Not floating point: hermitex");
    return std::hermite(o,x);
}

//******************************************************************************
//Binomial coefficient
template <typename T,typename U>
Double binomial(const T n, const U k){
    static_assert(std::is_arithmetic<T>::value,"Not arithmetic type: binomial n");
    if constexpr(std::is_integral<U>::value){
#ifdef TESTING_STATES_JLG
        if(0>k) throw std::domain_error("bonimial: k negative");
#else
        if(0>k) return 0;
#endif
        if(0==k) return 1;
        if(1==k) return n;
        if constexpr(std::is_integral<T>::value){
            if(k>n) return 0;
            if(n==k) return 1;
            if(n-1==k) return n;
            if(2*k>n) return binomial(n,n-k);
        }
    }else static_assert(std::is_arithmetic<U>::value,"Not arithmetic type: binomial k");
    Double r=0;
#pragma omp critical(binomialgamma)
    r=std::exp(std::lgamma(n+1)-std::lgamma(k+1)-std::lgamma(n-k+1));
    if constexpr(std::is_integral<T>::value and std::is_integral<U>::value)
        return std::floor(.5+r);
    else return r;
}

//******************************************************************************
//Laguerre polynomials
template<typename T>
void laguerre(const unsigned int o,const T a, matrix<T>&Lag){
    static_assert(std::is_floating_point<T>::value,"Not floating point: laguerre");

    Lag.resize(o+1);
    Lag[0]={1};
    if(0==o) return;
    Lag[1]={1+a,-1};
    
    for(std::size_t i=2;i<=o;++i) Lag[i].assign(i+1,0);
#pragma omp parallel for
    for(std::size_t n=2;n<=o;++n)
        for(std::size_t i=0;i<=n;++i){
            Lag[n][i]=(i%2?-1:1) * binomial(n+a,n-i) / tgamma(i+T(1));
        }
}

template<typename T>
void laguerrex(const unsigned int o,const T a,const T x, std::vector<T>&poly){
    static_assert(std::is_floating_point<T>::value,"Not floating point: laguerrex");

    poly.resize(o+1);
    poly[0]=1;
    if(0==o) return;
    poly[1]=1+a-x;

    for(std::size_t k=1;o>k;++k)
        poly[k+1]=((2*k+1+a-x)*poly[k]-(k+a)*poly[k-1])/(k+1);
    return;
}

template<typename T>
T laguerrex(const unsigned int o,const T a,const T x){
    static_assert(std::is_floating_point<T>::value,"Not floating point: laguerrex");

    if(0==o) return 1;
    if(1==o) return 1+a-x;

    T result=0;
    for(unsigned int i=o;i>0;i--){
        result*=-x/(i+1.);
        result+=binomial(o+a,o-i);
    }
    result*=-x;
    result+=binomial(o+a,int(o));
    return result;
}

//******************************************************************************
//Harmonic Osciallator Normalizatoin
template<typename T>
void harmoniclx(const unsigned int nmax,const unsigned int lmax,matrix<T>&Harm){
    static_assert(std::is_floating_point<T>::value,"Not floating point: harmoniclx");

    Harm.assign(nmax+1,std::vector<T>(lmax+1,0));
    Harm[0][0]=Mln2_4rootPI;

    for(unsigned int l=1;l<=lmax;++l)
        Harm[0][l]=Harm[0][l-1]-std::log(l+Double(0.5))/2;

#pragma omp parallel for ordered
    for(unsigned int n=1;n<=nmax;++n){
#pragma omp ordered
        Harm[n][0]=Harm[n-1][0]+(std::log(Double(n))-std::log(n+Double(0.5)))/2;
        for(unsigned int l=1;l<=lmax;++l){
            Harm[n][l]=Harm[n][l-1]-std::log(n+l+Double(0.5))/2;
        }
    }
}

/*Double harmoniclx(const unsigned int n,const unsigned int l){
    //if(0==n && 0==l) return (M_LN2+std::log(M_2_SQRTPI))/2;
    //else if (0==l) return harmoniclx(n-1,0)+(std::log(n)-std::log(n+0.5))/2;
    //else return harmoniclx(n,l-1)-std::log(n+l+0.5)/2;
    return (M_LN2+std::lgamma(n+1)-std::lgamma(n+l+1.5))/2;
}*/

//******************************************************************************
//Normalized Associated Legendre polynomials
/*template<typename T>
void nalegendrex(const unsigned int lmax, const unsigned int m,const T x,
                 std::vector<T>&poly){
    //up to lmax, given m
    static_assert(std::is_floating_point<T>::value,"Not floating point: nalegendrex");
#ifdef TESTING_STATES_JLG
    if(m>lmax) throw std::invalid_argument("nalegendrex: m>lmax");
#endif

    poly.assign(lmax+1,0);
    poly[m]=(T(m)-.5)*M_LN2+std::lgamma(T(m)+.5)
    +(std::log(2*T(m)+1)-std::lgamma(2*T(m)+1)-std::log(T(M_PI))
      +T(m)*std::log((1-x)*(1+x)))/2;

    if(lmax>m) poly[m+1]=poly[m]+std::log(Double(2*m+3))/2.;
    poly[m]=(m%2?-1:1)*std::exp(poly[m]);
    if(lmax==m) return;
    poly[m+1]=x*(m%2?-1.:1.)*std::exp(poly[m+1]);
    if(m+1==lmax) return;

    T of=sqrt(2.*m+3.);
    for(unsigned int l=m+2;l<=lmax;++l){
        T f=std::sqrt(Double((4*l*l-1)/(l*l-m*m)));
        poly[l]=(x*poly[l-1]-poly[l-2]/of)*f;
        of=f;
    }
    return;
}

template<typename T>
void nalegendrex(const unsigned int l,const T x, std::vector<T>&poly){
    static_assert(std::is_floating_point<T>::value,"Not floating point: nalegendrex");

    poly.assign(l+1,0);

    matrix<T> p;
    nalegendrex(l,x,p);
#pragma omp parallel for
    for(int m=0;m<=l;++m) poly[m]=p[l][m];
    return;
}

template<typename T>
void nalegendrex(const unsigned int lmax,const T x,matrix<T>&poly){
    static_assert(std::is_floating_point<T>::value,"Not floating point: nalegendrex");

    poly.resize(lmax+1);
#pragma omp parallel for
    for(std::size_t l=0;l<=lmax;++l) poly[l].resize(l+1,0);
#pragma omp parallel for
    for(std::size_t m=0;m<=lmax;++m){
        std::vector<T> p;
        nalegendrex(lmax,m,x,p);
        for(std::size_t l=m;l<=lmax;++l) poly[l][m]=p[l];
    }
    return;
}*/

template <typename T>
void nalegendrex(const unsigned int lmax,const T x,matrix<T>&poly) {
    static_assert(std::is_floating_point<T>::value,"Not floating point: nalegendrex");
    poly.resize(lmax+1);
    for(size_t l=0;l<=lmax;l++) poly[l].assign(l+1,0);
    if(-1>x or 1<x){
#ifdef TESTING_STATES_JLG
        throw(std::domain_error("nalegendrex: |x|>1"));
#else
        return;
#endif
    }

    poly[0][0]=M_s2;
    if(0==lmax) return;
    const Double s=std::sqrt(1-x*x);
    const bool b=((1-std::abs(x))>.01);//true if x not close to +-1
    poly[1][0]=Ms3_s2*x;
    poly[1][1]=-Ms3_2*s;
    if(1==lmax) return;

    Double t=std::acos(std::abs(x));
#pragma omp parallel for ordered
    for(size_t l=2;l<=lmax;l++){
#pragma omp ordered
        poly[l][0]=(std::sqrt(Double(2*l-1)*(2*l+1))*x*poly[l-1][0]
            -std::sqrt((2*l+1)/Double(2*l-3))*(l-1)*poly[l-2][0])/l;
        if(b)
            poly[l][1]=std::sqrt(l/Double(l+1))/s
            *(x*poly[l][0]-std::sqrt(Double(2*l+1)/(2*l-1))*poly[l-1][0]);
        else poly[l][1]=std::sqrt(Double(2*l+1)*(l+1)*l/8)
            *(1-Double(3*l*(l+1)-2)/24*t*t)*t*parity(1+(l-1)*(0>x));
        for(size_t m=2;m<=l;m++){
            if(b)
                poly[l][m]=-2*x*(m-1)/s/std::sqrt((l*Double(l+1)-(m-1)*m))*poly[l][m-1]
                -std::sqrt(Double(l*(l+1)-(m-1)*(m-2))/(l*(l+1)-(m-1)*m))*poly[l][m-2];
            else poly[l][m]=std::sqrt(Double(2*l+1)/2)
                *std::exp((std::lgamma(l+m+1)-std::lgamma(l-m+1))/2-std::lgamma(m+1))
                *(1-Double(3*l*(l+1)-m*(m+1))/(12*(m+1))*t*t)
                *std::pow(t/2,m)*parity(m+(l-m)*(0>x));
        }
    }
}

template <typename T>
void nalegendrex(const unsigned int l,const T x,std::vector<T>&poly) {
    static_assert(std::is_floating_point<T>::value,"Not floating point: nalegendrex");
    poly.assign(l+1,0);
    if(-1>x or 1<x){
#ifdef TESTING_STATES_JLG
        throw(std::domain_error("nalegendrex: |x|>1"));
#else
        return;
#endif
    }

    poly[0]=nalegendrex(l,0,x);
    if(0==l) return;
    poly[1]=nalegendrex(l,1,x);
    if(1==l) return;
    const Double s=std::sqrt(1-x*x);
    const Double t=std::acos(std::abs(x));
    const bool b=((1-std::abs(x))>.01);//true if x not close to +-1
    for(size_t m=2;m<=l;m++){
        if(b)
            poly[m]=-2*x*(m-1)/s/std::sqrt((l*Double(l+1)-(m-1)*m))*poly[m-1]
            -std::sqrt(Double(l*(l+1)-(m-1)*(m-2))/(l*(l+1)-(m-1)*m))*poly[m-2];
        else poly[m]=std::sqrt(Double(2*l+1)/2)
            *std::exp((std::lgamma(l+m+1)-std::lgamma(l-m+1))/2-std::lgamma(m+1))
            *(1-Double(3*l*(l+1)-m*(m+1))/(12*(m+1))*t*t)
            *std::pow(t/2,m)*parity(m+(l+m)*(0>x));
    }
}

template<typename T>
T nalegendrex(const int l,const int m,const T x){
    static_assert(std::is_floating_point<T>::value,"Not floating point: nalegendrex");
    if(-1>x or 1<x){
#ifdef TESTING_STATES_JLG
     throw(std::domain_error("nalegendrex: |x|>1"));
#else
        return 0;
#endif
    }
    if(0>l and 0>m) return nalegendrex(-l-1,-m,x);
    if(0>l) return (m%2?-1:1)*nalegendrex(-l-1,m,x);
    if(0>m) return (m%2?-1:1)*nalegendrex(l,-m,x);

    return Ms2PI*std::sph_legendre(l,m,std::acos(x));
}

//******************************************************************************
//Associated Legendre polynomials
/*template<typename T>
void alegendrex(const unsigned int lmax,int m,const T x,
                std::vector<T>&poly){
    static_assert(std::is_floating_point<T>::value,"Not floating point: alegendrex");
    int s=1;
    if(0>m){
        m=-m;
        s=-1;
    }
    nalegendrex(lmax,m,x,poly);
    if(0==m){
#pragma omp parallel for
        for(int l=0;l<=lmax;++l) poly[l]*=std::sqrt(Double(2)/(2*l+1));
    }else{
#pragma omp parallel for
        for(int l=m;l<=lmax;++l){
            poly[l]*=std::exp((std::lgamma(Double(l+s*m+1))
                               -std::lgamma(Double(l-s*m+1))
                               +M_LN2-std::log(Double(2*l+1)))/2)
            *((-1==s && m%2)?-1:1);
        }
    }
    return;
}*/

template<typename T>
void alegendrex(const unsigned int l,const T x, std::vector<T>&poly,const bool b){
    static_assert(std::is_floating_point<T>::value,"Not floating point: alegendrex");

    nalegendrex(l,x,poly);
    int s=(b?1:-1);
#pragma omp parallel for
    for(size_t m=0;m<=l;++m){
        poly[m]*=std::exp((std::lgamma(Double(l+s*m+1))
                           -std::lgamma(Double(l-s*m+1))
                           +M_LN2-std::log(Double(2*l+1)))/2)
        *((!b && m%2)?-1:1);
    }
    return;
}

template<typename T>
void alegendrex(const unsigned int lmax,const T x,
                matrix<T>&poly,const bool b){
    static_assert(std::is_floating_point<T>::value,"Not floating point: alegendrex");

    nalegendrex(lmax,x,poly);
    int s=(b?1:-1);
#pragma omp parallel for
    for(size_t l=0;l<=lmax;++l) for(size_t m=0;m<=l;++m){
        poly[l][m]*=std::exp((std::lgamma(Double(l+s*m+1))
                              -std::lgamma(Double(l-s*m+1))
                              +M_LN2-std::log(Double(2*l+1)))/2)
        *((!b && m%2)?-1:1);
    }
    return;
}

template<typename T>
T alegendrex(const unsigned int l,int m,const T x){
    static_assert(std::is_floating_point<T>::value,"Not floating point: alegendrex");

    if(0==m) return std::legendre(l,x);
    if(0>m) return std::assoc_legendre(l,-m,x)
        *std::exp(std::lgamma(l+m+1)-std::lgamma(l-m+1));
    return std::assoc_legendre(l,m,x)*(m%2?-1:1);
}

//******************************************************************************
//Legendre polynomials
template<typename T>
void legendre(const unsigned int l, matrix<T>&P){
    static_assert(std::is_floating_point<T>::value,"Not floating point: legendre");
    P.resize(l+1);
    for(size_t i=0;i<=l;++i) P[i].assign(i+1,0);
    P[0][0]=1;
    if(0==l) return;
    P[1][1]=1;

    for(size_t n=2;n<=l;++n){
        if(!(n%2)) P[n][0]=-P[n-2][0]*(n-1)/n;
        for(size_t i=2-(n%2);i<=n-2;i+=2){
            P[n][i]=(P[n-1][i-1]*(2*n-1)-P[n-2][i]*(n-1))/n;
        }
        P[n][n]=P[n-1][n-1]*(2*n-1)/n;
    }
    return;
}

template<typename T>
inline void legendrex(const unsigned int lmax,
                      const T x,
                      std::vector<T>&poly){
    static_assert(std::is_floating_point<T>::value,"Not floating point: legendrex");
    Dmatrix P;
    nalegendrex(lmax,x,P);
    poly.resize(lmax+1);
#pragma omp parallel for
    for(size_t l=0;l<=lmax;++l) poly[l]=std::sqrt(Double(2)/(2*l+1))*P[l][0];
}

template<typename T>
inline T legendrex(const unsigned int l,const T x){
    static_assert(std::is_floating_point<T>::value,"Not floating point: legendrex");
    return std::legendre(l,x);
}
#endif //#ifndef ____polynomials_t
