//
//  funch.h
//  
//
//  Created by Johnathan Gross
//

#ifndef funch_h
#define funch_h

#ifdef D6_JLG
#define D3_JLG
#define _complex_JLG
#endif

#include "type.h"
#include <array>
#include <functional>

typedef std::array<Double,3> marray;

typedef std::function<Double(const Double,//x
                             const marray&)> func1;
//x is assumed to be the length, 0 to infinity
#ifdef D3_JLG
typedef std::function<Double(const Double,//rp
                             const Double,//rl
                             const Double,//tl
                             const marray&)> func3;
//rp and rl are assumed to be the lengths, 0 to infinity
//tl is assumed to be the cosine of the angle between rho and lambda, -1 to 1

inline func3 f1to3(const func1&f,
                   const bool t=true){
    return [&f,t](const Double rp,
                  const Double rl,
                  const Double,
                  const marray&m){
        return f(t?rp:rl,m);
    };
}

inline func1 f3to1(const func3&f,
                   const bool t=true){
    return [&f,t](const Double x,
                    const marray&m){
        return t?f(x,0,1,m):f(0,x,1,m);
    };
}

#ifdef D6_JLG
typedef std::function<CDouble(const Double,//rp
                              const Double,//rl
                              const Double,//tp
                              const Double,//tl
                              const Double,//pp
                              const Double,//pl
                              const marray&)> func6;
typedef std::function<Double(const Double,//rp
                             const Double,//rl
                             const Double,//tp
                             const Double,//tl
                             const Double,//pp
                             const Double,//pl
                             const marray&)> func6r;
//rp and rl are assumed to be of the lengths, 0 to infinity
//tp and tl are assumed to be the cosines of the angle, -1 to 1
//pp and pl are assumed to be the angle, 0 to 2pi

inline func6 f6rto6(const func6r&f){
    return [&f](const Double rp,
                const Double rl,
                const Double tp,
                const Double tl,
                const Double pp,
                const Double pl,
                const marray&m){
        return std::complex<Double>(f(rp,rl,tp,tl,pp,pl,m));
    };
}
inline func6r f6tor(const func6&f){
    return [&f](const Double rp,
                const Double rl,
                const Double tp,
                const Double tl,
                const Double pp,
                const Double pl,
                const marray&m){
        return f(rp,rl,tp,tl,pp,pl,m).real();
    };
}
inline func6r f6toi(const func6&f){
    return [&f](const Double rp,
                const Double rl,
                const Double tp,
                const Double tl,
                const Double pp,
                const Double pl,
                const marray&m){
        return f(rp,rl,tp,tl,pp,pl,m).imag();
    };
}

inline func6 f1to6(const func1&f,
                   const bool t=true){
    return [&f,t](const Double rp,
                  const Double rl,
                  const Double,
                  const Double,
                  const Double,
                  const Double,
                  const marray&m){
        return f(t?rp:rl,m);
    };
}

inline func6 f3to6(const func3&f){
    return [&f](const Double rp,
                const Double rl,
                const Double tp,
                const Double tl,
                const Double pp,
                const Double pl,
                const marray&m){
        return f(rp,rl,
                 std::sqrt((1.-tp*tp)*(1.-tl*tl))*std::cos(pp-pl)+tp*tl,
                 m);
    };
}

inline func3 f6to3(const func6&f){
    return [&f](const Double rp,
                const Double rl,
                const Double tl,
                const marray&m){
        return f(rp,rl,1,tl,0,0,m).real();
    };
}

inline func1 f6to1(const func6&f,
                   const bool t=true){
    return [&f,t](const Double x,
                  const marray&m){
        return (t?f(x,0,1,1,0,0,m):f(0,x,1,1,0,0,m)).real();
    };
}
#endif
#endif

#endif //funch_h
