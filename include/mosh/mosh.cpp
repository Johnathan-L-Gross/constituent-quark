//
//  Moshinsky.cpp
//
//
//  Created by Johnathan Gross
//
//

#define W_JLG_map
#include "algebra/algebra.h"
#include "mosh.h"
#include <cmath>
#include <cstdlib>
#include <algorithm>
#ifdef TESTING_STATES_JLG
#include <stdexcept>
#endif

//inline const Double Ms3_2=0.866025403784438646763723170752936183471402626905190314027903L;
//already defined in polynomials.t
inline constexpr long double Ms2_spi=
0.79788456080286535587989211986876373695171726232986931533185165934131585179860367700250466781461387286060511772527L;

//inline int abs(int i){return i<0?-i:i;}
constexpr Double par(const int i){return i%2?-1.:1.;}

//static std::vector<Double> FA(0),DFA(0);
//FA(n)=FA[n]=n!=Gamma(n+1)
//DFA(n)=DFA[n+1]=n!!=Gamma(n/2+1)*2^(n/2)*(Pi/2)^(((-1)^n-1)/4)
//=Gamma(n/2+1)*2^(n/2)*(n%2?1:sqrt(2/Pi))
//=Gamma(n/2+1)*(n%2?1:{2/sqrt(PI)}/sqrt(2))

inline Double fa(int n){return std::tgamma(Double(n+1));}
//inline Double lfa(int n){return std::lgamma(n+1);}
inline Double dfa(int n){
    return std::pow(2,Double(n)/2)*std::tgamma(Double(n)/2+1)*(n%2?1:Ms2_spi);
}
//gamma(2n+1)2^(n/2)(n%2?1:sqrt(2/pi))
/*inline Double ldfa(int n){return n*M_LOG2E/2.+std::lgamma((n+1)/2.)
    -std::log(M_PI);
}*/

/*inline Double DD(const int I,const int J,const int K){
    //THE FUNCTION SUBPROGRAM OF DELTA
    if (I+J<K || abs(I-J)>K) return 0.;

    return std::sqrt(fa(I+J-K)*fa(J+K-I)/fa(I+J+K+1)*fa(K+I-J));
}*/

inline double CG(const unsigned int i,
                 const unsigned int j,
                 const unsigned int k){
    return CG(2*k,0,2*i,0,2*j,0,false);
}
//CG(I,J,K)(Gan et al.)==CG(2K,0,2I,0,2J,0)(mine)
/*[[gnu::const]]
inline Double CG(const int I,const int J,const int K){
    //THE FUNCTION SUBPROGRAM OF TNREE-ZERO-CG COEFFICIENT
    if((I+J+K)%2) return 0.;
    if(I+J<K || abs(I-J)>K) return 0.;
    const int M=(I+K+J)/2;

    return par(M-I+J)*std::sqrt(2*K+1)
    /fa(M-I)*DD(I,J,K)/fa(M-J)*fa(M)/fa(M-K);
}*/

inline double WIG(const unsigned int a,
                  const unsigned int b,
                  const unsigned int c,
                  const unsigned int d,
                  const unsigned int e,
                  const unsigned int f){
    return RW(2*a,2*b,2*c,2*d,2*e,2*f);
}
//WIG(A,B,C,D,E,F)==RW(2A,2B,2C,2D,2E,2F)
/*[[gnu::const]]
inline Double WIG(const int IA,const int IB,const int IC,
           const int ID,const int IE,const int JF){
    //THE FUNCTION SUBPROGRAM OF RACAH COEFFICIENT
    const Double D=DD(IA,IB,IE)*DD(IC,ID,IE)*DD(IA,IC,JF)*DD(IB,ID,JF);
    if(0.==D) return 0.;
    const int N1=std::max({IA+IB+IE,IC+ID+IE,IA+IC+JF,IB+ID+JF});
    const int N2=std::min({IA+IB+IC+ID,IA+ID+IE+JF,IB+IC+IE+JF});
    if(N2<N1) return 0.;
    Double SUM=0.;
#pragma omp parallel for reduction(+:SUM)
    for(int K=N1;K<=N2;K++){
        SUM+=par(K)*fa(K+1)
        /(fa(K-IA-IB-IE)*fa(K-IC-ID-IE)*fa(K-IA-IC-JF)*fa(K-IB-ID-JF))
        /(fa(IA+IB+IC+ID-K)*fa(IA+ID+IE+JF-K)*fa(IB+IC+IE+JF-K));
    }
    return SUM*par(IA+IB+IC+ID)*D;
}*/

inline double U(const unsigned int A,
                const unsigned int B,
                const unsigned int C,
                const unsigned int D,
                const unsigned int E,
                const unsigned int EP,
                const unsigned int F,
                const unsigned int FP,
                const unsigned int G){
    return W9J(2*A,2*B,2*E,2*C,2*D,2*EP,2*F,2*FP,2*G,false);
}
//U(A,B,C,D,E,EP,F,FP,G)==W9J(2A,2B,2E,2C,2D,2EP,2F,2FP,2G)
/*[[gnu::const]]
inline Double U(const int KA,const int KB,const int KC,const int KD,const int KE,int KEP,const int KF,const int KFP,const int KG){
    //THE FUNCTION SUBPROFRAM OF U COEFFICIENT
    if((KA+KB)<KE || abs(KA-KB)>KE) return 0.;
    if((KC+KD)<KEP || abs(KC-KD)>KEP) return 0.;
    const int KSIG=KA+KB+KC+KD+KE+KF+KG+KEP+KFP;
    const int N1=std::max({abs(KB-KC),abs(KE-KF),abs(KEP-KFP)});
    const int N2=std::max({KB+KC,KE+KF,KEP+KFP});
    if(N2<N1) return 0.;
    Double SUM=0.;
#pragma omp parallel for reduction(+:SUM)
    for(int I=N1;I<=N2;I++){
        const Double W1=WIG(KB,KC,KE,KF,I,KA);
        if(0.==W1) continue;
        const Double W2=WIG(KB,KC,KFP,KEP,I,KD);
        if(0.==W2) continue;
        SUM+=(2*I+1)*W1*W2*WIG(KE,KF,KEP,KFP,I,KG);
    }
    return par(KSIG)*SUM;
}*/

inline Double H(const int L1,const int L2,
                const int LA1,const int LA2,
                const int L,const int K,
                const int LP1,const int LP2){
    //THE FUNCTION SUBPROGRAN OF H
    if(LP1>L1){
#ifdef TESTING_STATES_JLG
        throw std::invalid_argument("l3 must be <= l1");
#else
        return 0;
#endif
    }
    if(LP2>L2){
#ifdef TESTING_STATES_JLG
        throw std::invalid_argument("l4 must be <= l2");
#else
        return 0;
#endif
    }
    const Double cl=std::sqrt(Double((2*L1+1)*(2*L2+1)));
    Double h=0;
    const int IX=std::abs(LP1-LP2);
    const int IY=LP1+LP2;
    const int IX1=std::abs((L1-LP1)-(L2-LP2));
    const int IY1=(L1-LP1)+(L2-LP2);
#pragma omp parallel for reduction(+:h)
    for(int I1=IX;I1<=IY;I1++){
        const Double CGG=par(I1)*CG(LP1,LP2,I1)*CG(I1,K,LA1)
        *std::sqrt(Double(2*I1+1));
        if(0==CGG) continue;
        Double hsum=0;
        for(int I2=IX1;I2<=IY1;I2++){
            const Double UU=U(LP1,L1-LP1,LP2,L2-LP2,L1,L2,I1,I2,L);
            if(0==UU) continue;
            hsum+=std::sqrt(Double(2*I2+1))*UU
            *CG(L1-LP1,L2-LP2,I2)*CG(I2,K,LA2)*WIG(I2,I1,LA2,LA1,L,K);
        }
        h+=CGG*hsum;
    }
    return h*cl;
}

inline Double F(const int L1,const int L2,
                const int LA1,const int LA2,
                const int L,const int IQ,
                const int LP1,const int LP2){
    //THE FUNCTION SUBPROGRAM OF SMALL F
    if(LP1>L1){
#ifdef TESTING_STATES_JLG
        throw std::invalid_argument("l3 must be <= l1");
#else
        return 0;
#endif
    }
    if(LP2>L2){
#ifdef TESTING_STATES_JLG
        throw std::invalid_argument("l4 must be <= l4");
#else
        return 0;
#endif
    }
    Double f=0.;
#pragma omp parallel for reduction(+:f)
    for(int K=IQ;0<=K;K-=2){
        //std::cerr<<','<<KK<<' ';
        f+=(2*K+1)*dfa(IQ-1)/dfa(IQ-K)*dfa(IQ)/dfa(IQ+K+1)
        *H(L1,L2,LA1,LA2,L,K,LP1,LP2);
        //std::cerr<<f<<'\n';
    }
    return f;
}

inline Double FF(const int L1,
                 const int L2,
                 const int M1,
                 const int LA1,
                 const int LA2,
                 const int L,
                 const int IQ,
                 const int IS1,
                 const Double A,
                 const Double B){
    //THE FUNCTION SUBPROGRAN OF BIG F
    Double ff=0.;
#pragma omp parallel for reduction(+:ff)
    for(int I=0;I<=L1;I++){// for(int J=0;J<=L2;J++){
        int J=2*(M1-IS1)+LA1-IQ-I;
        //std::cerr<<' '<<I<<','<<J<<'\n';
        if(0>J or J>L2) continue;
        ff+=par(I)*std::pow(A,IQ+I-J+L2)*std::pow(B,IQ+J-I+L1)
        *F(L1,L2,LA1,LA2,L,IQ,I,J)
        *std::sqrt(fa(2*L1)/fa(2*I)/fa(2*(L1-I))*fa(2*L2)/fa(2*J)/fa(2*(L2-J)));
        //std::cerr<<'='<<ff<<'\n';
    }
    return ff;
}

inline Double R(const int N2,
                const int IS1,
                const int IS2,
                const int IQ,
                const Double A,
                const Double B){
    //THE FUNCTION SUBPROGRAM OF R
    Double r=0;
    const Double X=std::pow(Double(2),IQ);//IQ*M_LN2;
#pragma omp parallel for reduction(+:r)
    for(int I1=0;I1<=IS1;I1++){
        const Double X1=fa(IS1-I1)*fa(I1);
        Double r1=0;
        for(int I2=0;I2<=IS2;I2++){
            const Double X2=fa(IS2-I2)*fa(I2);
            const int I3=N2-I1-I2;
            if(0>I3 or IQ<I3) continue;
            //Double r2=0;
            //for(int I3=0;I3<=IQ;I3++){
                //if((I1+I2+I3)!=N2) continue;
                //const Double X3=fa(IQ-I3)*fa(I3);
                //r2+=par(I3)/X3;
            //}
            r1+=par(I3)/(fa(IQ-I3)*fa(I3))
            *std::pow(A,2*(IS1-I1+I2))*std::pow(B,2*(IS2+I1-I2))/X2;
        }
        r+=r1/X1;
    }
    return r*X;
}

Double tm3(const int N1,const int L1,
           const int N2,const int L2,
           const int M1,const int LA1,
           const int M2,const int LA2,
           const int L,
           const Double X,const Double Y,const Double Z,
           const int K){
    //THE FUNCTION SUBPROGRAM OF THREE-BODY TM COEFFICIENT
    if(0>N1 or 0>L1 or 0>N2 or 0>L2 or 0>M1 or 0>LA1 or 0>M2 or 0>LA2 or
       0>L or 0>X or 0>Y or 0>Z){
#ifdef TESTING_STATES_JLG
        throw std::invalid_argument("nonnegative arguments");
#else
    return 0;
#endif
    }
    if(1!=K and 2!=K and 3!=K){
#ifdef TESTING_STATES_JLG
    throw std::invalid_argument("k must be 1, 2, or 3");
#else
    return 0.;
#endif
    }
    if((2*(N1+N2)+L1+L2)!=(2*(M1+M2)+LA1+LA2)){
#ifdef TESTING_STATES_JLG
    throw std::invalid_argument("not the same HO level");
#else
    return 0.;
#endif
    }
    if(L1+L2<L or LA1+LA2<L or std::abs(L1-L2)>L or std::abs(LA1-LA2)>L){
#ifdef TESTING_STATES_JLG
        throw std::invalid_argument("triangle inequality");
#else
        return 0.;
#endif
    }

    Double A=0.,B=0.;
    switch (K) {
        case 1:
            A=(-1)*std::sqrt(X*Z/((X+Y)*(Y+Z)));//error? X*Z -> (X+Z)
            B=std::sqrt(Y*(X+Y+Z)/((X+Y)*(Y+Z)));
            break;
        case 2:
            A=(-1)*std::sqrt(Y*Z/((X+Y)*(X+Z)));//error? Y*Z -> (Y+Z)
            B=(-1)*std::sqrt(X*(X+Y+Z)/((X+Y)*(X+Z)));
            break;
        default:
            //K is 3
            return (N1==N2)*(N2==M2)*(L1==LA1)*(L2==LA2);
    }
    //std::cerr<<A<<' '<<B<<' ';

    const Double COEF=par(N1+N2-M1-M2+L1+LA2+L)
    *std::sqrt((2*L1+1)*(2*L2+1)
               *dfa(2*M1+2*LA1+1)/dfa(2*N1+2*L1+1)
               *dfa(2*M2+2*LA2+1)/dfa(2*N2+2*L2+1)
               *dfa(2*M1)/dfa(2*N1)
               *dfa(2*M2)/dfa(2*N2))
    *fa(N1)*fa(N2);
    //std::cerr<<COEF<<' ';
    Double SUM=0.;
    int NN1=N1+N2;
#pragma omp parallel for reduction(+:SUM)
    for(int I=0;I<=NN1;I++) for(int M=0;M<=NN1-I;M++){
        //std::cerr<<I<<','<<M<<'\n';
        int N=NN1-I-M;
        const Double R8=R(N2,M,N,I,A,B);
        const Double F8=FF(L1,L2,M1,LA1,LA2,L,I,M,A,B);
        SUM+=par(I)*R8*F8;
        //std::cerr<<R8<<' '<<F8<<' '<<SUM<<'\n';
    }
    //std::cerr<<SUM<<'\n';
    return COEF*SUM;
}

#ifdef TM4_JLG
Double tm4(const int N1,const int L1,
           const int N2,const int L2,
           const int N3,const int L3,
           const int M1,const int LA1,
           const int M2,const int LA2,
           const int M3,const int LA3,
           const int L0,const int L0T,const int L,
           const Double X1,const Double X2,const Double X3,const Double X4,
           const int K){
    if(0>N1 or 0>L1 or 0>N2 or 0>L2 or 0>N3 or 0>L3 or
       0>M1 or 0>LA1 or 0>M2 or 0>LA2 or 0>M3 or 0>LA3 or
       0>L0 or 0>L0T or 0>L or 0>X1 or 0>X2 or 0>X3 or 0>X4){
        //THE SUBPROGRAM OF 4-BODY TM COEFFICIENT
#ifdef TESTING_STATES_JLG
        throw std::invalid_argument("nonnegative arguments");
#else
        return 0;
#endif
    }
    auto W1=[&](const int K1,const int K2,const int K3){
        return WIG(L1,K1,L,K2,K3,L0);};
    auto W2=[&](const int K1,const int K2,const int K3){
        return WIG(LA1,K1,L,K2,K3,L0T);};
    auto A1=[&](const Double B1,const Double B2,const Double B3,const int KK){
        return tm3(N2,L2,N3,L3,M2,LA2,M3,LA3,L0,B1,B2,B3,KK);};
    auto A2=[&](const int K1,const int K2,const int KK){
        return tm3(N1,L1,K1,K2,M2,LA2,M3,LA3,L0T,X1,X2,X3+X4,KK);};
    auto A3=[&](const int K1,const int K2,const int I,
                const Double Y3,const int KK){
        return tm3(N1,L1,K1,K2,M1,LA1,M2,LA2,I,X1,X2,Y3,KK);};
    auto A4=[&](const int K1,const int K2,
                const Double B1,const Double B2,const Double B3){
        return tm3(N2,L2,N3,L3,K1,K2,M3,LA3,L0,B1,B2,B3,2);};
    auto A5=[&](const int K1,const int K2,const int K3,const int KK){
        return tm3(N1,L1,N2,L2,M1,LA1,K1,K2,K3,X1,X2,X3,KK);};
    auto A6=[&](const int K1,const int K2,
                const Double B1,const Double B2,const int KK){
        return tm3(K1,K2,N3,L3,M2,LA2,M3,LA3,L0T,B1,B2,X4,KK);};
    auto A7=[&](const int K1,const int K2,const int K3,const int K4,const int KK){
        return tm3(N2,L2,N3,L3,K1,K2,K3,K4,L0,X3,X4,X1+X2,KK);};
    auto A8=[&](const int K1,const int K2,const int K3,
                const int K4,const int K5,const Double B3){
        return tm3(N1,L1,K1,K2,M1,LA1,K3,K4,K5,X1,X2,B3,2);};
    auto A9=[&](const int K1,const int K2,const int K3,const int K4,
                const Double B1,const Double B3){
        return tm3(K1,K2,K3,K4,M2,LA2,M3,LA3,L0T,B1,X2,B3,1);};
    auto A10=[&](const int K1,const int K2,const int K3,
                 const Double B3,const int KK){
        return tm3(N2,L2,K1,K2,M1,LA1,M2,LA2,K3,X3,X4,B3,KK);};
    auto A11=[&](const int K1,const int K2,const int K3,const int KK){
        return tm3(N1,L1,N3,L3,K1,K2,M3,LA3,K3,X1,X2,X3+X4,KK);};

    const Double SS=std::sqrt((2*L0+1)*(2*L0T+1));
    Double S=0.;
    switch(K){
        case 11123:
        case 22123:
            return (N1==M1)*(N2==M2)*(N3==M3)
            *(L1==LA1)*(L2==LA2)*(L3==LA3)*(L0==L0T);
        case 11124:
            if(!(N1==M1 && L1==LA1 && L0==L0T)) return 0.;
            return A1(X1+X2,X3,X4,2)*par(LA2);
        case 12123:
            if(!(N1==M1 && L1==LA1 && L0==L0T)) return 0.;
            return A1(X1+X2,X3,X4,1)*par(LA3);
        case 21124:
            if(!(N1==M1 && L1==LA1 && L0==L0T)) return 0.;
            return A1(X3,X4,X1+X2,1)*par(L3+LA2);
        case 21123:
            if(!(N1==M1 && L1==LA1 && L0==L0T)) return 0.;
            return A1(X3,X4,X1+X2,2)*par(L3);
        case 21342:
            if(!(N2==M1 && L2==LA1)) return 0.;
            return SS*W1(L3,L2,L0T)
            *tm3(N1,L1,N3,L3,M2,LA2,M3,LA3,L0T,X1,X2,X3+X4,1)
            *par(L+L0+L0T+L3+LA2);
        case 21341:
            if(!(N2==M1 && L2==LA1)) return 0.;
            return SS*W1(L3,L2,L0T)
            *tm3(N1,L1,N3,L3,M2,LA2,M3,LA3,L0T,X1,X2,X3+X4,2)
            *par(L+L0+L0T+L3);
        case 11341:
        case 11342:
#pragma omp parallel for reduction(+:S)
            for(int M=0;M<=N2+N3+(L2+L3+1)/2;M++){
                //for(int N=abs(L1-L0T);N<=L1+L0T;N++){
                //if((2*(N2+N3)+L2+L3) != (2*(M1+M)+LA1+N) ||
                //(2*(N1+M)+L1+N) != (2*(M2+M3)+LA2+LA3)) continue;
                const int N=2*(N2+N3-M1-M)+L2+L3-LA1;
                if (N!=2*(M2+M3-N1-M)+LA2+LA3-L1) continue;
                if(N<std::abs(L1-L0T) || N>L1+L0T) continue;
                Double C=0.;
                if(K==11341) C=A2(M,N,2);
                else C=A2(M,N,1);
                S+=W2(N,L1,L0)*C
                *tm3(N2,L2,N3,L3,M1,LA1,M,N,L0,X1+X2,X3,X4,1);
            }
            return SS*S*par(L+L0+L0T+((11342==K)?LA2:0));
        case 11231:
        case 11132:
            if((2*(N1+N2)+L1+L2)!=(2*(M1+M2)+LA1+LA2)) return 0.;
            if(!(N3==M3 && L3==LA3)) return 0.;
#pragma omp parallel for reduction(+:S)
            for(int I=std::abs(L1-L2);I<=L1+L2;I++){
                int C=0.;
                if(K==11231) C=tm3(N1,L1,N2,L2,M1,LA1,M2,LA2,I,X1,X2,X3,1);
                else C=tm3(N1,L1,N2,L2,M1,LA1,M2,LA2,I,X1,X2,X3,2)*par(LA1);
                S+=(2*I+1)*W1(L2,L3,I)*W2(LA2,L3,I)*C;
            }
            return SS*S;
        case 11142:
        case 11241:
        case 21142:
        case 21132:
        case 21231:
        case 21241:
            for(int M=0;M<=N2+N3+(L2+L3+1)/2;M++){
                //for(int N=abs(LA3-L0);N<=LA3+L0;N++){
                //if((2*(N1+M)+L1+N) != (2*(M1+M2)+LA1+LA2) ||
                //(2*(N2+N3)+L2+L3) != (2*(M+M3)+N+LA3)) continue;
                int N=2*(M1+M2-N1-M)+LA1+LA2-L1;
                if(N!=2*(N2+N3-M-M3)+L2+L3-LA3) continue;
                if(N<std::abs(LA3-L0) || N>LA3+L0) continue;
#pragma omp parallel for reduction(+:S)
                for(int I=std::abs(LA1-LA2);I<=LA1+LA2;I++){
                    Double C=0.;
                    Double CC=0.;
                    switch(K){
                        case 11142:
                            CC=A3(M,N,I,X4,2);
                            C=A4(M,N,X1+X2,X3,X4)*par(LA1+N);
                            break;
                        case 11241:
                            CC=A3(M,N,I,X4,1);
                            C=A4(M,N,X1+X2,X3,X4)*par(N);
                            break;
                        case 21142:
                            CC=A3(M,N,I,X4,2);
                            C=A4(M,N,X4,X3,X1+X2)*par(L3+LA1+L2);
                            break;
                        case 21241:
                            CC=A3(M,N,I,X4,1);
                            C=A4(M,N,X4,X3,X1+X2)*par(L3+L2);
                            break;
                        case 21132:
                            CC=A3(M,N,I,X3,2);
                            C=A4(M,N,X3,X4,X1+X2)*par(L3+LA1);
                            break;
                        case 21231:
                            CC=A3(M,N,I,X3,1);
                            C=A4(M,N,X3,X4,X1+X2)*par(L3);
                            break;
                    }
                    S+=(2*I+1)*W1(N,LA3,I)*W2(LA2,LA3,I)*CC*C;
                }
            }
            return SS*S;
        case 11134:
        case 11234:
        case 12132:
            for(int M=0;M<=N1+N2+(L1+L2+1)/2;M++){
                //for(int N=abs(L3-L0T);N<=L3+L0T;N++){
                //if((2*(N1+N2)+L1+L2) != (2*(M1+M)+LA1+N) ||
                //(2*(M+N3)+N+L3) != (2*(M2+M3)+LA2+LA3)) continue;
                int N=2*(N1+N2-M1-M)+L1+L2-LA1;
                if(N!=2*(M2+M3-M-N3)+LA2+LA3-L3) continue;
                if(N<std::abs(L3-L0T) || N>L3+L0T) continue;
#pragma omp parallel for reduction(+:S)
                for(int I=std::abs(L1-L2);I<=L1+L2;I++){
                    Double C=0.,CC=0.;
                    switch(K){
                        case 11134:
                            C=A5(M,N,I,2);
                            CC=A6(M,N,X1+X3,X2,2)*par(LA1+LA2);
                            break;
                        case 11234:
                            C=A5(M,N,I,1);
                            CC=A6(M,N,X2+X3,X1,2)*par(LA2);
                            break;
                        case 12132:
                            C=A5(M,N,I,2);
                            CC=A6(M,N,X1+X3,X2,1)*par(LA1+LA3);
                            break;
                    }
                    S+=(2*I+1)*W1(L2,L3,I)*W2(N,L3,I)*C*CC;
                }
            }
            return SS*S;
        case 11143:
        case 11243:
#pragma omp parallel for collapse(3) reduction(+:S)
            for(int J1=0;J1<=N2+N3+(L2+L3+1)/2;J1++)
                for(int J2=0;J2<=2*(N2+N3)+L2+L3;J2++)
                    for(int J3=0;J3<=N2+N3+(L2+L3+1)/2;J3++){
                        //for(int J4=abs(L0-J2);J4<=L0+J2;J4++)
                        //if((2*(N2+N3)+L2+L3) != (2*(J1+J3)+J2+J4)) continue;
                        const int J4=2*(N2+N3-J1-J3)+L2+L3-J2;
                        if(J4<std::abs(L0-J2) || J4>L0+J2) continue;
                        for(int I=std::abs(L1-J2);I<=L1+J2;I++)
                            for(int M=0;M<=M2+M3+(LA2+LA3+1)/2;M++){
                                //for(int N=abs(L0T-J4);N<=L0T+J4;N++){
                                //if((2*(M+J3)+N+J4) != (2*(M2+M3)+LA2+LA3) ||
                                //(2*(N1+J1)+L1+J2) != (2*(M1+M)+LA1+N))
                                //continue;
                                const int N=2*(M2+M3-M-J3)+LA2+LA3-J4;
                                if(N!=2*(N1+J1-M1-M)+L1+J2-LA1) continue;
                                if(N<std::abs(L0T-J4) || N>L0T+J4) continue;
                                Double U1=0.;
                                Double U2=0.;
                                if(K==11143){
                                    U1=X1;
                                    U2=X2;
                                }else{//K==11243
                                    U1=X2;
                                    U2=X1;
                                }
                                S+=(2.*I+1.)*W1(J2,J4,I)*W2(N,J4,I)
                                *tm3(N2,L2,N3,L3,J1,J2,J3,J4,L0,X1+X2,X3,X4,2)
                                *tm3(M,N,J3,J4,M2,LA2,M3,LA3,L0T,U1+X4,U2,X3,2)
                                *tm3(N1,L1,J1,J2,M1,LA1,M,N,I,U1,U2,X4,2)
                                *par(J2);
                            }
                    }
            if(K==11243) S*=par(L1);
            return SS*S*par(LA1+LA2);
        case 22132:
        case 22142:
#pragma omp parallel for collapse(3) reduction(+:S)
            for(int J1=0;J1<=N2+N3+(L2+L3+1)/2;J1++)
                for(int J2=0;J2<=2*(N2+N3)+L2+L3;J2++)
                    for(int J3=0;J3<=N2+N3+(L2+L3+1)/2;J3++){
                        //for(int J4=abs(J2-L0);J4<=J2+L0;J4++)
                        //if(2*(N2+N3)+L2+L3!=2*(J1+J3)+J2+J4) continue;
                        const int J4=2*(N2+N3-J1-J3)+L2+L3-J2;
                        if(J4<std::abs(J2-L0) || J4>J2+L0) continue;
                        for(int I=std::abs(L1-J2);I<=L1+J2;I++)
                            for(int M=0;M<=N1+J1+(L1+J2+1)/2;M++){
                                //for(int N=abs(J4-L0T);N<=J4+L0T;N++){
                                //if((2*(N1+J1)+L1+J2) != (2*(M1+M)+LA1+N) ||
                                //(2*(M+J3)+ N+J4) != (2*(M2+M3)+LA2+LA3))
                                //continue;
                                const int N=2*(N1+J1-M1-M)+L1+J2-LA1;
                                if(N!=2*(M2+M3-M-J3)+LA2+LA3-J4) continue;
                                if(N<std::abs(J4-L0T) || N>J4+L0T) continue;
                                Double AA=0.,BB=0.,C=0.;
                                if(K==22132){
                                    AA=A7(J1,J2,J3,J4,2);
                                    BB=A8(J1,J2,M,N,I,X3);
                                    C=A9(M,N,J3,J4,X1+X3,X4);
                                }else{//K==22142
                                    AA=A7(J1,J2,J3,J4,1);
                                    BB=A8(J1,J2,M,N,I,X4);
                                    C=A9(M,N,J3,J4,X1+X4,X3)*par(J2);
                                }
                                S+=(2.*I+1.)*W1(J2,J4,I)*W2(N,J4,I)*AA*BB*C;
                            }
                    }
            return SS*S*par(L3+LA1+LA3);
        case 12142:
#pragma omp parallel for reduction(+:S)
            for(int M=0;M<=N1+N2+(L1+L2+1)/2;M++){
                //for(int N=0;N<=2*(N1+N2)+L1+L2;N++){
                //if(2*(N1+N2)+L1+L2!=2*(M2+M)   +LA2+N ||
                //2*(M+N3)+N+L3!=2*(M1+M3)+LA1+LA3) continue;
                const int N=2*(N1+N2-M2-M)+L1+L2-LA2;
                if(N!=2*(M1+M3-M-N3)+LA1+LA3-L3) continue;
                if(N<0 || N>2*(N1+N2)+L1+L2) continue;
                for(int I=std::abs(L1-L2);I<=L1+L2;I++)
                    for(int J=std::abs(LA1-LA3);J<=LA1+LA3;J++){
                        S+=(2*I+1)*(2*J+1)*W1(L2,L3,I)*W2(LA3,LA2,J)
                        *WIG(LA2,N,L,L3,I,J)
                        *tm3(N1,L1,N2,L2,M2,LA2,M,N,I,X1,X2,X3,1)
                        *tm3(M,N,N3,L3,M1,LA1,M3,LA3,J,X2+X3,X1,X4,1)
                        *par(J);
                    }
            }
            return SS*S*par(L+L0T+LA3);
        case 21243:
        case 21134:
        case 21143:
        case 21234:
#pragma omp parallel for reduction(+:S)
            for(int I=std::abs(L1-L3);I<=L1+L3;I++){
                Double Wsum=0;
                for(int J=std::abs(LA1-LA2);LA1+LA2;J++){
                    Double Tsum=0;
                    for(int M=0;M<=N1+N3+(L1+L3+1)/2;M++){
                        //for(int N=0;N<=2*(N1+N3)+L1+L3;N++){
                        //if(2*(N1+N3)+L1+L3!=2*(M+M3)+N+LA3 ||
                        //2*(N2+M)+L2+N!=2*(M1+M2)+LA1+LA2) continue;
                        const int N=2*(N1+N3-M-M3)+L1+L3-LA3;
                        if(N!=2*(M1+M2-N2-M)+LA1+LA2-L2) continue;
                        if(N<0 || N>2*(N1+N3)+L1+L3) continue;
                        const Double W3=WIG(L2,N,L,LA3,J,I);
                        Double BB=0.;
                        Double C=0.;
                        switch(K){
                            case 21243:
                                BB=A10(M,N,J,X2,1)*par(I+L+L0+L3+LA1+N);
                                C=A11(M,N,I,1);
                                break;
                            case 21134:
                                BB=A10(M,N,J,X1,2)*par(I+L+L0+L3);
                                C=A11(M,N,I,2);
                                break;
                            case 21143:
                                BB=A10(M,N,J,X1,1)*par(I+L+L0+L3+LA1);
                                C=A11(M,N,I,2);
                                break;
                            case 21234:
                                BB=A10(M,N,J,X2,2)*par(N+I+L+L0+L3);
                                C=A11(M,N,I,1);
                                break;
                        }
                        Tsum+=W3*BB*C;
                    }
                    Wsum+=(2*J+1)*W2(LA2,LA3,J)*Tsum;
                }
                S+=(2*I+1)*W1(L3,L2,I)*Wsum;
            }
            return SS*S;
        default:
#ifdef TESTING_STATES_JLG
            throw std::invalid_argument("K not on approved list");
#else
            return 0.;
#endif
    }
}

Double tm4a(const int n1,const int l1,
            const int n2,const int l2,
            const int n3,const int l3,
            const int n1p,const int l1p,
            const int n2p,const int l2p,
            const int n3p,const int l3p,
            const int l0,const int l0p,const int L,
            const Double m1,const Double m2,const Double m3,const Double m4,
            const int K){
    //4 body TMC for K configuration only
    //l0 couples l1 and l2 instead of l2 and l3
    //only transformations that rotate 1,2,3 or move 4 onto r1
    if(0>n1 or 0>l1 or 0>n2 or 0>l2 or 0>n3 or 0>l3 or
       0>n1p or 0>l1p or 0>n2p or 0>l2p or 0>n3p or 0>l3p or
       0>l0 or 0>l0p or 0>L or 0>m1 or 0>m2 or 0>m3 or 0>m4){
#ifdef TESTING_STATES_JLG
        throw std::invalid_argument("nonnegative arguments");
#else
        return 0;
#endif
    }
    Double SS=std::sqrt(Double((2*l0+1)*(2*l0+1)));
    switch (K) {
        case 1234:
            return (n1==n1p)*(l1==l1p)*(n2==n2p)*(l2==l2p)*(n3==n3p)*(l3==l3p)
            *(l0==l0p);
        case 2314:
            if((n3!=n3p) or (l3!=l3p) or (l0!=l0p)) return 0;
            return tm3(n1,l1,n2,l2,n1p,l1p,n2p,l2p,l0,m1,m2,m3,1);
        case 3124:
            if((n3!=n3p) or (l3!=l3p) or (l0!=l0p)) return 0;
            return tm3(n1,l1,n2,l2,n1p,l1p,n2p,l2p,l0,m1,m2,m3,2);
        case 2413:
        case 4123:
        {
            Double Wsum=0;
#pragma omp parallel for reduction(+:Wsum)
            for(int l0pp=std::abs(l2-l3);l0pp<=l2+l3;l0pp++){
                Double Tsum=0;
                const int Hsum=2*(n2+n3-n3p)+l2+l3-l3p;
                for(int l2pp=Hsum%2;l2pp<=Hsum;l2pp+=2){
                    const int n2pp=(Hsum-l2pp)/2;

                    Tsum+=W6J(l1,l2pp,l0p,l3p,L,l0pp)
                    *tm3(n2,l2,n3,l3,n2pp,l2pp,n3p,l3p,l0pp,m1+m2,m3,m4,2)
                    *tm3(n1,l1,n2pp,l2pp,n1p,l1p,n2p,l2p,l0p,m1,m2,m4,
                         (2413==K)?1:2);
                }
                Wsum+=(2*l0pp+1)*W6J(l1,l2,l0,l3,L,l0pp)*Tsum;
            }
            return par(l2+l3+l3+L)*SS*Wsum;
        }
        case 3421:
        {
            Double L2sum=0;
#pragma omp parallel for reduction(+:L2sum)
            for(int l0pp=std::abs(l2-l3);l0pp<=l2+l3;l0pp++){
                Double L3sum=0;
                for(int l03=std::abs(l3p-l2p);l03<=l3p+l2p;l03++){
                    Double Tsum=0;
                    const int Hsum=2*(n2+n3-n1p)+l2+l3-l1p;
                    for(int l3pp=Hsum%2;l3pp<=Hsum;l3pp+=2){
                        const int n3pp=(Hsum-l3pp)/2;

                        Tsum+=W6J(l1,l3pp,l03,l1p,L,l0pp)
                        *tm3(n2,l2,n3,l3,n1p,l1p,n3pp,l3pp,l0pp,m1+m2,m3,m4,1)
                        *tm3(n1,l1,n3pp,l3pp,n2p,l2p,n3p,l3p,l03,m1,m2,m3+m4,1);
                    }
                    L3sum+=(2*l03+1)*W6J(l1p,l2p,l0p,l3p,L,l03)*Tsum;
                }
                L2sum+=par(l0pp)*(2*l0pp+1)*W6J(l1,l2,l0,l3,L,l0pp)*L3sum;
            }
            return par(l1p+l2+l3+L)*SS*L2sum;
        }
        default:
#ifdef TESTING_STATES_JLG
            throw std::invalid_argument("K not on list");
#else
            return 0.;
#endif
    }
}

Double tm4b(const int n1,const int l1,
            const int n2,const int l2,
            const int n3,const int l3,
            const int n1p,const int l1p,
            const int n2p,const int l2p,
            const int n3p,const int l3p,
            const int l0,const int l0p,const int L,
            const Double m1,const Double m2,const Double m3,const Double m4,
            const int K){
    //4 body TMC
    //l0 couples l1 and l2 instead of l2 and l3
    //sums over intermediate l0's
    if(0>n1 or 0>l1 or 0>n2 or 0>l2 or 0>n3 or 0>l3 or
       0>n1p or 0>l1p or 0>n2p or 0>l2p or 0>n3p or 0>l3p or
       0>l0 or 0>l0p or 0>L or 0>m1 or 0>m2 or 0>m3 or 0>m4){
#ifdef TESTING_STATES_JLG
        throw std::invalid_argument("nonnegative arguments");
#else
        return 0;
#endif
    }

    Double S=0;
#pragma omp parallel for reduction(+:S)
    for(int k=std::abs(l2-l3);k<=l2+l3;k++){
        Double T=0;
        for(int l=std::abs(l2p-l3p);l<=l2p+l3;l++){
            T+=RW(l1p,l2p,L,l3p,l0p,l)*
            tm4(n1,l1,n2,l2,n3,l3,n1p,l1p,n2p,l2p,n3p,l3p,k,l,L,m1,m2,m3,m4,K);
        }
        S+=RW(l1,l2,L,l3,l0,k)*T;
    }
    return S*std::sqrt(Double((2*l0+1)*(2*l0p+1)));
}
#endif //TM4_JLG


Double spin3(const int S,
             const int P1,
             const int P2,
             const int k)
noexcept(
#ifdef TESTING_STATES_JLG
         false
#else
         true
#endif
){
    if(!(1==k or 2==k or 3==k)){
#ifdef TESTING_STATES_JLG
        throw std::invalid_argument("k must be 1, 2, or 3");
#else
        return 0;
#endif
    }
    if(!(1==S or 3==S)){
#ifdef TESTING_STATES_JLG
        throw std::invalid_argument("S must be 1/2 or 3/2");
#else
        return 0;
#endif
    }
    if(!((1==P1 or 0==P1) and (1==P2 or 0==P2))){
#ifdef TESTING_STATES_JLG
        throw std::invalid_argument("P must be 0 or 1");
#else
        return 0;
#endif
    }
    if(3==S){
        if(!(0==P1 and 0==P2)){
#ifdef TESTING_STATES_JLG
            throw std::invalid_argument("P must be 0 if S=3/2");
#else
            return 0;
#endif
        } else return 1;
    }
    //now 1==S
    if(3==k) return P1==P2;
    if(P1==P2) return -0.5;
    return par(1+k+P1)*Ms3_2;
}
