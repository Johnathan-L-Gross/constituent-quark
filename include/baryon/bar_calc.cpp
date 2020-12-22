//
//  bar_calc.cpp
//
//
//  Created by Johnathan Gross
//
//

#include "bar_calc.h"
//#include "functions.h"
#define W_JLG_map
#include "algebra/algebra.h"
#ifdef TESTING_STATES_JLG
#include <iostream>
#endif
//#include "mout.h"

/*
typedef std::vector<CDmatrix> CD3;
typedef std::vector<CD3> CD4;
typedef std::vector<CD4> CD5;
typedef std::vector<CD5> CD6;
 */

inline constexpr int SMIN_B=1;
inline constexpr int SMAX_B=3;
inline constexpr long double Ms3_4=0.433012701892219323381861585376468091735701313452595157013951L;
inline constexpr long double Ms5_s3=1.290994448735805628393088466594133203610973901763863608862524L;
inline constexpr long double Ms2_s3=0.816496580927726032732428024901963797321982493552223376144L;
inline constexpr long double M1_s3=0.577350269189625764509148780501957455647601751270126876018L;
inline constexpr long double M1_s2=0.707106781186547524400844362104849039284835937688474036588L;
inline constexpr long double M8_PI=2.54647908947032537230214021396022979255135433184730317996267750494234876214762456144182084426004937529716548362812L;
inline constexpr long double MPI_2=1.5707963267948966192313216916397514420985846996875529104874722961539082031431044993140174126710585339910740432566411533235469223L;
//inline constexpr long double  Ms3_s2=1.224744871391589049098642037352945695982973740328335064216L;
//already defined in polynomials.t

//******************************************************************************
//declaration of functions not in baryon.h

void term0(sparse&,//ham
           const imatrix&,//states
           const int1&,//params
           Double,//a
           const marray&,//mass
           const func1&,//f
           const int,//k
           const bool);//
//1D scalar
#ifdef D3_JLG
void term0(sparse&,//ham
           const imatrix&,//states
           const int3&,//params
           Double,//ax
           Double,//ay
           const marray&,//mass
           const func3&,//f
           const int,//k
           const bool);//sym
//3D scalar
#endif
#ifdef D6_JLG
void term0(CDmatrix&,//ham
           const imatrix&,//states
           const int6&,//params
           Double,//ax
           Double,//ay
           const marray&,//mass
           const func6&,//f
           const int,//k
           const bool);//sym
//6D scalar
#endif

void term1(sparse&,//ham
           const imatrix&,//states
           const int1&,//params
           Double,//a
           const marray&,//mass
           const func1&,//f
           const int,//k
           const bool);//sym
//1D spin-spin
#ifdef D3_JLG
void term1(sparse&,//ham
           const imatrix&,//states
           const int3&,//params
           Double,//ax
           Double,//ay
           const marray&,//mass
           const func3&,//f
           const int,//k
           const int,//spin
           const bool);//sym
//3D spin-spin
#endif
#ifdef D6_JLG
void term1(CDmatrix&,//ham
           const imatrix&,//states
           const int6&,//params
           Double,//ax
           Double,//ay
           const marray&,//mass
           const func6&,//f
           const int,//k
           const int,//spin
           const bool);//sym
//6D spin-spin
#endif
/*allowed values for spin:
 11,12,13,21,22,23,31,32,33
 only 33,12 allowed for 1D
 */

void term2(sparse&,//ham
           const imatrix&,//states
           const int1&,//params
           const Double,//a
           const marray&,//mass
           const func1&,//f
           const bool);//sym
//1D tensor

void term3(sparse&,//ham
           const imatrix&,//states
           const int1&,//params
           const Double,//a
           const marray&,//mass
           const func1&,//f
           const marray&);//pars
//1D vector 2 body

void term4(sparse&,//ham
           const imatrix&,//states
           const int1&,//params
           const Double,//ax
           const Double,//ay
           const marray&,//mass
           const func1&,//f
           const marray&);//pars
//1D vector 3 body

void term5(sparse&,//ham
           const imatrix&,//states
           //const int1&,//params
           Double,//a
           const marray&,//mass
           const func1&,//f
           const int,//k
           const int,//spin
           const bool);//sym
//1D Dirac delta potential spin-spin


inline void term(int I,
                 sparse&ham,
                 const imatrix&states,
                 const int1&params,
                 const Double alpha_rho,
                 const Double alpha_lambda,
                 const marray&mass,
                 const func1&f,
                 const int k,
                 const int spin,
                 const marray&pars,
                 const bool sym){
    //1D selector

    //bool var=(0b10&k);
    //std::cerr<<"term 1"<<std::endl;
    const Double alpha=(k&0b10)?alpha_rho:alpha_lambda;
    switch (I) {
        case 0://scalar
            term0(ham,states,params,alpha,mass,f,k,sym);
            return;
        case 1://S.S
#ifdef TESTING_STATES_JLG
            if(!(0b10&k)) throw std::invalid_argument("must be rho\n");
            if(12!=spin and 21!=spin) throw std::invalid_argument("S.S must be 12");
#endif
            term1(ham,states,params,alpha,mass,f,k,sym);
            return;
        case 2://spin tensor
#ifdef TESTING_STATES_JLG
            if(!(0b10&k)) throw std::invalid_argument("var must be rho\n");
            if(12!=spin and 21!=spin) throw std::invalid_argument("S2 must be 12");
#endif
            term2(ham,states,params,alpha,mass,f,sym);
            return;
        case 3://vector 2 body
#ifdef TESTING_STATES_JLG
            if(!(0b10&k)) throw std::invalid_argument("var must be rho\n");
            if(sym) throw std::invalid_argument("2B vector not symmetric\n");
#endif
            term3(ham,states,params,alpha,mass,f,pars);
            return;
        case 4://vector 3 body
#ifdef TESTING_STATES_JLG
            if(!(0b10&k)) throw std::invalid_argument("var must be rho\n");
            if(sym) throw std::invalid_argument("3B vector not symmetric\n");
#endif
            term4(ham,states,params,alpha_rho,alpha_lambda,mass,f,pars);
            return;
        case 5://delta S.S
#ifdef TESTING_STATES_JLG
            if(!(0b10&k)) throw std::invalid_argument("must be rho\n");
#endif
            term5(ham,states,alpha,mass,f,k,spin,sym);
            return;
        default:
            return;
    }
}
#ifdef D3_JLG
inline void term(int I,
                 sparse&ham,
                 const imatrix&states,
                 const int3&params,
                 const Double alpha_rho,
                 const Double alpha_lambda,
                 const marray&mass,
                 const func3&f,
                 const int k,
                 const int spin,
                 const bool sym){
    //3D selector

    //std::cerr<<"term 3"<<std::endl;
    //const Double alphas=(k&0b10)?alphas_rho:alphas_lambda;
    switch (I) {
        case 0://scalar
            term0(ham,states,params,alpha_rho,alpha_lambda,mass,f,k,sym);
            return;
        case 1://S.S
            term1(ham,states,params,alpha_rho,alpha_lambda,mass,f,k,spin,sym);
            return;
#ifdef TESTING_STATES_JLG
        case 2://s tensor
            throw std::invalid_argument("3D can't do S^2.R^2\n");
        case 3://vector 2 body
            throw std::invalid_argument("3D can't do vector\n");
        case 4://vector 3 body
            throw std::invalid_argument("3D can't do vector\n");
        case 5://dirac S.S
            throw std::invalid_argument("3D can't do Dirac\n");
#endif
        default:
#ifdef TESTING_STATES_JLG
            throw std::invalid_argument("I not allowed\n");
#endif
            return;
    }
}
#endif
#ifdef D6_JLG
inline void term(int I,//6D selector
                 CDmatrix&ham,
                 const imatrix&states,
                 const int6&params,
                 const Double alpha_rho,
                 const Double alpha_lambda,
                 const marray&mass,
                 const func6&f,
                 const int k,
                 const int spin,
                 const bool sym){
#ifdef TESTING_STATES_JLG
    if(sym) std::cerr<<"6D may not be symmetric\n";
#endif
    //const Double alphas=(k&0b10)?alphas_rho:alphas_lambda;
    switch (I) {
        case 0://scalar
            term0(ham,states,params,alpha_rho,alpha_lambda,mass,f,k,sym);
            return;
        case 1://S.S
            term1(ham,states,params,alpha_rho,alpha_lambda,mass,f,k,spin,sym);
            return;
#ifdef TESTING_STATES_JLG
        case 2://spin tensor
            throw std::invalid_argument("6D can't do S^2.R^2\n");
        case 3://vector 2 body
            throw std::invalid_argument("6D can't do vector\n");
        case 4://vector 3 body
            throw std::invalid_argument("6D can't do vector\n");
        case 5://dirac S.S
            throw std::invalid_argument("3D can't do Dirac\n");
#endif
        default:
            return;
    }
}
#endif
//inline int abs(int x){return (0>x?-x:x);}


//******************************************************************************
//Hamiltonian calculations

void hamcalc(std::vector<sparse>&ham,//1D
             const imatrix&states,
             const int1&params,
             const Double omega,
             const marray&mass,
             const std::vector<func1>&f,
             //const bvector&var,
             const ivector&k,
             const ivector&tens,
             const ivector&spin,
             const ivector&rot,
             const std::vector<marray>&pars,
             const bool sym,
             const Double error){
    //std::cerr<<"hamcalc 1\n";
    std::size_t s=f.size();
#ifdef TESTING_STATES_JLG
    if(0==s) throw std::invalid_argument("no functions to evaluate");
    if(k.size()!=s or tens.size()!=s or spin.size()!=s or rot.size()!=s
       or pars.size()!=s)
        throw std::invalid_argument("not same size\n");
#else
    if(0==s){
        ham.resize(0);
        return;
    }
#endif
    //std::cerr<<"hamcalc 1"<<std::endl;
    //std::cerr<<"s="<<s<<'\n';


    ham.resize(s,sparse(states.size(),states.size()));
    for(std::size_t i=0;i<s;i++){
        term(tens[i],ham[i],states,params,
             sqrt(omega*mrho(mrot(rot[i],mass))),
             sqrt(omega*mlambda(mrot(rot[i],mass))),
             mrot(rot[i],mass),f[i],k[i],spin[i],pars[i],sym);
        if(error) ham[i].set_error(error,true);
    }
    //std::cerr<<ham[1][5][2]<<std::endl;
}

#ifdef D3_JLG
void hamcalc(std::vector<sparse>&ham,//3D
             const imatrix&states,
             const int3&params,
             const Double omega,
             const marray&mass,
             const std::vector<func3>&f,
             const ivector&k,
             const ivector&tens,
             const ivector&spin,
             const ivector&rot,
             const bool sym,
             const Double error){
    //std::cerr<<"hamcalc 3\n";
    std::size_t s=f.size();
#ifdef TESTING_STATES_JLG
    if(0==s) throw std::invalid_argument("no functions to evaluate");
    if(k.size()!=s or tens.size()!=s or spin.size()!=s or rot.size()!=s)
        throw std::invalid_argument("not same size\n");
#else
    if(0==s){
        ham.resize(0);
        return;
    }
#endif
    //std::cerr<<"s="<<s<<'\n';

    ham.resize(s,sparse(states.size(),states.size(),0));
    for(std::size_t i=0;i<s;i++){
        term(tens[i],ham[i],states,params,
             sqrt(omega*mrho(mrot(rot[i],mass))),
             sqrt(omega*mlambda(mrot(rot[i],mass))),
             mrot(rot[i],mass),f[i],k[i],spin[i],sym);
        if(error) ham[i].set_error(error,true);
    }
}
#endif//D3_JLG
#ifdef D6_JLG
void hamcalc(std::vector<CDmatrix>&ham,//6D
             const imatrix&states,
             const int6&params,
             const Double omega,
             const marray&mass,
             const std::vector<func6>&f,
             const ivector&k,
             const ivector&tens,
             const ivector&spin,
             const ivector&rot,
             const bool sym,
             const Double error){
    std::size_t s=f.size();
#ifdef TESTING_STATES_JLG
    if(0==s) throw std::invalid_argument("no functions to evaluate");
    if(k.size()!=s or tens.size()!=s or spin.size()!=s or rot.size()!=s)
        throw std::invalid_argument("not same size\n");
#else
    if(0==s){
        ham.resize(0);
        return;
    }
#endif
    //std::cerr<<"s="<<s<<'\n';

    ham.resize(s,CDmatrix(states.size(),CDvector(states.size(),0)));
    for(std::size_t i=0;i<s;i++){
        term(tens[i],ham[i],states,params,
             sqrt(omega*mrho(mrot(rot[i],mass))),
             sqrt(omega*mlambda(mrot(rot[i],mass))),
             mrot(rot[i],mass),f[i],k[i],spin[i],sym);
        if(error){
#pragma omp parallel for collapse(2)
            for(std::size_t m=0;m<ham[i].size();m++)
                for(std::size_t n=0;n<ham[i].size();n++){
                    bool br=true,bi=true;
                    if(std::abs(ham[i][m][n].real())<error) br=false;
                    if(std::abs(ham[i][m][n].imag())<error) bi=false;
                    if(br and bi);
                    else if(br xor bi) ham[i][m][n]=
                        CDouble((br?ham[i][m][n].real():0),
                                (bi?ham[i][m][n].imag():0));
                    else ham[i][m][n]=0;
                }
        }
    }
}
#endif//D6_JLG
//******************************************************************************
//term0

inline bool sel0(const imatrix&states,
                 const std::size_t m,
                 const std::size_t n){
    if(states[m][0]!=states[n][0]) return false;//Pf
    if(states[m][1]!=states[n][1]) return false;//J
    if(states[m][2]!=states[n][2]) return false;//Pj
    if(states[m][3]!=states[n][3]) return false;//M
    if(states[m][5]!=states[n][5]) return false;//L
    if(states[m][10]!=states[n][10]) return false;//S
    if(states[m][11]!=states[n][11]) return false;//Ps
    return true;
}

void term0(sparse&ham,//1D
           const imatrix&states,
           const int1&params,
           Double a,
           const marray&mass,
           const func1&f,
           const int k,
           const bool sym){
    const bool var=(k & 0b10);
    const bool mom=(k & 0b1);
    //std::cerr<<"term0 1 "<<var<<' '<<k<<' '<<mom<<std::endl;

    if(mom) a=1/a;

    Dvector fx(params.n());
#pragma omp parallel for
    for(size_t x=0;x<params.n();x++){
        //std::cerr<<x<<',';
        const Double xx=(params.polytype()?
                         params.x(x):std::sqrt(params.x(x)));
        fx[x]=f(xx/a,mass);
        //std::cerr<<'\n';
    }
    //std::cerr<<"fx"<<std::endl;

    const std::size_t s=states.size();
    ham=sparse(s,s);

#pragma omp parallel for schedule(dynamic)
    for(size_t m=0;m<s;m++)for(size_t n=(sym?m:0);n<s;n++){
        //std::cerr<<m<<' '<<n<<'\n';
        if(!sel0(states,m,n)) continue;
        if(states[m][6]!=states[n][6] and !var) continue;//np
        if(states[m][7]!=states[n][7]) continue;//lp
        if(states[m][8]!=states[n][8] and var) continue;//nl
        if(states[m][9]!=states[n][9]) continue;//ll

        Double hmn=0.;
        for(size_t x=0;x<params.n();x++){
            hmn+=params.dx(x)*fx[x]
            *params.rad(states[m][var?6:8],states[m][var?7:9],x)
            *params.rad(states[n][var?6:8],states[n][var?7:9],x);
            //std::cerr<<x<<' '<<hmn<<'\n';
        }
        if(mom)
            hmn*=parity(states[m][6]+states[m][8]-states[n][6]-states[n][8]);
#pragma omp critical(term01D)
        {
            ham.set_entry(m,n,hmn);
            if(sym and m!=n) ham.set_entry(n,m,hmn);
        }
/*#pragma omp critical
        if(ham[m][n]!=0 and m!=n) std::cerr<<var<<' '<<mom<<' '<<m<<' '<<n<<' '<<ham[m][n]<<std::endl;*/
    }//m,n
    //std::cerr<<"exit term 1\n";
}
#ifdef D3_JLG
void term0(sparse&ham,//3D
           const imatrix&states,
           const int3&params,
           Double ax,
           Double ay,
           const marray&mass,
           const func3&f,
           const int k,
           const bool sym){
    //std::cerr<<"term0 3 "<<f.size()<<'\n';
    
    const bool mom=(k & 0b1);
    if(mom){
        ax=1/ax;
        ay=1/ay;
    }

    std::vector<Dmatrix> fxyt(params.n(),
                              Dmatrix(params.n(),
                                      Dvector(params.l(),0.)));
    std::vector<Dmatrix> fxyu(params.n(),
                              Dmatrix(params.n(),
                                      Dvector(params.l(),0.)));
#pragma omp parallel for
    for(size_t x=0;x<params.n();x++){
        const Double xx=(params.polytype()?params.x(x):std::sqrt(params.x(x)));
        for(size_t y=0;y<params.n();y++){
            const Double yy=(params.polytype()?params.x(y):std::sqrt(params.x(y)));
            for(size_t t=0;t<params.l();t++){
                const Double tt=params.t(t);
                fxyt[x][y][t]=f(xx/ax,yy/ay,tt,mass);
                fxyu[x][y][t]=f(xx/ax,yy/ay,-tt,mass);
            }//t
        }//y
    }//x

    std::size_t s=states.size();
    ham=sparse(s,s);

#pragma omp parallel for schedule(dynamic)
    for(size_t m=0;m<s;m++)for(size_t n=(sym?m:0);n<s;n++){
        //std::cerr<<'\n'<<m<<','<<n;
        if(!sel0(states,m,n)) continue;
//std::cerr<<" calc3";

        //const int L1=states[m][5]*2;
        //const int L2=states[n][5]*2,
        const int L=states[m][5]*2,
        lp1=states[m][7]*2,
        lp2=states[n][7]*2,
        ll1=states[m][9]*2,
        ll2=states[n][9]*2;
        const int Lmax=std::min(lp1+lp2,ll1+ll2),
        Lmin=std::max(std::abs(lp1-lp2),std::abs(ll1-ll2));

        const Double mnfactor=
        std::sqrt((lp1+1)*(lp2+1)*(ll1+1)*(ll2+1)/Double(2))
        *parity((L+lp1+lp2)/2)
        *(mom?
          parity((states[m][6]+states[m][8]-states[n][6]-states[n][8]
                  +(states[m][7]+states[m][9]-states[n][7]-states[n][9])/2))
          :1);
        
        Dvector Kfactor((Lmax-Lmin)/2+1,0);
        for(int K=0;K<=(Lmax-Lmin)/2;K++){
            const int Lr=Lmin+K*2;
            Kfactor[K]=W3J(lp1,lp2,Lr,0,0,0)
            *W3J(ll1,ll2,Lr,0,0,0)
            *W6J(lp1,lp2,Lr,ll2,ll1,L)
            *sqrt(Lr+1);
        }
        
        Double tint=0;
        for(std::size_t t=0;t<params.l();t++){
            
            Double Ksum1=0,
            Ksum2=0;
            for(int K=0;K<=(Lmax-Lmin)/2;K++){
                Ksum1+=Kfactor[K]*params.Y(Lmin/2+K,t);
                Ksum2+=Kfactor[K]*params.Y2(Lmin/2+K,t);
            }//K
            if(0==Ksum1 and 0==Ksum2) continue;
            
            Double xint1=0,
            xint2=0;
            for(std::size_t x=0;x<params.n();x++){
                const Double dx=params.dx(x)
                *params.rad(states[m][6],states[m][7],x)
                *params.rad(states[n][6],states[n][7],x);
               
                Double yint1=0,
                yint2=0;
                for(std::size_t y=0;y<params.n();y++){
                    const Double dy=params.dx(y)
                    *params.rad(states[m][8],states[m][9],y)
                    *params.rad(states[n][8],states[n][9],y);
                    
                    if(0!=Ksum1) yint1+=fxyt[x][y][t]*dy;
                    if(0!=Ksum2) yint2+=fxyu[x][y][t]*dy;
                }//y
                if(0!=Ksum1) xint1+=yint1*dx;
                if(0!=Ksum2) xint2+=yint2*dx;
            }//x
            tint+=(Ksum1*xint1+Ksum2*xint2)*params.dt(t);
        }//t
#pragma omp critical(term03)
        {
            ham.set_entry(m,n,tint*mnfactor);
            if(sym and m!=n) ham.set_entry(n,m,tint*mnfactor);
        }
    }//m,n
    //std::cerr<<'\n';
}
#endif//D3_JLG
#ifdef D6_JLG
void term0(CDmatrix&ham,//6D
           const imatrix&states,
           const int6&params,
           Double ax,
           Double ay,
           const marray&mass,
           const func6&f,
           const int k,
           const bool sym){
    //std::cerr<<"term0 6"<<std::endl;
    
    const bool mom=(k & 0b1);
    if(mom){
        ax=1./ax;
        ay=1./ay;
    }

    const std::size_t st=states.size();
    ham=CDmatrix(st,CDvector(st,0));

    //if storing function values
    /*CD6 ff(0);
    bool b6=sqr(params.n()*params.l()*params.m())*4<1e8;
    if(b6){
        //Could take up too much space
        //200x200x(50x2)x(50x2)x10x10=4e10 entries
        //50x50x(10x2)x(10x2)x5x5=2.5e7 entries
        //std::cerr<<"creating ff"<<std::endl;
        ff=CD6(params.n(),//x
               CD5(params.n(),//y
                   CD4(2*params.l(),//s
                       CD3(2*params.l(),//t
                           CDmatrix(params.m(),//p
                                    CDvector(params.m(),0))))));//q
        //std::cerr<<"created ff"<<std::endl;
#pragma omp parallel for
        for(size_t x=0;x<params.n();x++){
            Double xx=params.x(x);//GLQx.at(GLQn)
            if(params.polytype()) xx*=xx;
            
            for(size_t y=0;y<params.n();y++){
                //cerr<<m<<','<<n<<' '<<ms<<' '<<x<<','<<y<<'\n';
                Double yy=params.x(y);//GLQx.at(GLQn)
                if(params.polytype()) yy*=yy;
                
                for(size_t s=0;s<2*params.l();s++){//GPQl
                    Double ss=params.t(s/2);//GPQx.at(GPQl)
                    if(s%2)ss=-ss;
                    //symmetric about ss=0
                    
                    for(size_t t=0;t<2*params.l();t++){
                        //std::cerr<<m<<','<<n<<' '<<ms<<' '<<x<<','<<y<<' '<<m1<<','<<m2<<' '<<s<<','<<t<<'\n';
                        Double tt=params.t(t/2);//GPQx.at(GPQl)
                        if(t%2)tt=-tt;
                        //symmetric about tt=0
                        
                        for(size_t p=0;p<params.m();p++){//GPQm
                            Double pp=params.phi(p);//GPQx.at(GPQm)
                            //if(p%2)pp=-pp;
                            //symmetric about pp=0
                            //pp=phix(pp);
                            //symmetric about pp=PI
                            
                            for(size_t q=0;q<params.m();q++){
                                //std::cerr<<'\n'<<m<<','<<n<<' '<<ms<<' '<<x<<','<<y<<' '<<m1<<','<<m2<<' '<<s<<','<<t<<' '<<p<<','<<q<<' ';
                                Double qq=params.phi(q);//GPQx.at(GPQm)
                                //if(q%2)qq=-qq;
                                //qq=phix(qq);
                                ff[x][y][s][t][p][q]=
                                f(xx/ax,yy/ay,ss,tt,pp,qq,mass);
                            }//q
                        }//p
                    }//t
                }//s
            }//y
        }//x
        std::cerr<<"written ff"<<std::endl;
    }
     */

    //std::cerr<<"6int\n";
#pragma omp parallel for schedule(dynamic)
    for(size_t m=0;m<st;m++)for(size_t n=(sym?m:0);n<st;n++){
        if(!sel0(states,m,n)) continue;
//#pragma omp critical(nmstart)
        //std::cerr<<m<<','<<n<<'\n';

        const int L=states[m][5]*2,
        lp1=states[m][7]*2,
        lp2=states[n][7]*2,
        ll1=states[m][9]*2,
        ll2=states[n][9]*2;
        const int l1max=std::min(lp1,ll1),
        l2max=std::min(lp2,ll2);
        
        Dvector CG1(l1max+1,0),
        CG2(l2max+1,0);
        for(int m1=0;m1<=l1max;m1++){
            const int mm=2*m1-l1max;
            CG1[m1]=CG(L,0,lp1,mm,ll1,-mm);
        }
        //dispvec(CG1,false,std::cerr);
        for(int m2=0;m2<=l2max;m2++){
            const int mm=2*m2-l2max;
            CG2[m2]=CG(L,0,lp2,mm,ll2,-mm);
        }
        //dispvec(CG2,false,std::cerr);
        
        CDouble sint=0;
        for(std::size_t s=0;s<params.l()*2;s++){
            const Double ss=(s%2?-params.t(s/2):params.t(s/2));
            
            CDouble tint=0;
            for(std::size_t t=0;t<params.l()*2;t++){
                const Double tt=(t%2?-params.t(t/2):params.t(t/2));
                
                CDouble pint=0;
                for(std::size_t p=0;p<params.m();p++){
                    const Double pp=params.phi(p);
                    
                    CDouble qint=0;
                    for(std::size_t q=0;q<params.m();q++){
                        const Double qq=params.phi(q);
                        
                        CDouble m1sum=0;
                        for(int m1=0;m1<=l1max;m1++){
                            if(!CG1[m1]) continue;
                            const int mp1=2*m1-l1max,
                            ml1=l1max-2*m1;
                            const int par1=(0>mp1?parity(mp1):1)
                            *(0>ml1?parity(ml1):1);
                            
                            m1sum+=CG1[m1]*par1
                            *(s%2?
                              params.Y2(lp1/2,std::abs(mp1)/2,s/2)
                              :params.Y(lp1/2,std::abs(mp1)/2,s/2))
                            *(t%2?
                              params.Y2(ll1/2,std::abs(ml1)/2,t/2)
                              :params.Y(ll1/2,std::abs(ml1)/2,t/2))
                            *params.phase(-mp1,p)
                            *params.phase(-ml1,q);
                        }//m1
                        if(CDouble(0)==m1sum) continue;
                        
                        CDouble m2sum=0;
                        for(int m2=0;m2<=l2max;m2++){
                            if(!CG2[m2]) continue;
                            const int mp2=2*m2-l2max,
                            ml2=l2max-2*m2;
                            const int par2=(0>mp2?parity(mp2):1)
                            *(0>ml2?parity(ml2):1);
                            
                            m2sum+=CG2[m2]*par2
                            *(s%2?
                              params.Y2(lp2/2,std::abs(mp2)/2,s/2)
                              :params.Y(lp2/2,std::abs(mp2)/2,s/2))
                            *(t%2?
                              params.Y2(ll2/2,std::abs(ml2)/2,t/2)
                              :params.Y(ll2/2,std::abs(ml2)/2,t/2))
                            *params.phase(mp2,p)
                            *params.phase(ml2,q);
                        }//m2
                        if(CDouble(0)==m2sum) continue;
                        
                        CDouble xint=0;
                        for(std::size_t x=0;x<params.n();x++){
                            const Double xx=(params.polytype()?
                                             params.x(x):
                                             std::sqrt(params.x(x)));
                            
                            CDouble yint=0;
                            for(std::size_t y=0;y<params.n();y++){
                                const Double yy=(params.polytype()?
                                                 params.x(y):
                                                 std::sqrt(params.x(y)));
                                
                                yint+=f(xx/ax,yy/ay,ss,tt,pp,qq,mass)
                                *params.dx(y)
                                *params.rad(states[m][8],states[m][9],y)
                                *params.rad(states[n][8],states[n][9],y);
                            }//y
                            xint+=yint*params.dx(x)
                            *params.rad(states[m][6],states[m][7],x)
                            *params.rad(states[n][6],states[n][7],x);
                        }//x
                        qint+=m1sum*m2sum*xint*params.dphi(q);
                    }//q
                    if(CDouble(0)!=qint) pint+=qint*params.dphi(p);
                }//p
                if(CDouble(0)!=pint) tint+=pint*params.dt(t/2);
            }//t
            if(CDouble(0)!=tint) sint+=tint*params.dt(s/2);
        }//s
        
        /*CDouble m1sum=0;
        for(int mp1=-l1max;mp1<=l1max;mp1+=2){
            //if(!CG1[m1]) continue;
            const int //mp1=2*m1-l1max,
            ml1=-mp1;//l1max-2*m1;
            const int par1=(0>mp1?parity(mp1/2):1)
            *(0>ml1?parity(ml1/2):1);
            
            CDouble m2sum=0;
            for(int mp2=-l2max;mp2<=l2max;mp2+=2){
                //if(!CG2[m2]) continue;
                const int //mp2=2*m2-l2max,
                ml2=-mp2;//l2max-2*m2;
                const int par2=(0>mp2?parity(mp2/2):1)
                *(0>ml2?parity(ml2/2):1);
                
                CDouble xint=0;
                for(std::size_t x=0;x<params.n();x++){
                    const Double xx=(params.polytype()?
                                     params.x(x):
                                     std::sqrt(params.x(x)));
                    
                    CDouble yint=0;
                    for(std::size_t y=0;y<params.n();y++){
                        const Double yy=(params.polytype()?
                                         params.x(y):
                                         std::sqrt(params.x(y)));
                        
                        CDouble sint=0;
                        for(std::size_t s=0;s<params.l()*2;s++){
                            const Double ss=(s%2?-params.t(s/2):params.t(s/2));
                            
                            CDouble tint=0;
                            for(std::size_t t=0;t<params.l()*2;t++){
                                const Double tt=(t%2?-params.t(t/2):params.t(t/2));
                                
                                CDouble pint=0;
                                for(std::size_t p=0;p<params.m();p++){
                                    const Double pp=params.phi(p);
                                    
                                    CDouble qint=0;
                                    for(std::size_t q=0;q<params.m();q++){
                                        const Double qq=params.phi(q);
                                        
                                        qint+=f(xx/ax,yy/ay,ss,tt,pp,qq,mass)
                                        *params.dphi(q)
                                        *params.phase(-ml1,q)
                                        *params.phase(ml2,q);
                                    }//q
                                    pint+=qint*params.dphi(p)
                                    *params.phase(-mp1,p)
                                    *params.phase(mp2,p);
                                }//p
                                tint+=pint*params.dt(t/2)
                                *(t%2?
                                  params.Y2(ll1/2,std::abs(ml1)/2,t/2)
                                  :params.Y(ll1/2,std::abs(ml1)/2,t/2))
                                *(t%2?
                                  params.Y2(ll2/2,std::abs(ml2)/2,t/2)
                                  :params.Y(ll2/2,std::abs(ml2)/2,t/2));
                            }//t
                            sint+=tint*params.dt(s/2)
                            *(s%2?
                              params.Y2(lp1/2,std::abs(mp1)/2,s/2)
                              :params.Y(lp1/2,std::abs(mp1)/2,s/2))
                            *(s%2?
                              params.Y2(lp2/2,std::abs(mp2)/2,s/2)
                              :params.Y(lp2/2,std::abs(mp2)/2,s/2));
                        }//s
                        yint+=sint*params.dx(y)
                        *params.rad(states[m][8],states[m][9],y)
                        *params.rad(states[n][8],states[n][9],y);
                    }//y
                    xint+=yint*params.dx(x)
                    *params.rad(states[m][6],states[m][7],x)
                    *params.rad(states[n][6],states[n][7],x);
                }//x
                m2sum+=xint*CDouble(CG(L,0,lp2,mp2,ll2,ml2))*CDouble(par2);
            }//m2
            m1sum+=m2sum*CDouble(CG(L,0,lp1,mp1,ll1,ml1))*CDouble(par1);
        }//m1*/
        
        ham[m][n]=sint
        *(mom?
          std::polar(Double(1),Double(MPI_2)*((states[m][6]+states[m][8]
                                               +states[n][6]+states[n][8])*2
                                              +states[m][7]+states[m][9]
                                              -states[n][7]-states[n][9])):
          CDouble(1));
        if(sym and m!=n) ham[n][m]=ham[m][n];
//#pragma omp critical(nmend)
        //std::cerr<<m<<','<<n<<" end\n";
    }//m,n
    //std::cerr<<'\n';
    //std::cerr<<"term0 6 end\n";
}
#endif//D6_JLG

//******************************************************************************
//term1 and term5

inline bool sel1(const imatrix&states,
                 const std::size_t m,
                 const std::size_t n){
    if(states[m][0]!=states[n][0]) return false;//Pf
    if(states[m][1]!=states[n][1]) return false;//J
    if(states[m][2]!=states[n][2]) return false;//Pj
    if(states[m][3]!=states[n][3]) return false;//M
    if(states[m][5]!=states[n][5]) return false;//L
    if(states[m][10]!=states[n][10]) return false;//S
    return true;
}

inline bool sel5(const imatrix&states,
                 const std::size_t m,
                 const std::size_t n){
    if(states[m][0]!=states[n][0]) return false;//Pf
    if(states[m][1]!=states[n][1]) return false;//J
    if(states[m][2]!=states[n][2]) return false;//Pj
    if(states[m][3]!=states[n][3]) return false;//M
    if(states[m][5]!=states[n][5]) return false;//L
    if(0!=states[m][7] or 0!=states[n][7]) return false;//lp
    if(states[m][8]!=states[n][8]) return false;//nl
    if(states[m][9]!=states[n][9]) return false;//ll
    if(states[m][10]!=states[n][10]) return false;//S
    return true;
}

#ifndef TESTING_STATES_JLG
constexpr
#endif
Double SdotS(const int spin,const int s1,const int p1,const int s2,const int p2){
    const bool P1=mod(p1,2);
    const bool P2=mod(p2,2);
    if((1!=s1 and 3!=s1) or (1!=s2 and 3!=s2)){
#ifdef TESTING_STATES_JLG
        std::cerr<<"not right spin: "<<s1<<' '<<s2<<std::endl;
        throw std::domain_error("S.S spins");
#endif
        return 0.;
    }
    if(s1!=s2) return 0.;
    switch (spin) {
        case 11:
        case 22:
        case 33:
            if(P1 != P2) return 0.;
            else return .75;
        case 12:
        case 21:
            if(P1 != P2) return 0.;
            else return (P1?-.75:.25);
        case 23:
        case 32:
            if(3==s1) return .25;
            else if(P1 != P2) return Ms3_4;
            else if(1==P1) return 0.;
            else return -.5;
        case 13:
        case 31:
            if(3==s1) return .25;
            else if(P1 != P2) return -Ms3_4;
            else if(1==P1) return 0.;
            else return -.5;

        default:
#ifdef TESTING_STATES_JLG
            std::cerr<<"not appropriate for S dot S: "
            <<spin<<' '<<s1<<','<<p1<<' '<<s2<<','<<p2<<'\n';
#endif
            return 0.;
    }
}

void term1(sparse&ham,//1D
           const imatrix&states,
           const int1&params,
           Double a,
           const marray&mass,
           const func1&f,
           const int k,
           const bool sym){
    //std::cerr<<"term1 1"<<std::endl;
    
    const bool mom=(k & 0b1);
    if(mom) a=1/a;

    Dvector fx(params.n(),0);
#pragma omp parallel for
    for(size_t x=0;x<params.n();x++){
        //std::cerr<<x<<':';
        const Double xx=(params.polytype()?
                         params.x(x):std::sqrt(params.x(x)));
        fx[x]+=f(xx/a,mass);
        //std::cerr<<xx<<' '<<fx[x]<<'\n';
    }
    //std::cerr<<"fx\n";

    const std::size_t s=states.size();
    ham=sparse(s,s);

#pragma omp parallel for schedule(dynamic)
    for(size_t m=0;m<s;m++)for(size_t n=(sym?m:0);n<s;n++){
        //std::cerr<<m<<' '<<n<<'\n';
        if(!sel1(states,m,n)) continue;
        if(states[m][7]!=states[n][7]) continue;//lp
        if(states[m][8]!=states[n][8]) continue;//nl
        if(states[m][9]!=states[n][9]) continue;//ll

        const Double SdS=SdotS(12,states[m][10],states[m][11],
                               states[n][10],states[n][11]);
        if(0==SdS) continue;

        Double hammn=0;
        for(size_t x=0;x<params.n();x++){
            hammn+=params.dx(x)*fx[x]
            *params.rad(states[m][6],states[m][7],x)
            *params.rad(states[n][6],states[n][7],x);
        }
        hammn*=SdS;
        if(mom) hammn*=parity(states[m][6]-states[n][6]);
#pragma omp critical(term11)
        {
            ham.set_entry(m,n,hammn);
            if(sym and m!=n) ham.set_entry(n,m,hammn);
        }
    }//m,n
    //std::cerr<<"exit term1\n";
}
#ifdef D3_JLG
void term1(sparse&ham,//3D
           const imatrix&states,
           const int3&params,
           Double ax,
           Double ay,
           const marray&mass,
           const func3&f,
           const int k,
           const int spin,
           const bool sym){
    //std::cerr<<"term1 3"<<std::endl;
    
    const bool mom=(k & 0b1);
    if(mom){
        ax=1/ax;
        ay=1/ay;
    }

    std::vector<Dmatrix> fxyt(params.n(),
                              Dmatrix(params.n(),
                                      Dvector(params.l(),0.)));
    std::vector<Dmatrix> fxyu(params.n(),
                              Dmatrix(params.n(),
                                      Dvector(params.l(),0.)));
#pragma omp parallel for
    for(size_t x=0;x<params.n();x++){
        const Double xx=(params.polytype()?params.x(x):std::sqrt(params.x(x)));
        for(size_t y=0;y<params.n();y++){
            const Double yy=(params.polytype()?
                             params.x(y):std::sqrt(params.x(y)));
            for(size_t t=0;t<params.l();t++){
                const Double tt=params.t(t);
                fxyt[x][y][t]=f(xx/ax,yy/ay,tt,mass);
                fxyu[x][y][t]=f(xx/ax,yy/ay,-tt,mass);
            }//t
        }//y
    }//x
    //std::cerr<<"f"<<std::endl;

    const std::size_t s=states.size();
    ham=sparse(s,s);

#pragma omp parallel for schedule(dynamic)
    for(size_t m=0;m<s;m++)for(size_t n=(sym?m:0);n<s;n++){
        if(!sel1(states,m,n))continue;
        //std::cerr<<m<<','<<n<<'\n';

        //const int L1=states[m][5]*2,L2=states[n][5]*2,
        const int L=states[m][5]*2,
        lp1=states[m][7]*2,
        lp2=states[n][7]*2,
        ll1=states[m][9]*2,
        ll2=states[n][9]*2;
        const int Lmax=std::min(lp1+lp2,ll1+ll2),
        Lmin=std::max(std::abs(lp1-lp2),std::abs(ll1-ll2));

        const Double mnfactor=sqrt(Double((lp1+1)*(lp2+1)*(ll1+1)*(ll2+1))/2)
        *SdotS(spin,states[m][10],states[m][11],states[n][10],states[n][11])
        *parity((L+ll1+ll2)/2,Double(1))
        *(mom?
        parity((states[m][6]+states[m][8]-states[n][6]-states[n][8]
                +(states[m][7]+states[m][9]-states[n][7]-states[n][9])/2))
        :1);
        if(0==mnfactor) continue;
        //std::cerr<<"mnfactor\n";

        Dvector Kfactor((Lmax-Lmin)/2+1,0);
        for(int K=0;K<=(Lmax-Lmin)/2;K++){
            const int Lr=Lmin+K*2;
            Kfactor[K]=W3J(lp1,lp2,Lr,0,0,0)
            *W3J(ll1,ll2,Lr,0,0,0)
            *W6J(lp1,lp2,Lr,ll2,ll1,L)
            *sqrt(Lr+1);
        }
        
        Double tint=0;
        for(std::size_t t=0;t<params.l();t++){
            
            Double Ksum1=0,
            Ksum2=0;
            for(int K=0;K<=(Lmax-Lmin)/2;K++){
                Ksum1+=Kfactor[K]*params.Y(Lmin/2+K,t);
                Ksum2+=Kfactor[K]*params.Y2(Lmin/2+K,t);
            }//K
            if(0==Ksum1 and 0==Ksum2) continue;
            
            Double xint1=0,
            xint2=0;
            for(std::size_t x=0;x<params.n();x++){
                const Double dx=params.dx(x)
                *params.rad(states[m][6],states[m][7],x)
                *params.rad(states[n][6],states[n][7],x);
               
                Double yint1=0,
                yint2=0;
                for(std::size_t y=0;y<params.n();y++){
                    const Double dy=params.dx(y)
                    *params.rad(states[m][8],states[m][9],y)
                    *params.rad(states[n][8],states[n][9],y);
                    
                    if(0!=Ksum1) yint1+=fxyt[x][y][t]*dy;
                    if(0!=Ksum2) yint2+=fxyu[x][y][t]*dy;
                }//y
                if(0!=Ksum1) xint1+=yint1*dx;
                if(0!=Ksum2) xint2+=yint2*dx;
            }//x
            tint+=(Ksum1*xint1+Ksum2*xint2)*params.dt(t);
        }//t
#pragma omp critical(term13)
        {
            ham.set_entry(m,n,tint*mnfactor);
            if(sym and m!=n) ham.set_entry(n,m,tint*mnfactor);
        }
    }//m,n
    //std::cerr<<'\n';
}
#endif//D3_JLG
#ifdef D6_JLG
void term1(CDmatrix&ham,//6D
           const imatrix&states,
           const int6&params,
           Double ax,
           Double ay,
           const marray&mass,
           const func6&f,
           const int k,
           const int spin,
           const bool sym){
    //std::cerr<<"term1 6 "<<f.size()<<'\n';
    
    const bool mom=(k & 0b1);
    if(mom){
        ax=1/ax;
        ay=1/ay;
    }

    const std::size_t st=states.size();
    ham=CDmatrix(st,CDvector(st,0.));

    //if storing function values
    /*CD6 ff(0);
    bool b6=sqr(params.n()*params.l()*params.m())*4<1e8;
    if(b6){
        //Could take up too much space
        //200x200x(50x2)x(50x2)x10x10=4e10 entries
        //50x50x(10x2)x(10x2)x5x5=2.5e7 entries
        //std::cerr<<"creating ff"<<std::endl;
        ff=CD6(params.n(),//x
               CD5(params.n(),//y
                   CD4(2*params.l(),//s
                       CD3(2*params.l(),//t
                           CDmatrix(params.m(),//p
                                    CDvector(params.m(),0))))));//q
        //std::cerr<<"created ff"<<std::endl;
#pragma omp parallel for
        for(size_t x=0;x<params.n();x++){
            Double xx=params.x(x);//GLQx.at(GLQn)
            if(params.polytype()) xx*=xx;
            
            for(size_t y=0;y<params.n();y++){
                //cerr<<m<<','<<n<<' '<<ms<<' '<<x<<','<<y<<'\n';
                Double yy=params.x(y);//GLQx.at(GLQn)
                if(params.polytype()) yy*=yy;
                
                for(size_t s=0;s<2*params.l();s++){//GPQl
                    Double ss=params.t(s/2);//GPQx.at(GPQl)
                    if(s%2)ss=-ss;
                    //symmetric about ss=0
                    
                    for(size_t t=0;t<2*params.l();t++){
                        //std::cerr<<m<<','<<n<<' '<<ms<<' '<<x<<','<<y<<' '<<m1<<','<<m2<<' '<<s<<','<<t<<'\n';
                        Double tt=params.t(t/2);//GPQx.at(GPQl)
                        if(t%2)tt=-tt;
                        //symmetric about tt=0
                        
                        for(size_t p=0;p<params.m();p++){//GPQm
                            Double pp=params.phi(p);//GPQx.at(GPQm)
                            //if(p%2)pp=-pp;
                            //symmetric about pp=0
                            //pp=phix(pp);
                            //symmetric about pp=PI
                            
                            for(size_t q=0;q<params.m();q++){
                                //std::cerr<<'\n'<<m<<','<<n<<' '<<ms<<' '<<x<<','<<y<<' '<<m1<<','<<m2<<' '<<s<<','<<t<<' '<<p<<','<<q<<' ';
                                Double qq=params.phi(q);//GPQx.at(GPQm)
                                //if(q%2)qq=-qq;
                                //qq=phix(qq);
                                ff[x][y][s][t][p][q]=
                                f(xx/ax,yy/ay,ss,tt,pp,qq,mass);
                            }//q
                        }//p
                    }//t
                }//s
            }//y
        }//x
        std::cerr<<"written ff"<<std::endl;
    }
    */

    //std::cerr<<"6int\n";
#pragma omp parallel for schedule(dynamic)
    for(size_t m=0;m<st;m++)for(size_t n=(sym?m:0);n<st;n++){
        if(!sel1(states,m,n)) continue;
#pragma omp critical
        std::cerr<<m<<','<<n<<'\n';
        //std::cerr<<" calc6";
        
        const Double SdS=SdotS(spin,states[m][10],states[m][11],
                               states[n][10],states[n][11]);
        if(0==SdS) continue;

        const int L=states[m][5]*2,
        lp1=states[m][7]*2,
        lp2=states[n][7]*2,
        ll1=states[m][9]*2,
        ll2=states[n][9]*2;
        
        
        const int l1max=std::min(lp1,ll1),
        l2max=std::min(lp2,ll2);
        
        Dvector CG1(l1max+1,0),
        CG2(l2max+1,0);
        for(int m1=0;m1<=l1max;m1++){
            const int mm=2*m1-l1max;
            CG1[m1]=CG(L,0,lp1,mm,ll1,-mm);
        }
        //dispvec(CG1,false,std::cerr);
        for(int m2=0;m2<=l2max;m2++){
            const int mm=2*m2-l2max;
            CG2[m2]=CG(L,0,lp2,mm,ll2,-mm);
        }
        //dispvec(CG2,false,std::cerr);
        
        CDouble sint=0;
        for(std::size_t s=0;s<params.l()*2;s++){
            const Double ss=(s%2?-params.t(s/2):params.t(s/2));
            
            CDouble tint=0;
            for(std::size_t t=0;t<params.l()*2;t++){
                const Double tt=(t%2?-params.t(t/2):params.t(t/2));
                
                CDouble pint=0;
                for(std::size_t p=0;p<params.m();p++){
                    const Double pp=params.phi(p);
                    
                    CDouble qint=0;
                    for(std::size_t q=0;q<params.m();q++){
                        const Double qq=params.phi(q);
                        
                        CDouble m1sum=0;
                        for(int m1=0;m1<=l1max;m1++){
                            if(!CG1[m1]) continue;
                            const int mp1=2*m1-l1max,
                            ml1=l1max-2*m1;
                            const int par1=(0>mp1?parity(mp1):1)
                            *(0>ml1?parity(ml1):1);
                            
                            m1sum+=CG1[m1]*par1
                            *(s%2?
                              params.Y2(lp1/2,std::abs(mp1)/2,s/2)
                              :params.Y(lp1/2,std::abs(mp1)/2,s/2))
                            *(t%2?
                              params.Y2(ll1/2,std::abs(ml1)/2,t/2)
                              :params.Y(ll1/2,std::abs(ml1)/2,t/2))
                            *params.phase(-mp1,p)
                            *params.phase(-ml1,q);
                        }//m1
                        if(CDouble(0)==m1sum) continue;
                        
                        CDouble m2sum=0;
                        for(int m2=0;m2<=l2max;m2++){
                            if(!CG2[m2]) continue;
                            const int mp2=2*m2-l2max,
                            ml2=l2max-2*m2;
                            const int par2=(0>mp2?parity(mp2):1)
                            *(0>ml2?parity(ml2):1);
                            
                            m2sum+=CG2[m2]*par2
                            *(s%2?
                              params.Y2(lp2/2,std::abs(mp2)/2,s/2)
                              :params.Y(lp2/2,std::abs(mp2)/2,s/2))
                            *(t%2?
                              params.Y2(ll2/2,std::abs(ml2)/2,t/2)
                              :params.Y(ll2/2,std::abs(ml2)/2,t/2))
                            *params.phase(mp2,p)
                            *params.phase(ml2,q);
                        }//m2
                        if(CDouble(0)==m2sum) continue;
                        
                        CDouble xint=0;
                        for(std::size_t x=0;x<params.n();x++){
                            const Double xx=(params.polytype()?
                                             params.x(x):
                                             std::sqrt(params.x(x)));
                            
                            CDouble yint=0;
                            for(std::size_t y=0;y<params.n();y++){
                                const Double yy=(params.polytype()?
                                                 params.x(y):
                                                 std::sqrt(params.x(y)));
                                
                                yint+=f(xx/ax,yy/ay,ss,tt,pp,qq,mass)
                                *params.dx(y)
                                *params.rad(states[m][8],states[m][9],y)
                                *params.rad(states[n][8],states[n][9],y);
                            }//y
                            xint+=yint*params.dx(x)
                            *params.rad(states[m][6],states[m][7],x)
                            *params.rad(states[n][6],states[n][7],x);
                        }//x
                        qint+=m1sum*m2sum*xint*params.dphi(q);
                    }//q
                    if(CDouble(0)!=qint) pint+=qint*params.dphi(p);
                }//p
                if(CDouble(0)!=pint) tint+=pint*params.dt(t/2);
            }//t
            if(CDouble(0)!=tint) sint+=tint*params.dt(s/2);
        }//s
        
        /*CDouble xint=0;
        for(size_t x=0;x<params.n();x++){
            const Double xx=(params.polytype()?
                             params.x(x):
                             std::sqrt(params.x(x)));//GLQx.at(GLQn)
            
            CDouble yint=0;
            for(size_t y=0;y<params.n();y++){
                const Double yy=(params.polytype()?
                                 params.x(y):
                                 std::sqrt(params.x(y)));//GLQx.at(GLQn)
                
                CDouble sint=0;
                for(size_t s=0;s<params.l()*2;s++){
                    const Double ss=(s%2?
                                     -params.t(s/2):
                                     params.t(s/2));//GPQx.at(GPQl)
                    
                    CDouble tint=0;
                    for(size_t t=0;t<params.l()*2;t++){
                        const Double tt=(t%2?
                                         -params.t(t/2):
                                         params.t(t/2));//GPQx.at(GPQl)
                        
                        CDouble pint=0;
                        for(size_t p=0;p<params.m();p++){
                            const Double pp=params.phi(p);//GPQx.at(GPQm)
                            
                            CDouble qint=0;
                            for(size_t q=0;q<params.m();q++){
                                const Double qq=params.phi(q);//GPQx.at(GPQm)
                                
                                CDouble m1sum=0;
                                for(int m1=0;m1<=l1max;m1++){
                                    if(!CG1[m1]) continue;
                                    const int mp1=2*m1-l1max,
                                    ml1=l1max-2*m1;
                                    const int par1=(0>mp1?parity(mp1):1)
                                    *(0>ml1?parity(ml1):1);
                                    
                                    m1sum+=CG1[m1]*par1
                                    *(s%2?
                                      params.Y2(lp1/2,std::abs(mp1)/2,s/2)
                                      :params.Y(lp1/2,std::abs(mp1)/2,s/2))
                                    *(t%2?
                                      params.Y2(ll1/2,std::abs(ml1)/2,t/2)
                                      :params.Y(ll1/2,std::abs(ml1)/2,t/2))
                                    *params.phase(-mp1,p)
                                    *params.phase(-ml1,q);
                                }//m1
                                if(CDouble(0)==m1sum) continue;
                                
                                CDouble m2sum=0;
                                for(int m2=0;m2<=l2max;m2++){
                                    if(!CG2[m2]) continue;
                                    const int mp2=2*m2-l2max,
                                    ml2=l2max-2*m2;
                                    const int par2=(0>mp2?parity(mp2):1)
                                    *(0>ml2?parity(ml2):1);
                                    
                                    m2sum+=CG2[m2]*par2
                                    *(s%2?
                                      params.Y2(lp2/2,std::abs(mp2)/2,s/2)
                                      :params.Y(lp2/2,std::abs(mp2)/2,s/2))
                                    *(t%2?
                                      params.Y2(ll2/2,std::abs(ml2)/2,t/2)
                                      :params.Y(ll2/2,std::abs(ml2)/2,t/2))
                                    *params.phase(mp2,p)
                                    *params.phase(ml2,q);
                                }//m2
                                if(CDouble(0)==m2sum) continue;
                                
                                qint+=m1sum*m2sum*params.dphi(q)
                                *f(xx/ax,yy/ay,ss,tt,pp,qq,mass);
                            }//q
                            pint+=qint*params.dphi(p);
                        }//p
                        tint+=pint*params.dt(t/2);
                    }//t
                    sint+=tint*params.dt(s/2);
                }//s
                yint+=sint*params.dx(y)
                *params.rad(states[m][8],states[m][9],y)
                *params.rad(states[n][8],states[n][9],y);
            }//y
            xint+=yint*params.dx(x)
            *params.rad(states[m][6],states[m][7],x)
            *params.rad(states[n][6],states[n][7],x);
        }//x*/
        
        ham[m][n]=sint*SdS
        *(mom?
          std::polar(Double(1),Double(MPI_2)*((states[m][6]+states[m][8]
                                               +states[n][6]+states[n][8])*2
                                              +states[m][7]+states[m][9]
                                              -states[n][7]-states[n][9])):
          CDouble(1));
        if(sym and m!=n) ham[n][m]=ham[n][m];
    }//m,n
    //std::cerr<<'\n';
    //std::cerr<<"term1 6 end\n";
}
#endif//D6_JLG

void term5(sparse&ham,//1D
           const imatrix&states,
           //const int1&params,
           Double a,
           const marray&mass,
           const func1&f,
           const int k,
           const int spin,
           const bool sym){
    //std::cerr<<"term5 1"<<std::endl;
    
    const bool mom=(k & 0b1);
    if(mom) a=1/a;
    const Double f0=f(0,mass);
    //std::cerr<<"fx\n";

    const std::size_t s=states.size();
    ham=sparse(s,s);

#pragma omp parallel for schedule(dynamic)
    for(size_t m=0;m<s;m++)for(size_t n=(sym?m:0);n<s;n++){
        //std::cerr<<m<<' '<<n<<'\n';
        if(!sel5(states,m,n)) continue;
        //if(states[m][11]!=states[n][11]) continue;//Ps

        const Double SdS=SdotS(spin,states[m][10],states[m][11],
                               states[n][10],states[n][11]);
        if(0==SdS) continue;
        
        Double h=f0*SdS*(a*a*a)*M8_PI*
        std::exp((std::lgamma(1.5+states[m][6])-
                  std::lgamma(1+states[m][6])+
                  std::lgamma(1.5+states[n][6])-
                  std::lgamma(1+states[n][6]))/2);
        
        if(mom) h*=parity(states[m][6]-states[n][6]);
#pragma omp critical(term5)
        {
            ham.set_entry(m,n,h);
            if(sym and m!=n) ham.set_entry(n,m,h);
        }
    }//m,n
}

//******************************************************************************
//term2
//S2ten term unnecessary, math redone so that it would return +-1

inline bool sel2(const imatrix&states,
                 const std::size_t m,
                 const std::size_t n){
    if(states[m][0]!=states[n][0]) return false;//Pf
    if(states[m][1]!=states[n][1]) return false;//J
    if(states[m][2]!=states[n][2]) return false;//Pj
    if(states[m][3]!=states[n][3]) return false;//M
    if(!istriangle(states[m][5]*2,states[n][5]*2,4)) return false;//L triangle
    if(!istriangle(states[m][7]*2,states[n][7]*2,4)) return false;//lr triangle
    if(states[m][8]!=states[n][8]) return false;//nl
    if(states[m][9]!=states[n][9]) return false;//ll
    if(3!=states[m][10] and 3!=states[n][10])return false;//S b/c S'=3 or S=3
    if(1==states[m][11] or 1==states[n][11]) return false;//Ps b/c rho=0
    return true;
}

void term2(sparse&ham,//1D spin tensor
           const imatrix&states,
           const int1&params,
           const Double a,
           const marray&mass,
           const func1&f,
           const bool sym){
    //std::cerr<<"term2 1"<<std::endl;
    //int N=states.back()[4];
    //N=2(np+nl)+lp+ll so max possible np is N/2 if even or (N-1)/2 if odd
    //So N/2 with integer division works either way
    //Uses GHQ

 	//const Double mp=mass[0]*mass[1]/(mass[0]+mass[1]);
 	//const Double ml=(mass[0]+mass[1])*mass[2]/(mass[2]+mass[0]+mass[1]);

    //const Double a=omega*mp;

    Dvector fx(params.n());
#pragma omp parallel for
    for(size_t x=0;x<params.n();x++){
        const Double xx=(params.polytype()?params.x(x):std::sqrt(params.x(x)));
        fx[x]=f(xx/a,mass);
    }

    std::size_t s=states.size();
    ham=sparse(s,s);

    //std::cerr<<"start int "<<s<<std::endl;
#pragma omp parallel for schedule(dynamic)
    for(size_t m=0;m<s;m++)for(size_t n=(sym?m:0);n<s;n++){
        if(!sel2(states,m,n)) continue;

//#pragma omp critical(disp1)
        //std::cerr<<10000+100*m+n<<std::endl;
        Double mnfactor=
        parity(states[m][7]+states[n][9]+(states[m][1]-1)/2,
        //(-1)^(J+S'+l_rho'+l_lambda)<S'||sqrt(6/5)S_2||S>
        //(-1)^(J+S')<S'||sqrt(6/5)S_2||S> == (-1)^(J-1/2)
        std::sqrt(5*(2*states[m][5]+1)*(2*states[n][5]+1)*(2*states[n][7]+1)));
        //sqrt(5/6)sqrt(6(2L'+1)(2L+1)(2l_rho+1))
        if(0.==mnfactor) continue;
        mnfactor*=CG(2*states[m][7],0,2*states[n][7],0,4,0);//lr',lr,2
        if(0.==mnfactor) continue;
        mnfactor*=W6J(2*states[m][5],2*states[n][5],4,//L',L,2
                      states[n][10],states[m][10],states[n][1]);//S,S',J
        if(0.==mnfactor) continue;
        mnfactor*=W6J(2*states[n][7],2*states[m][7],4,//lr,lr',2
                      2*states[m][5],2*states[n][5],2*states[m][9]);//L',L,ll
        if(0.==mnfactor) continue;

//#pragma omp critical(disp2)
        //std::cerr<<20000+100*m+n<<std::endl;
        Double hmn=0.;
        for(size_t x=0;x<params.n();x++){
            hmn+=params.dx(x)*fx[x]
            *params.rad(states[m][6],states[m][7],x)
            *params.rad(states[n][6],states[n][7],x);
        }//x
//#pragma omp critical(disp3)
        //std::cerr<<30000+100*m+n<<std::endl;
#pragma omp critical(term2)
        {
            ham.set_entry(m,n,hmn*mnfactor);
            if(sym and m!=n) ham.set_entry(n,m,hmn*mnfactor);
        }
//#pragma omp critical(disp4)
        //std::cerr<<40000+100*m+n<<std::endl;
    }//m,n
    //std::cerr<<states[5][1]<<' '<<states[2][1]<<' '<<ham[5][2]<<std::endl;
}
//No 3D or 6D terms because only dealing with rho integral.
//Other pairwise done through Moshinski transformation

//******************************************************************************
//term3 and term4

inline bool sel3(const imatrix&states,
                 const std::size_t m,
                 const std::size_t n){
    if(states[m][0]!=states[n][0]) return false;//Pf
    if(states[m][1]!=states[n][1]) return false;//J
    if(states[m][2]!=states[n][2]) return false;//Pj
    if(states[m][3]!=states[n][3]) return false;//M
    if(!istriangle(states[m][5]*2,states[n][5]*2,2)) return false;//L triangle
    if(states[m][7]!=states[n][7] or 0==states[m][7]) return false;//lr
    if(states[m][8]!=states[n][8]) return false;//nl
    if(states[m][9]!=states[n][9]) return false;//ll
    return true;
}

inline bool sel4(const imatrix&states,
                 const std::size_t m,
                 const std::size_t n){
    if(states[m][0]!=states[n][0]) return false;//Pf
    if(states[m][1]!=states[n][1]) return false;//J
    if(states[m][2]!=states[n][2]) return false;//Pj
    if(states[m][3]!=states[n][3]) return false;//M
    if(istriangle(states[m][5]*2,states[n][5]*2,2)) return false;//L triangle
    if(std::abs(states[m][7]-states[n][7])!=1) return false;//lr'=lr+-1
    if(std::abs(states[m][9]-states[n][9])!=1) return false;//ll'=ll+-1
    return true;
}

inline Double Svec(const unsigned int s1,const unsigned int p1,
                   const unsigned int s2,const unsigned int p2,
                   const marray&pars){
    const Double a=pars[0];
    const Double b=pars[1];
    const Double g=pars[2];

    if((1!=s1 and 3!=s1) or (1!=s2 and 3!=s2)){
#ifdef TESTING_STATES_JLG
        std::cerr<<"not right spin: "<<s1<<' '<<s2<<std::endl;
        throw std::domain_error("Svec spins");
#endif
        return 0.;
    }

    switch (1000*s1+100*p1+10*s2+p2) {
        case 3030:
            return Ms5_s3*(a+b+g);
        case 3011:
            return b-a;
        case 3010:
            return M1_s3*(2*g-a-b);
        case 1130:
            return a-b;
        case 1111:
            return Ms3_s2*g;
        case 1110:
            return M1_s2*(b-a);
        case 1030:
            return M1_s3*(a+b-2*g);
        case 1011:
            return M1_s2*(b-a);
        case 1010:
            return Ms2_s3*(a+b-0.5*g);
        default:
#ifdef TESTING_STATES_JLG
            std::cerr<<"not right spin\n"<<
            s1<<','<<p1<<','<<s2<<','<<p2<<','<<a<<','<<b<<','<<g<<std::endl;
            throw std::domain_error("Svec spins");
#endif
            return 0.;
    }
}

void term3(sparse&ham,
           const imatrix&states,
           const int1&params,
           const Double a,
           const marray&mass,
           const func1&f,
           const marray&pars){
    //const Double mp=mass[0]*mass[1]/(mass[0]+mass[1]);
    //const Double ml=(mass[0]+mass[1])*mass[2]/(mass[2]+mass[0]+mass[1]);
    //const Double ax=omega*mp;
    //const Double ay=omega*ml;

    Dvector fx(params.n());
#pragma omp parallel for
    for(size_t x=0;x<params.n();x++){
        const Double xx=(params.polytype()?params.x(x):std::sqrt(params.x(x)));
        fx[x]=f(xx/a,mass)*a/xx;
    }

    std::size_t s=states.size();
    ham=sparse(s,s);

#pragma omp parallel for collapse(2) schedule(dynamic)
    for(size_t m=0;m<s;m++)for(size_t n=0;n<s;n++){
        //std::cerr<<m<<' '<<n<<'\n';
        if(!sel3(states,m,n)) continue;

        const int L1=states[m][5]*2;
        const int L2=states[n][5]*2;
        const int lp=states[m][7]*2;
        //const int lp2=states[n][7]*2;
        const int ll=states[m][9]*2;
        //const int ll2=states[n][9]*2;

        const Double mnfactor=Svec(states[m][10],states[m][11],
                                   states[n][10],states[n][11],pars)
        *W6J(L1,L2,2,
             states[n][10],states[m][10],states[m][1])//S,S',J
        *W6J(lp,lp,2,
             L1,L2,ll)
        *std::sqrt((L1+1)*(L2+1)*lp*(lp+2)*(lp+1)/4)
        *parity((states[n][1]+states[m][10]+lp+ll)/2+1,1.);
        if(0.==mnfactor) continue;

        Double hmn=0.;
        for(size_t x=0;x<params.n();x++){
            hmn+=params.dx(x)*fx[x]
            *params.rad(states[m][6],states[m][7],x)
            *params.rad(states[n][6],states[n][7],x);
        }//x
#pragma omp critical(term3)
        ham.set_entry(m,n,hmn*mnfactor);
    }//m,n
}

void term4(sparse&ham,
           const imatrix&states,
           const int1&params,
           const Double ax,
           const Double ay,
           const marray&mass,
           const func1&f,
           const marray&pars){
    const std::size_t s=ham.size();

    //const Double mp=mass[0]*mass[1]/(mass[0]+mass[1]);
    //const Double ml=(mass[0]+mass[1])*mass[2]/(mass[2]+mass[0]+mass[1]);
    //const Double ax=omega_rho*mp;
    //const Double ay=omega_lambda*ml;

    Dvector fx(params.n());
#pragma omp parallel for
    for(size_t x=0;x<params.n();x++){
        Double xx=params.x(x);//GLQx.at(GLQn)
        if(!params.polytype()) xx=std::sqrt(xx);
        fx[x]=f(xx/ax,mass);
    }

    static std::map<std::array<int,3>,Double> pfactor;
#pragma omp parallel for collapse(2) schedule(dynamic)
    for(size_t m=0;m<s;m++)for(size_t n=0;n<=s;n++){
        if(!sel4(states,m,n)) continue;

        const std::array<int,3> pf({states[m][9],states[m][8],states[n][8]});
        bool pbool=false;
#pragma omp critical(pfactor)
        {
            if(0==pfactor.count(pf)) pfactor[pf]=0.;
            else pbool=true;
        }
        if(pbool) continue;

        Double pfac=0.;
        for(size_t x=0;x<params.n();x++){
            const Double xx=(params.polytype()?
                             params.x(x):std::sqrt(params.x(x)));
            pfac+=params.dx(x)*xx
            *params.rad(states[m][8],states[m][9],x)
            *params.rad(states[n][8],states[n][9],x);;
        }//x
#pragma omp critical(pfactor2)
        pfactor[pf]=pfac;//*std::sqrt(ay);
        //alpha factor stays out of pfactor so that it is independent of scale
        //requires ay factor in integral
    }//m,n

#pragma omp parallel for collapse(2) schedule(dynamic)
    for(size_t m=0;m<s;m++)for(size_t n=0;n<s;n++){
        if(!sel4(states,m,n)) continue;

        if(0==pfactor[{states[m][9],states[m][8],states[n][8]}]) continue;

        const int L1=states[m][5]*2;
        const int L2=states[n][5]*2;
        const int lp1=states[m][7]*2;
        const int lp2=states[n][7]*2;
        const int ll1=states[m][9]*2;
        const int ll2=states[n][9]*2;

        const Double mnfactor=ay
        *pfactor[{states[m][9],states[m][8],states[n][8]}]
        *Svec(states[m][10],states[m][11],
              states[n][10],states[n][11],pars)
        *W6J(L1,L2,2,//L',L,1
             states[n][10],states[m][10],states[m][1])//S,S',J
        *W9J(2,2,2,
             lp1,ll1,L1,
             lp2,ll2,L2)
        *std::sqrt(6*(L1+1)*(L2+1)*(ll1>ll2?ll1:ll2)*(lp1>lp2?lp1:lp2)/4)
        *parity((states[m][1]+states[m][10]+L2)/2+states[m][8]+states[n][8]+(lp2>lp1),1.);
        if(0.==mnfactor) continue;

        Double hmn=0.;
        for(size_t x=0;x<params.n();x++){
            hmn+=params.dx(x)*fx[x]
            *params.rad(states[m][6],states[m][7],x)
            *params.rad(states[n][6],states[n][7],x);
        }//x
#pragma omp critical(term4)
        ham.set_entry(m,n,hmn*mnfactor);
    }//m,n
}
