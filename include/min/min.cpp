//
//  min.cpp
//
//
//  Created by Johnathan Gross
//
//

#include "min.h"
#include <iostream>
#include <stdexcept>
#include <string>

inline constexpr long double GOLD=1.61803398874989484820458683436563811772030917980576286213544862270526046281890244970720720418939113748475408807538L;
inline constexpr long double CGOLD=0.38196601125010515179541316563436188227969082019423713786455137729473953718109755029279279581060886251524591192461L;

inline constexpr Double EPS=std::numeric_limits<Double>::epsilon();

template<class T>
constexpr T sqr(T x){
    static_assert(std::is_arithmetic<T>::value);
    return x*x;
}

/*
arma::sp_mat sparse2arma(const sparse&S){
    arma::umat loc(2,S.size());
    arma::vec val(S.size());
    auto it=S.cbegin();
    for(arma::uword count=0;count<S.size();count++){
        loc(0,count)=it->first[0];
        loc(1,count)=it->first[1];
        val(count)=it->second;
        it++;
    }
    return arma::sp_mat(loc,val,S.m(),S.n(),true,false);
}

sparse arma2sparse(const arma::sp_mat&S){
    sparse M(S.n_rows,S.n_cols);
    for(auto it=S.begin();it!=S.end();++it){
        M.set_entry(it.row(),it.col(),*it);
    }
    return M;
}
 */

std::array<Double,6> bracket(const std::function<Double(const Double)>&func,
                             Double xa,
                             Double xb,
                             bool pos,
                             Double fa,
                             Double fb){
    //std::cerr<<"bracket"<<std::endl;
    if(pos){
#ifndef TESTING_STATES_JLG
        if(0>xa) xa=1e-3;
        if(0>xb) xb=2e-3;
#else
        if(0>xa) throw std::invalid_argument("negative guess");
        if(0>xb) throw std::invalid_argument("negative guess");
#endif
    }
    if(std::isnan(fa)) fa=func(xa);
    if(std::isnan(fb)) fb=func(xb);
    //std::cerr<<"f"<<std::endl;
    if(fb>fa){//swaps so always moving in direction from a to b
        std::swap(xa,xb);
        std::swap(fa,fb);
    }
    Double xc=xb+GOLD*(xb-xa);
    if(pos and 0>xc){
        std::cerr<<"c="<<xc<<std::endl;
        xc=1e-3;
    }
    Double fc=func(xc);

    size_t steps=0;
    while(fb>fc){
        //std::cerr<<"step "<<steps<<std::endl;
        Double s=(fa-fb)*(xc-xb),
        ss=(fa-fb)*sqr(xc-xb),
        t=(fb-fc)*(xb-xa),
        tt=(fb-fc)*sqr(xb-xa);
        Double u=xb+(ss+tt)/
        (std::abs(s-t)>1e-10?2*(s-t):std::copysign(2e-10,s-t)),
        //avoid divide by 0
        ulim=xb+10*GOLD*(xc-xb),//parabolic extrapolation only valid near min
        fu;
        if(pos and 0>ulim) ulim=std::min(std::min(xb,xc)*1e-3,std::abs(ulim));
        steps++;

        if(pos and 0>u){
            //std::cerr<<"xa="<<xa<<" xb="<<xb<<" xc="<<xc
            //<<" fa="<<fa<<" fb="<<fb<<" fc="<<fc
            //<<" u="<<u<<std::endl;
            u=std::min(std::min(xb,xc)*1e-3,std::abs(u));
        }
        if((xb-u)*(u-xc)>0){//u between b and c
            fu=func(u);
            if(fu<fc){//b,u,c bracket minimum
                xa=xb;
                fa=fb;
                xb=u;
                fb=fu;
                break;
            }else if(fu>fb){//a,b,u bracket local minimum
                //fc<fb, so there is an even lower value
                //good enough for *a* minimum
                //should be unnecessary if only one minimum

                //std::cerr<<"recursive bracket "<<steps<<std::endl;
                //return bracket(func,u,xc,pos,fu,fc);
                xc=u;
                fc=fu;
                break;
            }else{//fc<=fu<=fb
                u=xc+GOLD*(xc-xb);
                fu=func(u);
            }
        }else if((xc-u)*(u-ulim)>0){//u between c and ulim
            fu=func(u);
            if(fu<fc){
                xb=xc;
                fb=fc;
                xc=u;
                fc=fu;

                u=u+GOLD*(u-xb);
                if(pos and 0>u){
                    std::cerr<<"u="<<u<<std::endl;
                    u=std::min(std::min(xc,ulim)*1e-3,std::abs(u));
                }
                fu=func(u);
            }//else b,c,u bracket
        }else if((u-ulim)*(ulim-xc)>=0){//u beyond ulim
            u=ulim;
            fu=func(u);
        }else{//u uphill from b, u=b, u=c, reject parabola
            u=xc+GOLD*(xc-xb);
            if(pos and 0>u){
                std::cerr<<"u="<<u<<std::endl;
                u=std::min(std::min(xb,xc)*1e-3,std::abs(u));
            }
            fu=func(u);
        }
        xa=xb;
        fa=fb;
        xb=xc;
        fb=fc;
        xc=u;
        fc=fu;
    }

    if(xa>xc){
        std::swap(xa,xc);
        std::swap(fa,fc);
    }
    //std::cerr<<steps<<std::endl;
    //std::cerr<<"bracket end"<<std::endl;
    return {xa,xb,xc,fa,fb,fc};
}

std::array<Double,6> Brentmin(const std::function<Double(const Double)>&func,
                              Double a,
                              Double b,
                              Double prec,
                              bool rel,
                              bool pos,
                              bool disp,
                              Double dx,
                              Double fa,
                              Double fb){
    //std::cerr<<"Brentmin"<<std::endl;
    auto br=bracket(func,a,b,pos,fa,fb);
    //std::cerr<<"br "<<br[0]<<' '<<br[3]<<' '<<br[1]<<' '<<br[4]
    //<<' '<<br[2]<<' '<<br[5]<<std::endl;
    a=br[0];
    fa=br[3];
    b=br[2];
    fb=br[5];
    Double x,w,v,u,u1,fx,fw,fv,fu,fu1;
    x=w=v=br[1];//lowest value,second lowest value,previous w
    u=u1=br[2];//current evaluation point,previous u
    fx=fw=fv=br[4];
    fu=fu1=br[5];
    Double d=0,e=0;//last step and next to last step
    Double a1=a,b1=b;//previous a and b
    std::size_t steps=0,eqcount=0;
    std::array<Double,6> store{a,x,b,fa,fx,fb};
    do{
        //std::cerr<<steps<<' '<<x<<' '<<fx<<' '
        //<<(w-x)<<' '<<(fw-fx)<<std::endl;
//#ifdef TESTING_STATES_JLG
        //if(100<steps) throw std::logic_error("too long");
//#endif

    loopstart:
        Double xm=(a+b)/2;
        Double dx1=dx*std::abs(x)+EPS*1e-3;
        Double dx2=2*dx1;
        Double dy=prec*(rel?std::max(std::abs(fx),EPS):1);
//#ifdef TESTING_STATES_JLG
         //std::cerr<<x<<' '<<std::abs(x-xm)+(b-a)/2<<' '<<dx2<<' '
         //<<fx<<' '<<std::max(fa,fb)-fx<<' '<<dy<<'\n';
//#endif
        if(std::max(fa,fb)<=(dy+fx) or
           //check if within relative/absolute precision in y
           std::abs(x-xm)<=(dx2-(b-a)/2)){
            //check if within relative machine precision in x
            //(b-a)<=dx2){
            //check if bounds close together
            //std::cerr<<steps<<std::endl;
            if(disp) std::cerr<<"Brentmin end"<<std::endl;
            if(eqcount) std::cerr<<"end eq "<<eqcount<<std::endl;
            return {a,x,b,fa,fx,fb};
        }

        if(eqcount) std::cerr<<"e "<<eqcount<<std::endl;
        switch(eqcount){
            case 0:
                break;
            case 1:
                u1=u;
                fu1=fu;
                a1=a;
                b1=b;
                steps++;
                if(disp) std::cerr<<steps<<':';

                u=u-2*d;
                if(u<=a) u=(a+x)/2;
                if(u>=b) u=(b+x)/2;
                e=d;
                d=u-x;
                if(disp) std::cerr<<u<<' '<<d<<' '<<u-u1<<std::endl;
                goto funceval;
            case 2:{
                Double fxx=func(x);
                std::cerr<<a<<' '<<x<<' '<<b<<' '<<u<<':'
                <<fa<<' '<<fb<<' '<<fxx<<' '<<fx<<' '<<fxx-fx<<' '<<dy<<std::endl;
                if(fxx!=fx){
                    store={a,x,b,fa,fx,fb};
                    br=bracket(func,a1,b1,pos);
                    a=br[0];
                    fa=br[3];
                    b=br[2];
                    fb=br[5];
                    x=w=v=br[1];
                    fx=fw=fv=br[4];
                    fu=fu1=br[5];
                    a1=a;
                    b1=b;
                    eqcount++;
                    goto loopstart;
                }else{
                    std::cerr<<a<<' '<<x<<' '<<b<<' '
                    <<fa<<' '<<fx<<' '<<fb<<std::endl;
                    return {a,x,b,fa,fx,fb};
                }
            }
            case 3:
                for(auto i:store) std::cerr<<i<<' ';
                std::cerr<<'\n'<<a<<' '<<x<<' '<<b<<' '
                <<fa<<' '<<fx<<' '<<fb<<std::endl;
                if(x!=store[1]){
                    d=e=0;
                    eqcount=0;
                    break;
                }else{
                    if(a>store[0]){
                        a=store[0];
                        fa=store[3];
                    }
                    if(b<store[2]){
                        b=store[2];
                        fb=store[5];
                    }
                    if(fx>store[4]) fx=store[4];
                    std::cerr<<"can't get more accurate"<<std::endl;
                    return {a,x,b,fa,fx,fb};
                }
            default:
                throw std::logic_error("eqcount out of range");
        }

        u1=u;
        fu1=fu;
        a1=a;
        b1=b;
        steps++;
        if(disp) std::cerr<<steps<<':';
        if(std::abs(e)>dx1){
            //next to last step bigger than error
            Double s=(fx-fw)*(x-v),
            ss=(fx-fw)*sqr(x-v),
            t=(fx-fv)*(x-w),
            tt=(fx-fv)*sqr(x-w),
            temp=e;
            e=d;
            Double p=(ss-tt),q=(2*(s-t));
            if(0<q)p=-p;
            else q=-q;

            if(std::abs(p)>=std::abs(q*temp/2) or
               //step bigger than second to last step/2
               p<=q*(a-x) or p>=q*(b-x)){
                //step goes below a or above b
                d=CGOLD*(e=(x>=xm?a-x:b-x));
            }else{
                d=p/q;
                u=x+d;
                if((u-a)<dx2 or (b-u)<dx2)
                    //u close to a or b
                    d=std::copysign(dx1,xm-x);
            }
        }else{
            //next to last step was small
            d=CGOLD*(e=(x>=xm?a-x:b-x));
            
        }
        u=(std::abs(d)>=dx1?x+d:x+std::copysign(dx1,d));
        d=u-x;
        //u=x+d;
        //if step is small, take minimum step
        if(disp) std::cerr<<u<<' '<<d<<' '<<u-u1<<std::endl;
        if(u==u1){
            u=(u+x)/2;
            d=u-x;
            std::cerr<<u1<<' '<<u<<' '<<d<<std::endl;
        }
        if(u==x) throw std::logic_error("u==x");
        if(u<a or b<u){
            std::cerr<<a<<' '<<x<<' '<<b<<' '<<u<<' '<<u1<<' '
            <<u-a<<' '<<b-u<<'\n'
            <<x-xm<<' '<<u-x<<' '<<d<<' '<<dx<<' '<<dx1<<' '<<dx2<<'\n'
            <<fa<<' '<<fx<<' '<<fb<<' '<<fu1<<' '<<dy<<std::endl;
            //throw std::logic_error("out of bounds");
            throw std::logic_error("out of bounds");
            //return {a,x,b,fa,fx,fb};
        }

    funceval:
        try{
            fu=func(u);
        }catch(std::exception&exp){
            std::cerr<<x<<' '<<a-x<<' '<<b-x<<' '<<w-x<<' '<<v-x<<' '<<u-x<<'\n';
            std::cerr<<dx<<' '<<dx1<<' '<<dx2<<' '<<d<<' '<<e<<'\n';
            std::cerr<<fx<<' '<<fa-fx<<' '<<fb-fx<<' '<<fw-fx<<' '<<fv-fx<<' '<<fv-fw<<' '<<fu-fx<<'\n';
            std::cerr<<prec<<' '<<dy<<std::endl;
            std::cerr<<exp.what()<<std::endl;
            throw exp;
        }catch(...){
            std::cerr<<"fu error"<<std::endl;
            throw;
        }
        if(disp) std::cerr<<fu<<' '<<fu-fx<<' '<<fu-fu1<<std::endl;

        if(fu==fu1) {
            eqcount++;
            std::cerr<<"equal output "<<eqcount<<':'
            <<x<<' '<<u<<' '<<fu<<std::endl;
        }else eqcount=0;
        if(fu<=fx){
            if(u>=x){
                a=x;
                fa=fx;
            }else{
                b=x;
                fb=fx;
            }
            v=w;
            fv=fw;
            w=x;
            fw=fx;
            x=u;
            fx=fu;
            //std::cerr<<"min "<<x<<' '<<fx<<std::endl;
        }else{//fu>fx
            if(u<x){
                a=u;
                fa=fu;
            }else{
                b=u;
                fb=fu;
            }
            if(fu<=fw or w==x){
                v=w;
                fv=fw;
                w=u;
                fw=fu;
            }else if(fu<=fv or v==x or v==w){
                v=u;
                fv=fu;
            }//else fu>fv
        }
        if(disp) std::cerr<<a<<' '<<x<<' '<<b<<' '<<b-a<<std::endl;
        if(a1==a and b1==b) throw std::logic_error("bounds not moved");
        if(a<a1) throw std::logic_error("a lower");
        if(b>b1) throw std::logic_error("b higher");
        if(pos and !x){
            std::cerr<<a<<' '<<x<<' '<<b<<' '
            <<u<<' '<<u1<<' '<<v<<' '<<w<<'\n'
            <<fa<<' '<<fx<<' '<<fb<<' '
            <<fu<<' '<<fu1<<' '<<fv<<' '<<fw<<'\n'
            <<d<<' '<<dx<<' '<<dx1<<' '<<dx2<<' '<<dy<<std::endl;
            throw std::logic_error("x is 0");
        }
    }while(true);
    std::cerr<<"Brentmin loop end"<<std::endl;
    throw std::runtime_error("out of loop");
}
