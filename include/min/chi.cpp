//
//  min.cpp
//
//
//  Created by Johnathan Gross
//
//

#include "chi.h"
#include "min.h"
#include "mout.h"
#include <iostream>
#include <stdexcept>

inline constexpr Double EPS=std::numeric_limits<Double>::epsilon();

template<class T>
constexpr T sqr(T x){
    static_assert(std::is_arithmetic<T>::value);
    return x*x;
}

Double chisquared(const Dvector&model,
                  const Dvector&data,
                  const Dvector&errors){
    /*Calculates the chi^2 of model to a set of data values with errors
     */
    std::size_t s=model.size();
#ifdef TESTING_STATES_JLG
    if(data.size()!=s){
        std::cerr<<"model: "<<s<<" data: "<<data.size()<<std::endl;
        throw std::invalid_argument("data not same size as model");
    }
    if(errors.size()!=s){
        std::cerr<<"model: "<<s<<" errors: "<<errors.size()<<std::endl;
        throw std::invalid_argument("errors not same size as model");
    }
#endif
    Double x=0;
#pragma omp parallel for reduction(+:x)
    for(size_t i=0;i<s;i++)
        x+=sqr((model[i]-data[i])/errors[i]);
    return x;
}
Double chisquared(const Dvector&model,
                  const Dvector&data,
                  const Double error){
    std::size_t s=model.size();
#ifdef TESTING_STATES_JLG
    if(data.size()!=s){
        std::cerr<<"model: "<<s<<" data: "<<data.size()<<std::endl;
        throw std::invalid_argument("data not same size as model");
    }
#endif
    Double x=0;
#pragma omp parallel for reduction(+:x)
    for(size_t i=0;i<s;i++)
        x+=sqr(model[i]-data[i]);
    return x/sqr(error);
}

Double deltaC(const Dvector&o,
              const Dvector&d,
              const Dvector&e){
    Double s=o.size();
#ifdef TESTING_STATES_JLG
    if(d.size()!=s) throw std::invalid_argument("data and out not same size");
    if(e.size()!=s) throw std::invalid_argument("errors and out not same size");
#endif
    
    Double A=0,B=0;
#pragma omp parallel for reduction(+:A,B)
    for(size_t i=0;i<d.size();i++){
        A+=std::pow(e[i],-2);
        B+=(o[i]-d[i])/sqr(e[i]);
    }
    return -B/A;
}

void newC(Dvector&o,
          Double&C,
          const Dvector&d,
          const Dvector&e){
    const Double c2=deltaC(o,d,e);
    C+=c2;
#pragma omp parallel for
    for(size_t i=0;i<o.size();i++) o[i]+=c2;
}

inline auto lmin(const std::function<Double(const Double)>&f,
                 const Double error,
                 const Double rel,
                 const bool disp=false){
    //std::cerr<<"lmin"<<std::endl;
    return Brentmin(f,0,1,error,rel,false,disp);
}

chimin::chimin(const cfunc3&f,
               const Dvector&d,
               const Dvector&s,
               const Dvector&e,
               const Dvector&i,
               const Dvector&v,
               const Double oer,
               const Double ier,
               const bool of,
               const bool r):
data(d),sig(s),eparam(e),iparam(i),vals(v),
nd(data.size()),ne(eparam.size()),neo(eparam.size()-of),ni(iparam.size()),
df(data.size()-eparam.size()),
oerror(oer),ierror(ier),off(of),rel(r),
func(f){
    if(ne>nd){
#ifdef TESTING_STATES_JLG
        throw std::invalid_argument("overfitting");
#else
        df=0;
#endif
    }
    if(1==sig.size()) sig=Dvector(nd,s[0]);
#ifdef TESTING_STATES_JLG
    else if(sig.size()!=nd)
        throw std::invalid_argument("error not same size as data");
#endif
    if(0==vals.size()) vals=Dvector(nd,0);
    else if(vals.size()!=nd){
#ifdef TESTING_STATES_JLG
        throw std::invalid_argument("vals not same size as data");
#else
        vals=Dvector(nd,0);
#endif
    }
    chisq=chisquared(vals,data,sig);
}

chimin::chimin(const cfunc2&f,
               const Dvector&d,
               const Dvector&s,
               const Dvector&e,
               const Dvector&i,
               const Dvector&v,
               const Double oer,
               const Double ier,
               const bool of,
               const bool r):
chimin(c2to3(f),d,s,e,i,v,oer,ier,of,r){
#ifdef TESTING_STATES_JLG
    if(ni)
        throw std::invalid_argument("funciton doesn't use internal parameters");
#endif
}

Double chimin::eval(const Dvector&e,Dvector&i,Dvector&o){
    func(e,i,o);
    return chisquared(o,data,sig);
}

bool is_square(const Dmatrix&M){
    const size_t ne=M.size();
    bool b=true;
#pragma omp parallel for reduction(&&:b)
    for(std::size_t i=0;i<ne;i++) b=M.size()==ne;
    return b;
}
bool mod_all(const Dmatrix&M,bvector&V){
    const size_t n=M.size();
    V=bvector(n,true);
    if(!is_square(M)) return false;
    bool b=true;
#pragma omp parallel for reduction(&&:b)
    for(size_t j=0;j<n;j++){
        bool bb=false;
        for(size_t i=0;i<n;i++) if(M[i][j]){
            bb=true;
            break;
        }
        if(!bb){
#pragma omp critical(Vbool)
            V[j]=false;
            b=false;
        }
    }
    return b;
}

void chimin::minimize(Dvector&P,
                      Dmatrix&xmat,
                      bool disp){
    //std::cerr<<"chimin::minimize"<<std::endl;
#ifdef TESTING_STATES_JLG
    if(P.size()!=ne) throw std::invalid_argument("P not size of parameters");
#endif
    if(0==xmat.size()){
        xmat=Dmatrix(neo,Dvector(neo,0));
#pragma omp parallel for
        for(std::size_t i=0;i<neo;i++){
            if(P[i]) xmat[i][i]=P[i]*.1;
            else xmat[i][i]=.1;
        }
    }
#ifdef TESTING_STATES_JLG
    else if(xmat.size()!=neo)
        throw std::invalid_argument("xmat not size of parameters");
    else{//xmat.size()==neo
        bvector ex;
        if(!mod_all(xmat,ex)){
            if(ex==bvector(neo,true))
                throw std::invalid_argument("xmat not square");
            else{
                dispvec(ex,false,std::cerr);
                throw std::invalid_argument("xmat doesn't modify");
            }
        }
    }
#endif
    //std::cerr<<"chi minimize steup\n";
    
    auto ff=[this](const Dvector&pp,
                         Dvector&in){
        /*Calculates func(|pp>) and returns chisquared
         */
        //std::cerr<<"ff"<<std::endl;
        Dvector out;
        try{
            func(pp,in,out);
        }catch(std::invalid_argument&e){
            std::cerr<<"invalid argument function eval ff: "<<e.what()<<'\n';
            dispvec(pp,false,std::cerr);
            throw;
        }catch(std::logic_error&e){
            std::cerr<<"logic error func eval ff: "<<e.what()<<std::endl;
            throw;
        }catch(...){
            std::cerr<<"other func eval ff"<<std::endl;
            throw;
        }
        return chisquared(out,data,sig);
    };

    auto fx=[this,&ff](const Dvector&pp,
               const Dvector&x,
               Dvector&in){
        /*Returns function that calculates func(|pp>+y*|x>) and returns
         chisquared
         Effect: stores out of func in o
         */
        //std::cerr<<"fx"<<std::endl;
        return [this,&pp,&x,&in,&ff](const Double y){
            //std::cerr<<"fx lambda"<<std::endl;
            Dvector z(pp);
            if(0!=y){
#pragma omp parallel for
                for(size_t i=0;i<neo;i++){
                    z[i]+=y*x[i];
                    //if(0>z[i]) z[i]=1e-3;
                }
            }
            try{
                return ff(z,in);
            }catch(std::invalid_argument&e){
                std::cerr<<"invalid argument ff in fx: "<<e.what()<<'\n';
                dispvec(z,false,std::cerr);
                throw std::invalid_argument("ia ff in fx");
            }catch(std::logic_error&e){
                std::cerr<<"logic error ff in fx: "<<e.what()<<'\n';
                dispvec(z,false,std::cerr);
                throw std::logic_error("le ff in fx");
            }catch(...){
                std::cerr<<"ff in fx other error"<<std::endl;
                throw;
            }
        };
    };
    
    if(off and 0==neo){//offset only parameter
        Dvector oi(nd);
        try{
            func(P,iparam,oi);
        }catch(...){
            std::cerr<<"offset error"<<std::endl;
            throw;
        }
        newC(oi,P[0],data,sig);
        eparam=P;
        vals=oi;
        chisq=chisquared(vals,data,sig);
        if(disp){
            std::cerr<<"\ndone offset only:"<<chisq<<std::endl;
            dispvec(eparam,false,std::cerr);
            dispvec(iparam,false,std::cerr);
            dispvec(vals,false,std::cerr);
        }
        return;
    }

    Dvector p0=P,//P0
    w=iparam,//internal parameters
    oi=Dvector(nd,0);//outputs
    Double f0;
    try{
        f0=ff(p0,w);//f0
    }catch(std::invalid_argument&e){
        std::cerr<<"invalid argument ff in f0: "<<e.what()<<'\n';
        dispvec(p0,false,std::cerr);
        throw std::logic_error("ff in f0");
    }catch(std::logic_error&e){
        std::cerr<<"logic error ff in f0: "<<e.what()<<'\n';
        dispvec(p0,false,std::cerr);
        throw std::logic_error("ff in f0");
    }catch(...){
        std::cerr<<"ff in f0 other error"<<std::endl;
        throw;
    }

    std::size_t steps=0;

    //std::cerr<<"chi loop\n";
    do{
        Dvector pi=p0;//current position
        //pt=p0;//previous position
        Double fi=f0,//current value
        ft=f0,//previous value
        bigf=0,//largest fall
        cc=pi[neo];//previous offset
        std::size_t bigi=0;//direction of biggest fall
        std::array<Double,6> l;//lmin output

        if(disp){
            try{
                func(pi,w,oi);
                fi=chisquared(oi,data,sig);
            }catch(...){
                std::cerr<<"minimize 1 error"<<std::endl;
                throw;
            }
            std::cerr<<"\nstep "<<steps<<':'<<fi<<'\n';
            dispvec(pi,false,std::cerr);
            dispvec(w,false,std::cerr);
            dispvec(oi,false,std::cerr);
            if(!(steps%10)) dispvec(data,false,std::cerr);
        }
        steps++;

        //minimize in each direction
        for(std::size_t i=0;i<neo;i++){
            std::size_t ii=i;
            //pt=pi;
            ft=fi;
            Dvector wt=w,ww=w;
            //std::cerr<<"lmin 1"<<std::endl;
            for(;ii<neo;ii++){
                l=lmin(fx(pi,xmat[ii],w),ierror,rel);
                fi=l[4];
                
                if(fi<ft+2*oerror) break;
                
                std::cerr<<"failed to minimize "<<i<<' '<<ii<<std::endl;
                Double ftt;
                try{
                    ftt=ff(pi,ww);
                }catch(...){
                    std::cerr<<"minimize 2 error"<<std::endl;
                    throw;
                }
                if (ftt>fi or ft>ftt-2*oerror) break;
                ww=wt;
            }
            if(i!=ii){
                if(i<neo) swap(xmat[i],xmat[ii]);
                else swap(xmat[i],xmat[neo-1]);
                w=ww;
            }
#pragma omp parallel for
            for(std::size_t j=0;j<neo;j++){
                pi[j]+=l[1]*xmat[i][j];
                Double temp=xmat[i][j]*l[1];
                if(std::abs(temp)<(pi[j]*EPS) and 0!=xmat[i][j])
                    xmat[i][j]=std::copysign(std::max(std::abs(xmat[i][j]),
                                                      std::abs(pi[j])/100),
                                             temp);
                else xmat[i][j]=temp;
            }
            if(off){
                try{
                    func(pi,w,oi);
                }catch(...){
                    std::cerr<<"minimize 3 error"<<std::endl;
                    throw;
                }
                cc=pi[neo];
                newC(oi,pi[neo],data,sig);
                fi=chisquared(oi,data,sig);
            }
            if(disp){
                Dvector xx=xmat[i];
                if(!off){
                    try{
                        func(pi,w,oi);
                        fi=chisquared(oi,data,sig);
                    }catch(...){
                        std::cerr<<"minimize 4 error"<<std::endl;
                        throw;
                    }
                }else{
                    xx.resize(ne,pi[neo]-cc);
                }
                std::cerr<<"\ndir "<<i<<':'<<fi<<std::endl;
                dispvec(xx,false,std::cerr);
                dispvec(pi,false,std::cerr);
                dispvec(oi,false,std::cerr);
            }

            //check largest step
            if((ft-fi)>bigf){
                bigf=ft-fi;
                bigi=i;
            }
        }

        //terminiation condition
        if((f0-fi)<oerror*(rel?std::max(f0,EPS):1)*(df?df:1) or fi<=df
           or 1==pi.size() or pi==p0){
            eparam=pi;
            try{
                func(pi,w,vals);
            }catch(...){
                std::cerr<<"minimize 5 error"<<std::endl;
                throw;
            }
            iparam=w;
            chisq=chisquared(vals,data,sig);
            if(disp){
                std::cerr<<"\ndone "<<steps<<':'<<chisq<<std::endl;
                dispvec(eparam,false,std::cerr);
                dispvec(iparam,false,std::cerr);
                dispvec(vals,false,std::cerr);
            }
            break;
        }

        Dvector pe(ne),//extrapolation point
        xstep(neo);//step size
#pragma omp parallel for
        for(std::size_t i=0;i<neo;i++){
            pe[i]=2*pi[i]-p0[i];
            xstep[i]=pi[i]-p0[i];
        }
        Double fe;
        try{
            fe=ff(pe,w);
        }catch(...){
            std::cerr<<"minimize 6 error"<<std::endl;
            throw;
        }

        if(fe<f0){
            if(2*(f0-2*fi+fe)*sqr(f0-fi-bigf)<bigf*sqr(f0-fe)){
                //std::cerr<<"lmin 2"<<std::endl;
                Dvector wt=w;
                l=lmin(fx(pi,xstep,w),ierror,rel);

                if(l[4]>fi){
                    try{
                        ft=ff(pi,wt);
                    }catch(...){
                        std::cerr<<"minimize 7 error"<<std::endl;
                        throw;
                    }
                    if(ft>l[4]+2*oerror) continue;
                }else fi=l[4];
#pragma omp parallel for
                for(std::size_t i=0;i<ne;i++) pi[i]+=l[1]*xstep[i];
                xmat[bigi]=xmat[neo-1];
                xmat[neo-1]=xstep;
                if(off){
                    try{
                        func(pi,w,oi);
                    }catch(...){
                        std::cerr<<"minimize 8 error"<<std::endl;
                        throw;
                    }
                    cc=pi[neo];
                    newC(oi,pi[neo],data,sig);
                    fi=chisquared(oi,data,sig);
                }
                p0=pi;
                f0=fi;
                if(disp){
                    Dvector xx=xmat[neo-1];
                    if(!off){
                        try{
                            func(p0,w,oi);
                        }catch(...){
                            std::cerr<<"minimize 9 error"<<std::endl;
                            throw;
                        }
                    }else{
                        xx.resize(ne,p0[neo]-cc);
                    }
                    f0=chisquared(oi,data,sig);
                    std::cerr<<"\next min "<<f0<<std::endl;
                    dispvec(xx,false,std::cerr);
                    dispvec(p0,false,std::cerr);
                    dispvec(oi,false,std::cerr);
                }
            }else if(fe<fi){
                p0=pe;
                f0=fe;
                if(disp){
                    try{
                        func(p0,w,oi);
                    }catch(...){
                        std::cerr<<"minimize 10 error"<<std::endl;
                        throw;
                    }
                    f0=chisquared(oi,data,sig);
                    std::cerr<<"\next "<<f0<<std::endl;
                    dispvec(p0,false,std::cerr);
                    dispvec(oi,false,std::cerr);
                }
            }else{
                p0=pi;
                f0=fi;
            }
        }else{
            p0=pi;
            f0=fi;
        }
        
        //check if xmat singular
        bool b=false;
        for(std::size_t i=1;i<neo;i++) if(xmat[i]==xmat[0]){
            b=true;
            break;
        }
        if(b){
            if(disp){
                std::cerr<<"xmat singular"<<std::endl;
                dispmat(xmat,false,std::cerr);
            }
            xmat=Dmatrix(neo,Dvector(neo,0));
#pragma omp parallel for
            for(std::size_t i=0;i<neo;i++) xmat[i][i]=p0[i]*.1;
        }
    }while(true);
}

