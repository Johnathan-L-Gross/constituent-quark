//
//  energy.cpp
//
//
//  Created by Johnathan Gross
//
//

#include "energy.h"
#include "min/min.h"
#ifdef TESTING_STATES_JLG
#include "mout.h"
#include <iostream>
#endif
#include <armadillo>
#include <stdexcept>

inline auto index_comp(const Dvector&v){
    /*
     returns lambda function that compares indexes using values
     of vector
     */
    return [&v](const std::size_t i,const std::size_t j){return v[i]<v[j];};
};

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

/*------------------------------------------------------------------------------
 state sort
------------------------------------------------------------------------------*/
void state_sort(tvector&ilist,const Dvector&en,const bool flag){
    if(flag) std::stable_sort(ilist.begin(),ilist.end(),index_comp(en));
    else std::sort(ilist.begin(),ilist.end(),index_comp(en));
}

/*------------------------------------------------------------------------------
 declaring eigenvalue functions
------------------------------------------------------------------------------*/
#ifdef _poly_comp_JLG
void eigen(const arma::sp_cx_mat&,//ham
           arma::cx_vec&,//eigv
           arma::cx_mat&,//eigm
           const arma::uword,//k
           const bool=0);//flag
/*
 Calculates eigenvalues and eigenvectors of ham, most general
 Effect: writes eigenvalues to eigv and eigenvectors to eigm
 */
void eigen(const arma::sp_cx_mat&,//ham
           arma::vec&,//eigv
           arma::cx_mat&,//iegm
           const arma::uword,//k
           const bool=0);//flag
//overload for complex Hermetian matricies
void eigen(const arma::sp_mat&,//ham
           arma::cx_vec&,//eigv
           arma::cx_mat&,//iegm
           const arma::uword,//k
           const bool=0);//flag
//overload for real matricies
#endif
void eigen(const arma::sp_mat&,//ham
           arma::vec&,//eigv
           arma::mat&,//iegm
           const arma::uword);//k
//overload for real symmetric matricies


#ifdef _poly_comp_JLG
/*------------------------------------------------------------------------------
complex
------------------------------------------------------------------------------*/
void eigen(const arma::sp_cx_mat&ham,
           arma::cx_vec&eigv,
           arma::cx_mat&eigm,
           const arma::uword k,
           const bool flag){
    eigv.zeros(ham.n_rows);
    eigm.zeros(ham.n_rows,ham.n_cols);
#ifdef TESTING_STATES_JLG
    //std::cerr<<"complex\n";
    if(!ham.is_square()) throw std::invalid_argument("Matrix not square");
#endif //TESTING_STATES_JLG
    if(!flag){
        //test if symmetric or real
        bool sym=ham.is_hermitian(),r=true;
        for(auto it=ham.begin();it!=ham.end();++it){
            if((*it).imag()){
                r=false;
                break;
            }
        }

        if(sym and r){//real and symmetric
            arma::vec eigv2(real(eigv));
            arma::mat eigm2(real(eigm));
            arma::sp_mat ham2(real(ham));
            eigen(ham2,eigv2,eigm2,k);
            eigv=arma::conv_to<arma::cx_vec>::from(eigv2);
            eigm=arma::conv_to<arma::cx_mat>::from(eigm2);
            return;
        }
    }
    eigs_gen(eigv,eigm,ham,k,"sr"); //not real nor symmetric
    //#warning "do not use complex, asymmetric eigen"
}

/*------------------------------------------------------------------------------
 Hermitian
 -----------------------------------------------------------------------------*/
void eigen(const arma::sp_cx_mat&ham,
           arma::vec&eigv,
           arma::cx_mat&eigm,
           arma::uword k,
           const unsigned char flag){
    eigv.zeros(ham.n_rows);
    eigm.zeros(ham.n_rows,ham.n_cols);
#ifdef TESTING_STATES_JLG
    //std::cerr<<"Hermitian\n";
    if(!ham.is_hermitian()) throw std::invalid_argument("Matrix not symmetric");
    if(!ham.is_square()) throw std::invalid_argument("Matrix not square");
#endif
    if(!flag){
        //test if real (and therefore symmetric)
        bool r=true;
        for(auto it=ham.begin();it!=ham.end();++it){
            if((*it).imag()){
                r=false;
                break;
            }
        }

        if(r){//real
            arma::sp_mat ham2(real(ham));
            arma::mat eigm2(real(eigm));
            eigen(ham2,eigv,eigm2,k);
            eigm=arma::conv_to<arma::cx_mat>::from(eigm2);
            return;
        }
    }
    arma::cx_vec eigv2=arma::conv_to<arma::cx_vec>::from(eigv);
    eigs_gen(eigv2,eigm,ham,k,"sr");
    eigv=real(eigv2);
}

/*------------------------------------------------------------------------------
 real
 -----------------------------------------------------------------------------*/
void eigen(const arma::sp_mat&ham,
           arma::cx_vec&eigv,
           arma::cx_vec&eigm,
           arma::uword k,
           const unsigned char flag){
    eigv.zeros(ham.n_rows);
    eigm.zeros(ham.n_rows,ham.n_cols);
#ifdef TESTING_STATES_JLG
    //std::cerr<<"real";
    if(!ham.is_square()) throw std::invalid_argument("Matrix not square");
#endif
    if(!flag and ham.is_symmetric()){
        arma::vec eigv2(real(eigv));
        arma::mat eigm2(real(eigm));
        eigen(ham,eigv2,eigm2,k);
        eigv=arma::conv_to<arma::cx_vec>::from(eigv2);
        return;
    }
    arma::sp_cx_mat ham2=arma::conv_to<arma::sp_cx_mat>::from(ham);
    eigs_gen(eigv,eigm,ham2,k,"sr");
}

#endif //_poly_comp_JLG

/*------------------------------------------------------------------------------
 real symmetric
 -----------------------------------------------------------------------------*/
void eigen(const arma::sp_mat&ham,
           arma::vec&eigv,
           arma::mat&eigm,
           arma::uword k){
    eigv.zeros(ham.n_rows);
    eigm.zeros(ham.n_rows,ham.n_cols);
#ifdef TESTING_STATES_JLG
    //std::cerr<<"real symmetric\n";
    if(!ham.is_square()) throw std::invalid_argument("Matrix not square");
    if(!ham.is_symmetric())throw std::invalid_argument("Matrix not symmetric");
#endif
    eigs_sym(eigv,eigm,ham,k,"sa");
    if(!eigv.is_sorted()) throw std::logic_error("unsorted");
}

/*------------------------------------------------------------------------------
 energy
 -----------------------------------------------------------------------------*/
Double energy(const sparse&ham,
              const std::size_t N,
              const std::size_t k){
    //Nth energy level
    //const size_t sham=ham.size();
    if(N>=ham.m()) throw std::out_of_range("not enough energy levels");
    arma::uword kk=std::min(k,ham.m()-1);
    if(N>kk) throw std::invalid_argument("N>k");
    for(auto it=ham.cbegin();it!=ham.cend();it++)if(std::isnan(it->second)){
        throw std::invalid_argument("entry is NaN");
    }

    arma::vec eval(ham.m(),arma::fill::zeros);
    arma::mat evec(ham.m(),ham.n(),arma::fill::zeros);
    auto h=sparse2arma(ham);
#ifdef TESTING_STATES_JLG
    try{
        eigen(h,eval,evec,kk);
    }catch(std::invalid_argument&e) {
        std::cerr<<"energy invalid argument error\n";
        //dispmat(ham,false,std::cerr);
        //dispmat(ham-ham.transpose(),false,std::cerr);
        throw e;
    }catch(...){
        std::cerr<<"energy other error"<<std::endl;
        throw;
    }
#else
    eigen(h,eval,evec,kk);
#endif

    return eval[N];
}

Double energy(const sparse&ham,
              Dvector&eval,
              const std::size_t N,
              const std::size_t k){
    //Nth energy level, stores eigenvalues
    if(N>=ham.m()) throw std::out_of_range("not enough energy levels");
    arma::uword kk=std::min(k,ham.m()-1);
    if(N>kk) throw std::invalid_argument("N>k");
    for(auto it=ham.cbegin();it!=ham.cend();it++)if(it->second!=it->second){
        throw std::invalid_argument("entry is NaN");
    }

    arma::vec eval2(ham.m(),arma::fill::zeros);
    arma::mat evec(ham.m(),ham.n(),arma::fill::zeros);
    auto h=sparse2arma(ham);
    eigen(h,eval2,evec,kk);

    eval=arma::conv_to<Dvector>::from(eval2);
    return eval[N];
}

void energy(const sparse&ham,
            sparse&evec,
            Dvector&eval,
            const std::size_t k){
    for(auto it=ham.cbegin();it!=ham.cend();it++)if(it->second!=it->second){
        throw std::invalid_argument("entry is NaN");
    }
    //stores eigenvalues and eigenvectors
    arma::uword kk=std::min(k,ham.m()-1);
    arma::vec eval2(ham.m(),arma::fill::zeros);
    arma::mat evec2(ham.m(),ham.n(),arma::fill::zeros);
    auto h=sparse2arma(ham);
    eigen(h,eval2,evec2,kk);
    eval=arma::conv_to<Dvector>::from(eval2);
    evec=arma2sparse(arma::sp_mat(evec2));
    
    //largest element of eigenvector is positive
#pragma omp parallel for
    for(size_t i=0;i<kk;i++){
        Double m=std::abs(evec(0,i));
        std::size_t j=0;
        for(std::size_t jj=1;jj<evec.n();jj++){
            if(evec.count(jj,i)){
                const Double mm=std::abs(evec(jj,i));
                if(mm>m){
                    m=mm;
                    j=jj;
                }
            }
        }
        if(evec(j,i)<0)for(std::size_t jj=0;jj<evec.n();jj++){
            if(evec.count(jj,i)) evec.set_entry(jj,i,-evec(jj,i));
        }
    }
}

/*------------------------------------------------------------------------------
 minimization
 -----------------------------------------------------------------------------*/
std::array<Double,4> minen(std::function<void(sparse&,Double)>&Hfunc,
                           Double par,
                           Double step,
                           const std::size_t N,
                           const std::size_t k,
                           Double error,
                           const bool rel,
                           const bool pos,
                           const bool disp){
    //std::cerr<<"minen"<<std::endl;
    auto f=[&Hfunc,N,k](const Double p){
        //std::cerr<<"minen f"<<std::endl;
        sparse h;
        Hfunc(h,p);
        //std::cerr<<p<<std::endl;
        //dispmat(h,false,std::cerr);
        try{
            return energy(h,N,k);
            //std::cerr<<"minen f end"<<std::endl;
            //return en;
        }catch(std::exception&e){
            std::cerr<<e.what()<<' '<<p<<std::endl;
            throw e;
        }catch(...){
            std::cerr<<"minen error"<<std::endl;
            throw;
        }
    };
    
    error=std::max(std::numeric_limits<Double>::min(),error);
    //std::cerr<<"minen"<<std::endl;
    std::array<Double,6> s=Brentmin(f,par,par+step,error,rel,pos,disp);
    //std::cerr<<"b"<<std::endl;
    Double e=std::max(s[3],s[5])-s[4],
    d=std::max(s[1]-s[0],s[2]-s[1]);
#ifdef TESTING_STATES_JLG
    if(error>0 and e>error*(rel?std::abs(s[4]):1))
        std::cerr<<"Minimization did not produce low enough error"<<std::endl;
#endif
    //std::cerr<<"minen end"<<std::endl;
    return {s[1],d,s[4],e};
}
