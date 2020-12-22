//
//  int.cpp
//  
//
//  Created by Johnathan Gross
//
//

#include "int.h"
#ifdef TESTING_STATES_JLG
#include <stdexcept>
#include <iostream>
#endif

inline const long double M_s2PI=
0.39894228040143267793994605993438186847585863116493465766592582967065792589930183850125233390730693643030255886263L;

int1::int1(const bool polytype,
           const unsigned int N,
           const Dvector&XX,
           const Dvector&XW):
N_(N),polytype_(polytype),n_(XX.size()),XX_(XX),XW_(XW),
lag_(XX.size(),Dmatrix(N+1,Dvector(N/2+1,0))),
rad_(XX.size(),Dmatrix(N/2+1,Dvector(N+1,0))),dx_(XX.size(),0){
    setup();
}

#ifdef D3_JLG
/*int3::int3(const bool polytype,const unsigned int N,const Dvector& XX,
           const Dvector&XW,const Dvector&LX,const Dvector&LW):
int1(polytype,N,XX,XW),l_(LX.size()),LX_(LX),LW_(LW),Y_(l_,Dmatrix(N_+1)),
Y2_(l_,Dmatrix(N_+1)){
    setup();
}*/

int3::int3(const int1& A,
           const unsigned int L,
           const Dvector& LX,
           const Dvector& LW):
int1(A),L_(L),l_(LX.size()),LX_(LX),LW_(LW),
Y_(LX.size(),Dmatrix(L+1)),Y2_(LX.size(),Dmatrix(L+1)){
    setup();
}
#endif //D3_JLG

#ifdef D6_JLG
/*int6::int6(const bool polytype,const unsigned int N,const Dvector& XX,
           const Dvector& XW,const Dvector& LX,const Dvector& LW,
           const Dvector& PX,const Dvector& PW):
int3(polytype,N,XX,XW,LX,LW),m_(2*PX.size()),m2_(PX.size()),PX_(PX),PW_(PW),
dphi_(m_,0),phi_(m_,0),phase_(m_,CDvector(N_+1,1.0)),phase2_(m_,CDvector(N_+1,1.0)){
    setup();
}

int6::int6(const int1& A,const Dvector& LX,const Dvector& LW,const Dvector& PX,
           const Dvector& PW):
int3(A,LX,LW),m_(2*PX.size()),m2_(PX.size()),PX_(PX),PW_(PW),dphi_(m_,0),
phi_(m_,0),phase_(m_,CDvector(N_+1,1.0)),phase2_(m_,CDvector(N_+1,1.0)){
    setup();
}*/

int6::int6(const int3&A,
           const Dvector& PX,
           const Dvector& PW):
int3(A),m_(2*PX.size()),m2_(PX.size()),PX_(PX),PW_(PW),
dphi_(2*PX.size(),0),phi_(2*PX.size(),0),
phase_(2*PX.size(),CDvector(L_+1,1.0)),phase2_(2*PX.size(),CDvector(L_+1,1.0)){
    setup();
}
#endif

void int1::setup(){
    //std::cerr<<"setup1\n";
#ifdef TESTING_STATES_JLG
    if(XX_.size()!=XW_.size())
        throw std::length_error("1: abscissa and weights not the same length");
#endif

    harmoniclx(N_/2,N_,harm_);

#pragma omp parallel for
    for(std::size_t x=0;x<n_;x++){
        Double xx,xx2;
        if(polytype_){
            xx=XX_[x];//GHQx.at(GHQn);
            xx2=sqr(xx);
        }else{
            xx2=XX_[x];//GLQx.at(GLQn);
            xx=std::sqrt(xx2);
        }
        for(unsigned int l=0;l<=N_;++l){
            laguerrex(N_/2,Double(l)+0.5,xx2,lag_[x][l]);
            //ideally only need to go to (N-l)/2
#ifdef TESTING_STATES_JLG
            if(lag_[x][l].size()!=N_/2+1)
                throw std::logic_error("1: lag_ wrong size");
#endif
            for(std::size_t i=0;i<lag_[x][l].size();++i){
                const int s=sign(lag_[x][l][i]);
                if(0!=s){
                    Double r=std::log(std::abs(lag_[x][l][i]))+harm_[i][l];
                    if(polytype_) r+=std::log(xx)*l;
                    else r+=std::log(xx2)*l/2.;
                    rad_[x][i][l]=s*std::exp(r);
                    //rad2_[x][l][i]=(l%2?-rad_[x][l][i]:rad_[x][l][i]);
                    //cerr<<i<<' '<<l<<' '<<rad[l][i]<<'\n';
                }//else rad[x][l][i]=0;
            }//i
        }//l
        if(polytype_) dx_[x]=XW_[x]*xx2;//GHQw.at(GHQn)[x]*xx2;
        else dx_[x]=XW_[x]*xx/2.;//GLqw.at(GLQn)[x]*xx/2;
    }//x
#pragma omp parallel for collapse(2)
    for(std::size_t i=0;i<harm_.size();++i)
        for(std::size_t j=0;j<harm_[0].size();++j)
            harm_[i][j]=std::exp(harm_[i][j]);
}

#ifdef D3_JLG
void int3::setup(){
    //std::cerr<<"steup3\n";
#ifdef TESTING_STATES_JLG
    if(LX_.size()!=LW_.size())
        throw std::length_error("3: abscissa and weights not the same length");
#endif

    //int1::setup();//x
    //cos(theta)
    for(std::size_t x=0;x<l_;x++){
        nalegendrex(L_,LX_[x],Y_[x]);//GPQx.at(GPQl)
        nalegendrex(L_,-LX_[x],Y2_[x]);//GPQx.at(GPQl)
        //nalegendrex(N_,0,LX_[x],ang_[x]);
        //nalegendrex(N_,0,-LX_[x],ang2_[x]);
    }
}
#endif

#ifdef D6_JLG
void int6::setup(){
    //std::cerr<<"setup6\n";
#ifdef TESTING_STATES_JLG
    if(PX_.size()!=PW_.size())
        throw std::length_error("6: abscissa and weights not the same length");
    if(L_>=m2_) std::cerr<<"Warning: 6D integration may be inacurate for large m\n"<<std::endl;
#endif
    //int3::setup();//x and cos(theta)
    //phi
#pragma omp parallel for
    for(std::size_t x=0;x<m2_;x++){
        phi_[m2_-x-1]=phix(-PX_[x]);
        phi_[m2_+x]=phix(PX_[x]);
        dphi_[m2_-x-1]=M_PI*PW_[x];
        dphi_[m2_+x]=M_PI*PW_[x];
    }
#pragma omp parallel for collapse(2)
    for(std::size_t x=0;x<m_;x++)for(std::size_t mm=0;mm<=L_;++mm){
        phase_[x][mm]=std::polar(Double(M_s2PI),phi_[x]*mm);
        phase2_[x][mm]=std::polar(Double(M_s2PI),-phi_[x]*mm);
    }
}

#endif
