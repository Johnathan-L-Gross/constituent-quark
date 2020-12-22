//
//  wf.cpp
//  
//
//  Created by Johnathan Gross
//
//

#include "wf.h"
#include "states.h"
#include <iostream>
#ifdef TESTING_STATES_JLG
#include <stdexcept>
#endif

constexpr std::size_t numberlength(int i){
    //calculates the displayed length of an integer number
    std::size_t n=1,ten=10;
    unsigned int ii=0;
    if(i<0){
        n++;
        ii=-i;
    }
    else ii=i;

    while(ii>=ten){
        ten*=10;
        n++;
    }
    return n;
}
constexpr std::size_t numberlength(std::size_t i){
    //overload for std::size_t
    std::size_t n=1,ten=10;
    while(i>=ten){
        ten*=10;
        n++;
    }
    return n;
}

//*****************************************************************************
//wf1
wf1::wf1(const imatrix&s,
         const marray&mass_):
states(s),
Pf(s[0][0]),J(s[0][1]),P(s[0][2]),M(s[0][3]),
Nmin(s[0][4]),Nmax(s[s.size()-1][4]),mass(mass_),sym(true),size(s.size()){
    //std::cerr<<"wf1 s\n";

    bool x=(s[0][0]+s[0][7]+s[0][11])%2;
    for(auto y:s) if(x!=(y[0]+y[7]+y[11])%2){
        sym=false;
        break;
    }
#ifdef TESTING_STATES_JLG
    if(mass[0]!=mass[1] and sym){
        std::cerr<<"Different mass1 and mass2 with definite symmetry "
        <<mass[0]<<' '<<mass[1]<<' '<<mass[2]<<std::endl;
    }
#endif
    
    transformss(states,mass,1,mosh1);
    transformss(states,mass,2,mosh2);
}

wf1::wf1(const int Pf_,
         const int J_,
         const int P_,
         const int M_,
         const int Nmin_,
         const int Nmax_,
         const marray&mass_,
         const bool sym_):
Pf(Pf_),J(J_),P(P_),M(M_),Nmin(Nmin_),Nmax(Nmax_),mass(mass_),sym(sym_){
    //std::cerr<<"wf1 all"<<std::endl;

    #ifdef TESTING_STATES_JLG
        if(mass[0]!=mass[1] and sym){
            std::cerr<<"Different mass1 and mass2 with definite symmetry "
            <<mass[0]<<' '<<mass[1]<<' '<<mass[2]<<std::endl;
        }
    #endif
    
    statematrix(Pf,J,P,M,Nmin,Nmax,states,sym);
    size=states.size();
    /*std::cerr<<"statematrix complete"<<std::endl;
    for(auto s:*states){
        for(auto i:s) std::cerr<<i<<' ';
        std::cerr<<'\n';
    }*/

    transformss(states,mass,1,mosh1);
    transformss(states,mass,2,mosh2);
    //std::cerr<<"transformations complete"<<std::endl;
    //std::cerr<<"wf1 all complete"<<std::endl;
}

/*
void swap(wf1&first,wf1&second){
    //std::cerr<<"wf1 swap\n";
    using std::swap;
    swap(first.states,second.states);
    swap(first.Pf,second.Pf);
    swap(first.J,second.J);
    swap(first.P,second.P);
    swap(first.M,second.M);
    swap(first.Nmin,second.Nmin);
    swap(first.Nmax,second.Nmax);
    swap(first.mass,second.mass);
    swap(first.sym,second.sym);
    swap(first.size,second.size);
    swap(first.mosh1,second.mosh1);
    swap(first.mosh2,second.mosh2);
    swap(first.dstates,second.dstates);
    swap(first.dmosh1,second.dmosh1);
    swap(first.dmosh2,second.dmosh2);
}
 */

void wf1::display(const std::size_t s,std::ostream&out)const{
    const size_t ss=std::min(s,size);
    const size_t numlen=numberlength(ss);
    out<<"#:Pf,J,Pj,M,N,L,np,lp,nl,ll,S,Ps\n";
    for(size_t i=0;i<ss;i++){
        out<<i;
        size_t ilen=numberlength(i);
        for(size_t j=ilen;j<numlen;j++) out<<' ';
        out<<':';
        for(size_t k=0;k<states[i].size();k++)
        if(0==k or 2==k) out<<(states[i][k]%2?'-':'+')<<' ';
        else if(11==k) out<<(states[i][10]==3?
                             states[i][k]%2?'A':'S':
                             states[i][k]%2?'r':'l');
        else out<<states[i][k]<<' ';
        out<<'\n';
    }
}

sparse wf1::P12()const{
    Dvector x(states.size());
#pragma omp parallel for
    for(std::size_t i=0;i<x.size();i++)
        x[i]=((states[i][0]+states[i][7]+states[i][11])%2?1:-1);
    //(-1)^(1+Pf+lp+Ps)
    return sparse(x);
}
sparse wf1::P12ss()const{
    Dvector x(states.size());
#pragma omp parallel for
    for(std::size_t i=0;i<x.size();i++)
        x[i]=((states[i][7]+states[i][11])%2?-1:1);
        //(-1)^(lp+Ps)
    return sparse(x);
}
sparse wf1::P13()const{
    return (Pf?-1:1)*mosh2*P12();
    //mosh2 does not account for color (-1 to account for that)
    //or flavor (assumes either fully symmetric or fully antisymmetric)
}
sparse wf1::P13ss()const{
    return mosh2*P12ss();
}
sparse wf1::P23()const{
    return (Pf?-1:1)*mosh1*P12();
    //mosh1 does not account for color (-1 to account for that)
    //or flavor (assumes either fully symmetric or fully antisymmetric)
}
sparse wf1::P23ss()const{
    return mosh1*P12ss();
}

//*****************************************************************************
//wf2
wf2::wf2(const wf1&wa,
         const wf1&wb):
wfa(wa),wfb(wb){
    //std::cerr<<"wf2 ab"<<std::endl;
#ifdef TESTING_STATES_JLG
    if(wfa.mass!=wfb.mass) throw std::invalid_argument("not same mass");
    if(wfa.size>wfb.size){
        std::cerr<<"A bigger than B"<<std::endl;
        std::swap(wfa,wfb);
    }
#endif

    {//scope j
        trans=tvector(0);
        size_t j=0;
        for(size_t i=0;i<wfb.size;i++) if(wfa.states[j]==wfb.states[i]){
            trans.push_back(i);
            j++;
            if(wfa.size==j){
                intern=true;
                break;
            }
        }
    }
    //for(auto x:(*trans)) std::cerr<<x<<' ';
    //std::cerr<<std::endl;
    
    if(intern){
        const size_t s1=wfa.size,s2=wfb.size;
        mosh1=sparse(s1,s2);
        mosh2=sparse(s1,s2);
        for(size_t m=0;m<s1;m++)for(size_t n=0;n<s2;n++){
            mosh1.set_entry(m,n,wfb.mosh1[{trans[m],n}]);
            mosh2.set_entry(m,n,wfb.mosh2[{trans[m],n}]);
        }
    }else trans=tvector(0);
}

/*
void swap(wf2&first,wf2&second){
    //std::cerr<<"wf2 swap\n";
    using std::swap;
    swap(first.wfa,second.wfa);
    swap(first.wfb,second.wfb);
    swap(first.intern,second.intern);
    swap(first.trans,second.trans);
    swap(first.mosh1,second.mosh1);
    swap(first.mosh2,second.mosh2);
    swap(first.dwfa,second.dwfa);
    swap(first.dwfb,second.dwfb);
    swap(first.dtrans,second.dtrans);
    swap(first.dmosh1,second.dmosh1);
    swap(first.dmosh2,second.dmosh2);
}
*/
