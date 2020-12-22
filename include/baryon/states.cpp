//
//  states.cpp
//
//
//  Created by Johnathan Gross
//
//

#include "states.h"
#include "mosh/mosh.h"

#include <iostream>
//#include <utility>
#include <algorithm>

#ifdef TESTING_STATES_JLG
#include <stdexcept>
#endif

inline bool istriangleL(const int i,const int j,const int k){
    return i+j>=k and std::abs(i-j)<=k;
}

inline bool compstate(const ivector&a,const ivector&b){
    //true if "a<b"
    if(mod(a[0],2)>mod(b[0],2)) return true;//Pf Lambda before Sigma
    else if(mod(a[0],2)<mod(b[0],2)) return false;
    if(a[1]<b[1]) return true;//J
    else if(a[1]>b[1]) return false;
    if(mod(a[2],2)<mod(b[2],2)) return true;//P + before -
    else if(mod(a[2],2)>mod(b[2],2)) return false;
    if(a[3]<b[3]) return true;//M
    else if(a[3]>b[3]) return false;
    if(a[4]<b[4]) return true;//N
    else if(a[4]>b[4]) return false;
    if(a[5]<b[5]) return true;//L
    else if(a[5]>b[5]) return false;
    if(a[7]+a[9]<b[7]+b[9]) return true;//lp+ll
    else if(a[7]+a[9]>b[7]+b[9]) return false;
    if(a[6]<b[6]) return true;//np
    else if(a[6]>b[6]) return false;
    if(a[7]<b[7]) return true;//lp
    else if(a[7]>b[7]) return false;
    if(a[10]<b[10]) return true;//S
    else if(a[10]>b[10]) return false;
    if(mod(a[0]+a[7]+a[11],2)<mod(b[0]+b[7]+b[11],2)) return true;//f-s-l sym
    else if(mod(a[0]+a[7]+a[11],2)>mod(b[0]+b[7]+b[11],2)) return false;
    //otherwise states are equal
    return false;
}

//******************************************************************************
//******************************************************************************
void statematrix(const int Pf,
                 const int J,
                 const int P,
                 const int M,
                 const int Nmin,
                 const int Nmax,
                 imatrix&slist,
                 const bool sym){
    {
        int error=0;
        if(J<0) error=1;
        else if(!(J%2)) error=2;
        else if(0!=P and 1!=P) error=3;
        else if(!(M%2)) error=4;
        else if(std::abs(M)>J) error=5;
        else if(0>Nmax) error=6;
        else if(Nmin>Nmax) error=7;
#ifdef TESTING_STATES_JLG
        switch (error) {
            case 1:
                throw std::invalid_argument("J<0");
            case 2:
                throw std::invalid_argument("J must be odd");
            case 3:
                throw std::invalid_argument("P must be 0 or 1");
            case 4:
                throw std::invalid_argument("M should be odd");
            case 5:
                throw std::invalid_argument("abs(M)>J");
            case 6:
                throw std::invalid_argument("Nmax<0");
            case 7:
                throw std::invalid_argument("Nmin>Nmax");
        }
#else
        if(error){
            slist.resize(0);
            std::cerr<<"error code "<<error<<std::endl;
            return;
        }
#endif
    }
    slist.resize(0);
    
    //parity preserving minimum value
    auto PP=[P](int x){return (x+P)%2?x+1:x;};

    const int Lmin=std::max(P,(J-SMAX_B)/2);
    const int Nm=std::max(PP(Lmin),PP(Nmin));
    //Nmin already used
    const int Lmax=(J+SMAX_B)/2;

//#ifdef TESTING_STATES_JLG
    //if(Nm>Nmax) throw std::logic_error("No allowed N");
//#endif

    //int N,L,np,lp,nl,ll,S,Ps,lpl,lmx,n;
    //Pf,J,P,M,N,L,np,lp,nl,ll,S,Ps
    //lpl,lmx,n not used to label state

    //N (+2)
//#pragma omp parallel for if(sort) schedule(dynamic)
    for(int N=Nm;N<=Nmax;N+=2){
        const int LM=std::min(Lmax,N);
        //L (+1)
        for(int L=Lmin;L<=LM;L++){
            const bool is1=istriangleL(J,1,2*L),
            is3=istriangleL(J,3,2*L);
            if(!is1 and !is3) continue;
            //lp+ll (+2)
            for(int lpl=PP(L);lpl<=N;lpl+=2){
                const int lmx=(L+lpl)/2;
                const int lmi=lpl-lmx;
#ifdef TESTING_STATES_JLG
                if(lmi<0) throw std::logic_error("lmi<0");
#endif
                const int n=(N-lpl)/2;
#ifdef TESTING_STATES_JLG
                if(2*n!=N-lpl) throw std::logic_error("2n!=N-(lp+ll)");
#endif
                //np (+1)
                for(int np=0;np<=n;np++){
                    const int nl=n-np;
                    //lp (+1)
                    for(int lp=lmi;lp<=lmx;lp++){
                        const int ll=lpl-lp;
#ifdef TESTING_STATES_JLG
                        if(!istriangleL(lp,ll,L))
                            throw std::logic_error("lp,ll,L not triangle");
#endif
                        const int Ps=mod(Pf+lp,2);
                        //#pragma omp critical
                        //{
                        if(is1){
                            slist.push_back({Pf,J,P,M,N,L,np,lp,nl,ll,1,Ps});
                            if(!sym) slist.push_back({Pf,J,P,M,N,L,np,lp,nl,ll,
                                1,mod(Ps+1,2)});
                        }
                        if(is3 and (!Ps or !sym))
                            slist.push_back({Pf,J,P,M,N,L,np,lp,nl,ll,3,0});
                        //}
                    }//lp
                }//np
            }//lp+ll
        }//L
    }//N
#ifdef TESTING_STATES_JLG
    if(0==slist.size()) throw std::logic_error("empty list");
#endif

    //std::sort(slist.begin(),slist.end(),compstate);

    return;
}

void combinestates(const std::vector<imatrix>&sin,
                   imatrix&sout){
    sout.resize(0);
    for(imatrix s:sin)for(ivector x:s) sout.push_back(x);
    std::sort(sout.begin(),sout.end(),compstate);
    for(auto it=sout.begin();it!=sout.end();){
        if(*it==(*(it-1))) it=sout.erase(it);
        else it++;
    }
}

void changeM(const int M,
             imatrix&statelist){
#pragma omp parallel for
    for(std::size_t i=0;i<statelist.size();i++){
#ifdef TESTING_STATES_JLG
        if (std::abs(M)>statelist[i][1]){
            std::cerr<<"M too large"<<std::endl;
            statelist[i][3]=statelist[i][1]*sign(M);
            continue;
        }
#endif
        statelist[i][3]=M;
    }
}

void transformss(const imatrix&s1,
                 const imatrix&s2,
                 const marray&mass,
                 const int k,
                 sparse&trans){
    //std::cerr<<"transform"<<std::endl;
    trans=sparse(s1.size(),s2.size(),trans.err());
#pragma omp parallel for collapse(2)
    //causes problems with parallel code in Moshinsky.cpp
    for(size_t m=0;m<s1.size();m++) for(size_t n=0;n<s2.size();n++){
        //std::cerr<<m<<','<<n<<'\n';
        if(s1[m][0]!=s2[n][0]) continue;//Pf
        if(s1[m][1]!=s2[n][1]) continue;//J
        if(s1[m][2]!=s2[n][2]) continue;//Pj
        if(s1[m][3]!=s2[n][3]) continue;//M
        if(s1[m][4]!=s2[n][4]) continue;//N
        if(s1[m][5]!=s2[n][5]) continue;//L
        //n1,l1,n2,l2 determine Moshinsy brackets
        if(s1[m][10]!=s2[n][10]) continue;//S
        if(3==s1[m][10])if(s1[m][11] or s2[n][11]) continue;
        //impossible S,Ps
        //std::cerr<<"calc"<<std::endl;
        Double t=spin3(s1[m][10],s1[m][11],s2[n][11],k)
        *tm3(s1[m][6],s1[m][7],s1[m][8],s1[m][9],
             s2[n][6],s2[n][7],s2[n][8],s2[n][9],s1[m][5],
             mass[0],mass[1],mass[2],k);
#pragma omp critical(trans)
        trans.set_entry(m,n,t);
        //std::cerr<<"calc complete"<<std::endl;
    }//m and n
    //std::cerr<<"transform complete"<<std::endl;
    trans.clean();
}

void transform(const imatrix&s1,
               const imatrix&s2,
               const marray&mass,
               const int k,
               sparse&trans){
    transformss(s1,s2,mass,k,trans);
    for(auto it=trans.begin();it!=trans.end();it++)
        if(!s1[it->first[0]][0])
            trans.set_entry(it->first[0],it->first[1],-it->second);
}
