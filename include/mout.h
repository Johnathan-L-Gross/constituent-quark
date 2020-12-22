//
//  mout.h
//  
//
//  Created by Johnathan Gross on 1/22/19.
//

#ifndef matrix_h
#define matrix_h

#ifdef SPARSE_JLG
#include "sparse/sparse.h"
#else
#include "type.h"
#endif
#include <iostream>
#include <array>
#include <algorithm>

template<class T>
inline void dispvec(const std::vector<T>&v,
                    bool orient=false,
                    std::ostream&o=std::cout){
    const char c=(orient?'\n':' ');
    for(std::size_t i=0;i<v.size();i++)
        o<<v[i]<<c;
    o<<std::endl;
}
inline void dispvec(const bvector&v,
                    bool orient=false,
                    std::ostream&o=std::cout){
    for(std::size_t i=0;i<v.size();i++){
        o<<v[i];
        if(orient) o<<'\n';
    }
    o<<std::endl;
}
template<class T,std::size_t N>
inline void disparray(const std::array<T,N>&a,
                      bool orient=false,
                      std::ostream&o=std::cout){
    const char c=(orient?'\n':' ');
    for(std::size_t i=0;i<N;i++)
        o<<a[i]<<c;
    o<<std::endl;
}

template<class T,size_t N>
inline void dispvec(const std::vector<T>&v,
                    std::array<size_t,N> l,
                    bool orient=false,
                    std::ostream&o=std::cout){
    //no bounds check
    const char c=(orient?'\n':' ');
    for(std::size_t i:l)
        o<<v[i]<<c;
    o<<std::endl;
}
template<std::size_t N>
inline void dispvec(const bvector&v,
                    std::array<size_t,N> l,
                    bool orient=false,
                    std::ostream&o=std::cout){
    //no bounds check
    for(std::size_t i:l){
        o<<v[i];
        if(orient) o<<'\n';
    }
    o<<std::endl;
}

template<class T>
inline void dispvec(const std::vector<T>&v,
                    std::size_t l,
                    bool orient=false,
                    std::ostream&o=std::cout){
    //no bounds check
    const char c=(orient?'\n':' ');
    for(std::size_t i=0;i<l;i++)
        o<<v[i]<<c;
    o<<std::endl;
}
inline void dispvec(const bvector&v,
                    std::size_t l,
                    bool orient=false,
                    std::ostream&o=std::cout){
    //no bounds check
    for(std::size_t i=0;i<l;i++){
        o<<v[i];
        if(orient) o<<'\n';
    }
    o<<std::endl;
}

template<class T>
inline void dispmat(const matrix<T>&d,
                    const bool inv=false,
                    std::ostream&o=std::cout){
    if(!inv){
        for(std::size_t i=0;i<d.size();i++){
            for(std::size_t j=0;j<d[i].size();j++)
                o<<d[i][j]<<' ';
            o<<'\n';
        }
    }else{
        std::size_t maxj=d[0].size();
#pragma omp parallel for reduction(max:maxj)
        for(std::size_t i=1;i<d.size();i++) maxj=d[i].size();
        for(std::size_t j;j<maxj;j++){
            for(std::size_t i;i<d.size();i++){
                if(d[i].size()>j) o<<d[i][j]<<' ';
                else o<<" x ";
            }
            o<<'\n';
        }
    }
    o<<std::endl;
}
inline void dispmat(const bmatrix&d,
                    const bool inv=false,
                    std::ostream&o=std::cout){
    if(!inv){
        for(std::size_t i=0;i<d.size();i++){
            for(std::size_t j=0;j<d[i].size();j++)
                o<<d[i][j];
            o<<'\n';
        }
    }else{
        std::size_t maxj=d[0].size();
#pragma omp parallel for reduction(max:maxj)
        for(std::size_t i=1;i<d.size();i++) maxj=d[i].size();
        for(std::size_t j;j<maxj;j++){
            for(std::size_t i;i<d.size();i++){
                if(d[i].size()>j) o<<d[i][j];
                else o<<'x';
            }
            o<<'\n';
        }
    }
    o<<std::endl;
}

#ifdef SPARSE_JLG
inline void dispmat(const sparse&d,
                    const bool inv=false,
                    std::ostream&o=std::cout){
    for(std::size_t i=0;i<(inv?d.n():d.m());i++){
        for(std::size_t j=0;j<(inv?d.m():d.n());j++)
            o<<d[(inv?sparsekey{j,i}:sparsekey{i,j})]<<' ';
        o<<'\n';
    }
    o<<std::endl;
}
#endif

template<class T,std::size_t N>
inline void dispcol(const matrix<T>&d,
                    const std::array<std::size_t,N> l,
                    std::ostream&o=std::cout){
    //no bounds check
    for(std::size_t i=0;i<d.size();i++){
        for(std::size_t j:l)
            o<<d[i][j]<<' ';
        o<<'\n';
    }
    o<<std::endl;
}
template<std::size_t N>
inline void dispcol(const bmatrix&d,
                    const std::array<std::size_t,N> l,
                    std::ostream&o=std::cout){
    //no bounds check
    for(std::size_t i=0;i<d.size();i++){
        for(std::size_t j:l)
            o<<d[i][j];
        o<<'\n';
    }
    o<<std::endl;
}

template<class T>
inline void dispcol(const matrix<T>&d,
                    const std::size_t j,
                    const bool orient=true,
                    std::ostream&o=std::cout){
    //no bounds check
    const char c=(orient?'\n':' ');
    for(std::size_t i=0;i<d.size();i++)
        o<<d[i][j]<<c;
    o<<std::endl;
}
inline void dispcol(const bmatrix&d,
                    const std::size_t j,
                    const bool orient=true,
                    std::ostream&o=std::cout){
    //no bounds check
    for(std::size_t i=0;i<d.size();i++){
        o<<d[i][j];
        if(orient) o<<'\n';
    }
    o<<std::endl;
}

template<class T,size_t N>
inline void disprow(const matrix<T>&d,
                    const std::array<size_t,N> l,
                    std::ostream&o=std::cout){
    //no bounds check
    for(std::size_t i:l){
        for(std::size_t j=0;j<d[i].size();j++)
            o<<d[i][j]<<' ';
        o<<'\n';
    }
    o<<std::endl;
}
template<std::size_t N>
inline void disprow(const bmatrix&d,
                    const std::array<size_t,N> l,
                    std::ostream&o=std::cout){
    //no bounds check
    for(std::size_t i:l){
        for(std::size_t j=0;j<d[i].size();j++)
            o<<d[i][j];
        o<<'\n';
    }
    o<<std::endl;
}

template<class T>
inline void disprow(const matrix<T>&d,
                    const std::size_t i,
                    std::ostream&o=std::cout){
    //no bounds check
    for(std::size_t j=0;j<d[i].size();j++)
        o<<d[i][j]<<' ';
    o<<std::endl;
}
inline void disprow(const bmatrix&d,
                    const std::size_t i,
                    std::ostream&o=std::cout){
    //no bounds check
    for(std::size_t j=0;j<d[i].size();j++)
        o<<d[i][j];
    o<<std::endl;
}

template<class T>
inline void disprowcol(const matrix<T>&d,
                       const std::vector<std::size_t> m,
                       const std::vector<std::size_t> n,
                       std::ostream&o=std::cout){
    //no bounds check
    for(std::size_t i:m){
        for(std::size_t j:n)
            o<<d[i][j]<<' ';
        o<<'\n';
    }
    o<<std::endl;
}
inline void disprowcol(const bmatrix&d,
                       const std::vector<std::size_t> m,
                       const std::vector<std::size_t> n,
                       std::ostream&o=std::cout){
    //no bounds check
    for(std::size_t i:m){
        for(std::size_t j:n)
            o<<d[i][j];
        o<<'\n';
    }
    o<<std::endl;
}

template<class T>
inline void dispsparse(const matrix<T>&d,
                       const bool inv=false,
                       T lim=0,
                       std::ostream&o=std::cout){
    if(lim<0) lim*=-1;
    if(!inv){
        for(std::size_t i=0;i<d.size();i++)
            for(std::size_t j=0;j<d[i].size();j++)
                if(std::abs(d[i][j])>lim) o<<i<<' '<<j<<' '<<d[i][j]<<std::endl;
    }else{
        std::size_t i=0,j=0;
        for(;j<d[i].size();j++){
            for(;i<d.size();i++){
                if(std::abs(d[i][j])>lim) o<<j<<' '<<i<<' '<<d[i][j]<<std::endl;
            }
            i=0;
        }
    }
    o<<std::endl;
}

#ifdef SPARSE_JLG
inline void dispsparse(const sparse&d,
                       const bool inv=false,
                       std::ostream&o=std::cout){
    sparse dd=(inv?d.transpose():d);
    for(auto it=dd.begin();it!=dd.end();it++)
        o<<it->first[0]<<' '<<it->first[1]<<' '<<it->second<<'\n';
    o<<std::endl;
}

inline void wavedisp(const sparse&W,
              std::size_t s=-1,//number of states
              std::size_t k=-1,//number of coefficients
              std::ostream&o=std::cout){
    s=std::min(s,W.n());
    k=std::min(k,W.m());
    //std::cerr<<s<<' '<<k<<std::endl;
    for(std::size_t i=0;i<s;i++){
        o<<i<<": ";
        tvector T(W.m());
#pragma omp parallel for
        for(std::size_t j=0;j<W.m();j++) T[j]=j;
        std::sort(T.begin(),T.end(),
                  [&W,i](const std::size_t&m,const std::size_t&n){
                      return std::abs(W(m,i))>std::abs(W(n,i));
                  });
        for(size_t j=0;j<k;j++)
            if(W.count(T[j],i))
                o<<T[j]<<' '<<W(T[j],i)<<", ";
            else break;
        o<<"\b\b \n";
    }
    o<<std::endl;
}
#endif //SPARSE_JLG
#endif /* matrix_h */
