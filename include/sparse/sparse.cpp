//
//  sparse.cpp
//  
//
//  Created by Johnathan Gross
//

#include "sparse.h"
#include<algorithm>
#include<cmath>
#ifdef TESTING_STATES_JLG
#include <iostream>
#endif
#include <stdexcept>

//constructor
sparse::sparse(const Dmatrix&A,const Double err):
M(A.size()),N(A[0].size()),error(std::max(err,Double(0))),MAP({}){
//#pragma omp parallel for schedule(dynamic)
    for(size_t i=0;i<A.size();i++){
#ifdef TESTING_STATES_JLG
        if(A[i].size()>N) N=A[i].size();
#endif
        for(size_t j=0;j<A[i].size();j++){
            const Double a=A[i][j];
            if(0==error){ if(0==a) continue;}
            else if(std::abs(a)<error) continue;
//#pragma omp critical(sparseconstmat)
            MAP[{i,j}]=a;
        }
    }
}
sparse::sparse(const Dvector&A,const Double err):
M(A.size()),N(A.size()),error(std::max(err,Double(0))),MAP({}){
//#pragma omp parallel for schedule(dynamic)
    for(size_t i=0;i<A.size();i++){
        const Double a=A[i];
        if(0==error){if(0==a) continue;}
        else if(std::abs(a)<error) continue;
//#pragma omp critical(sparseconstvec)
        MAP[{i,i}]=a;
    }
}

//value access
bool sparse::is_diag()const noexcept{
    if(M!=N) return false;
    for(auto it=MAP.begin();it!=MAP.end();it++)
        if(it->first[0]!=it->first[1]) return false;
    return true;
}
bool sparse::is_triu()const noexcept{
    if(M!=N) return false;
    for(auto it=MAP.begin();it!=MAP.end();it++)
        if(it->first[0]>it->first[1]) return false;
    return true;
}
bool sparse::is_tril()const noexcept{
    if(M!=N) return false;
    for(auto it=MAP.begin();it!=MAP.end();it++)
        if(it->first[0]<it->first[1]) return false;
    return true;
}
bool sparse::symmetric()const noexcept{
    if(M!=N) return false;
    for(auto it=MAP.begin();it!=MAP.end();it++){
        try{
            if(std::abs(MAP.at({it->first[1],it->first[0]})-it->second)>error){
#ifdef TESTING_STATES_JLG
                std::cerr<<it->first[0]<<' '<<it->first[1]<<':'<<it->second
                <<' '<<MAP.at({it->first[1],it->first[0]})
                <<' '<<MAP.at({it->first[1],it->first[0]})-it->second
                <<std::endl;
#endif
                return false;
            }
        }catch(std::out_of_range&){
            if(std::abs(it->second)>error){
#ifdef TESTING_STATES_JLG
                std::cerr<<it->first[0]<<' '<<it->first[1]<<':'<<it->second
                <<std::endl;
#endif
                return false;
            }
        }
    }
    return true;
}
bool sparse::has_nan()const noexcept{
    for(auto it=MAP.begin();it!=MAP.end();it++)
        if(std::isnan(it->second)) return true;
    return false;
}
bool sparse::has_inf()const noexcept{
    for(auto it=MAP.begin();it!=MAP.end();it++)
        if(std::isinf(it->second)) return true;
    return false;
}
bool sparse::finite()const noexcept{
    for(auto it=MAP.begin();it!=MAP.end();it++)
        if(!std::isfinite(it->second)) return false;
    return true;
}

//"element access"
Double sparse::at(const sparsekey k)const{
    if(k[0]<M && k[1]<N) return operator[](k);
    else throw std::out_of_range("sDmatrix out of range");
}

//max/min
sparsemap::const_iterator sparse::max()const{
    if(0==size()) return MAP.end();
    auto im=MAP.begin();
    Double m=im->second;
    for(auto it=++MAP.begin();it!=MAP.end();it++)
        if(m<it->second and
           std::numeric_limits<Double>::infinity()>it->second){
            im=it;
            m=im->second;
        }
    return im;
}
sparsemap::const_iterator sparse::min()const{
    if(0==size()) return MAP.end();
    auto im=MAP.begin();
    Double m=im->second;
    for(auto it=++MAP.begin();it!=MAP.end();it++)
        if(m>it->second and
           -std::numeric_limits<Double>::infinity()>it->second){
            im=it;
            m=im->second;
        }
    return im;
}
sparsemap::const_iterator sparse::biggest()const{
    if(0==size()) return MAP.end();
    auto im=MAP.begin();
    Double m=std::abs(im->second);
    for(auto it=++MAP.begin();it!=MAP.end();it++)
        if(m<std::abs(it->second) and
           std::numeric_limits<Double>::infinity()>std::abs(it->second)){
            im=it;
            m=std::abs(im->second);
        }
    return im;
}
sparsemap::const_iterator sparse::smallest()const{
    if(0==size()) return MAP.end();
    auto im=MAP.begin();
    Double m=std::abs(im->second);
    for(auto it=++MAP.begin();it!=MAP.end();it++){
        if(0==it->second) return it;
        if(m>std::abs(it->second)){
            im=it;
            m=std::abs(im->second);
        }
    }
    return im;
}

//manipulate contents
void sparse::set_error(const Double e,const bool c){
    Double t=error;
    if(error!=e) error=e;
    if(t<e or c) clean();
}
void sparse::set_entry(const std::size_t i,const std::size_t j,const Double x){
#ifdef TESTING_STATES_JLG
    if(i>M or j>N) throw std::invalid_argument("index out of range");
#endif
    MAP[{i,j}]=x;
}
void sparse::add_to_entry(const std::size_t i,const std::size_t j,const Double x){
#ifdef TESTING_STATES_JLG
    if(i>M or j>N) throw std::invalid_argument("index out of range");
#endif
    if(MAP.count({i,j})) MAP[{i,j}]+=x;
    else MAP[{i,j}]=x;
}
void sparse::resize(const std::size_t m,const std::size_t n){
    if(m<M or n<N) for(auto it=MAP.begin();it!=MAP.end();){
        if(it->first[0]>m or it->first[1]>n) it=MAP.erase(it);
        else it++;
    }
        M=m;
        N=n;
}
void sparse::clean(){
    for(auto it=MAP.begin();it!=MAP.end();){
        if(std::abs(it->second)<error) it=MAP.erase(it);
        else if(0==it->second) it=MAP.erase(it);
        else it++;
    }
}

//return modified matrices
sparse sparse::transpose()const{
    sparse result(N,M);
    for(auto it=MAP.begin();it!=MAP.end();it++)
        result.MAP[{it->first[1],it->first[0]}]=it->second;
    return result;
}
Dmatrix sparse::matrix()const{
    Dmatrix result(M,Dvector(N,0));
//#pragma omp parallel for
    for(auto it=MAP.begin();it!=MAP.end();it++){
        result[it->first[0]][it->first[1]]=it->second;
    }
    return result;
}

//submatrix
sparse sparse::sub(const std::size_t ri,
                   const std::size_t rf,
                   const std::size_t ci,
                   const std::size_t cf)const{
#ifdef TESTING_STATES_JLG
    if(ri>rf) throw std::invalid_argument("first row after last row\n");
    if(ci>cf) throw std::invalid_argument("first collumn after last collumn\n");
    if(rf>M) throw std::out_of_range("row out of range\n");
    if(cf>N) throw std::out_of_range("collumn out of range\n");
#else
    if(ri>rf) return sparse(0,0,error);
    if(ci>cf) return sparse(0,0,error);
    if(rf>M) return sub(ri,M,ci,cf);
    if(cf>N) return sub(ri,rf,ci,N);
#endif
    sparse temp(rf-ri,cf-ci,error);
    if(temp.M and temp.N) for(auto i=MAP.begin();i!=MAP.end();i++){
        if(i->first[0]<ri) continue;
        if(i->first[0]>=rf) continue;
        if(i->first[1]<ci) continue;
        if(i->first[1]>=cf) continue;
        temp.MAP[{i->first[0]-ri,i->first[1]-ci}]=i->second;
    }
    return temp;
}

//operator overloads
sparse operator-(const sparse&A){
    sparse result(A);
//#pragma omp parallel for
    for(auto it=result.MAP.begin();it!=result.MAP.end();it++)
        result.MAP[it->first]=-it->second;
    return result;
}
sparse operator+(const sparse&A,const sparse&B){
    if(A.M!=B.M and A.N!=B.N)
        throw(std::invalid_argument("Matrices not same size"));
    sparse result(A);
//#pragma omp parallel for
    for(auto it=B.MAP.begin();it!=B.MAP.end();it++){
        if(result.MAP.count(it->first)) result.MAP[it->first]+=it->second;
        else{
//#pragma omp critical(spareadd)
            result.MAP[it->first]=it->second;
        }
    }
    if(B.error>result.error) result.error=B.error;
    return result;
}
sparse operator-(const sparse&A,const sparse&B){
    if(A.M!=B.M and A.N!=B.N) throw(std::invalid_argument("Matrices not same size"));
    if(A.MAP==B.MAP) return sparse(A.M,A.N,std::max(A.error,B.error));
    sparse result(A);
//#pragma omp parallel for
    for(auto it=B.MAP.begin();it!=B.MAP.end();it++){
        if(result.MAP.count(it->first)) result.MAP[it->first]-=it->second;
        else{
//#pragma omp critical(spareadd)
            result.MAP[it->first]=-it->second;
        }
    }
    if(B.error>result.error) result.error=B.error;
    return result;
}
sparse operator*(const sparse&A,const Double&d){
    if(0==d) return sparse(A.M,A.N,A.error);
    sparse result(A);
    for(auto it=result.MAP.begin();it!=result.MAP.end();it++)
        it->second*=d;
    return result;
}
sparse operator/(const sparse&A,const Double&d){
    sparse result(A);
    for(auto it=result.MAP.begin();it!=result.MAP.end();it++)
        it->second/=d;
    return result;
}
/*sparse operator*(const sparse&A,const sparse&B){
    if(A.N!=B.M) throw(std::invalid_argument("Inner dimensions not same size"));
    sparse result(A.M,B.N,std::max(A.error,B.error));
    sparse BT=B.transpose();
    size_t K=A.N;

#pragma omp parallel for
    for(size_t i=0;i<A.M;i++){
        const auto Ai=A.MAP.lower_bound({i,0});
        const auto Af=A.MAP.upper_bound({i,K});
        if(Af==Ai) continue;
        for(size_t j=0;j<B.N;j++){
            auto x=Ai;
            auto y=BT.MAP.lower_bound({j,0});
            const auto Bf=BT.MAP.upper_bound({j,K});
            if(Bf==y) continue;
            Double AB=0;
            do{
                if(x->first[1]<y->first[1])
                    x++;
                else if(y->first[1]<x->first[1])
                    y++;
                else{
                    AB+=x->second*y->second;
                    x++;
                    y++;
                }
            }while(x!=Af and y!=Bf);
#pragma omp critical(sparemult)
            result.MAP[sparsekey{i,j}]=AB;
        }
    }

    return result;
}*/
sparse operator/(const sparse&A,const sparse&B){
    if(A.N!=B.N) throw(std::invalid_argument("Inner dimensions not same size"));
    sparse result(A.M,B.M,std::max(A.error,B.error));
    size_t K=A.N;

#pragma omp parallel for
    for(size_t i=0;i<A.M;i++){
        const auto Ai=A.MAP.lower_bound({i,0});
        const auto Af=A.MAP.upper_bound({i,K});
        if(Af==Ai) continue;
        for(size_t j=0;j<B.M;j++){
            auto x=Ai;
            auto y=B.MAP.lower_bound({j,0});
            const auto Bf=B.MAP.upper_bound({j,K});
            if(Bf==y) continue;
            Double AB=0;
            do{
                if(x->first[1]<y->first[1])
                    x++;
                else if(y->first[1]<x->first[1])
                    y++;
                else{
                    AB+=x->second*y->second;
                    x++;
                    y++;
                }
            }while(x!=Af and y!=Bf);
#pragma omp critical(sparemult)
            result.MAP[sparsekey{i,j}]=AB;
        }
    }

    return result;
}
sparse operator*(const Dvector&B,const sparse&A){
    sparse result(A);
    for(auto x=result.begin();x!=result.end();x++)
        x->second*=B[x->first[0]];
    return result;
}
sparse operator*(const sparse&A,const Dvector&B){
    sparse result(A);
    for(auto x=result.begin();x!=result.end();x++)
        x->second*=B[x->first[1]];
    return result;
}

sparse& operator+=(sparse&A,const sparse&B){
    if(A.M!=B.M and A.N!=B.N)
        throw(std::invalid_argument("Matrices not same size"));
    //#pragma omp parallel for
    for(auto it=B.MAP.begin();it!=B.MAP.end();it++){
        if(A.MAP.count(it->first)) A.MAP[it->first]+=it->second;
        else{
            //#pragma omp critical(spareadd)
            A.MAP[it->first]=it->second;
        }
    }
    if(B.error>A.error) A.error=B.error;
    return A;
}
sparse& operator-=(sparse&A,const sparse&B){
    if(A.M!=B.M and A.N!=B.N)
        throw(std::invalid_argument("Matrices not same size"));
    if(&A==&B or B.MAP==A.MAP){
        A.zero();
        return A;
    }
    //#pragma omp parallel for
    for(auto it=B.MAP.begin();it!=B.MAP.end();it++){
        if(A.MAP.count(it->first)) A.MAP[it->first]-=it->second;
        else{
            //#pragma omp critical(spareadd)
            A.MAP[it->first]=-it->second;
        }
    }
    if(B.error>A.error) A.error=B.error;
    return A;
}
sparse& operator*=(sparse&A,const Double&d){
    if(0==d){
        A.zero();
        return A;
    }
    for(auto it=A.MAP.begin();it!=A.MAP.end();it++)
        it->second*=d;
    return A;
}
sparse& operator/=(sparse&A,const Double&d){
    for(auto it=A.MAP.begin();it!=A.MAP.end();it++)
        it->second/=d;
    return A;
}
sparse& operator*=(sparse&A,const Dvector&B){
    for(auto x=A.begin();x!=A.end();x++)
        x->second*=B[x->first[0]];
    return A;
}
sparse& operator/=(sparse&A,const Dvector&B){
    for(auto x=A.begin();x!=A.end();x++)
        x->second*=B[x->first[1]];
    return A;
}
bool operator==(const sparse&A,const sparse&B){
    if(A.M!=B.M or A.N!=B.N) return false;

    const Double e=std::min(A.error,B.error);
    auto x=A.MAP.cbegin(),
    Af=A.MAP.cend();
    auto y=B.MAP.cbegin(),
    Bf=B.MAP.cend();

    while(x!=Af and y!=Bf){
        if(x->first==y->first){
            if(x->second==y->second){
                x++;
                y++;
                continue;
            }else if(e and std::abs(x->second-y->second)<=e){
                x++;
                y++;
                continue;
            }else return false;
        }else{//x->first!=y->first
            if(x->first<y->first){
                if(!A.error) return false;
                else if(std::abs(x->second)<A.error){
                    x++;
                    continue;
                }else return false;
            }else if(x->first>y->first){
                if(!B.error) return false;
                else if(std::abs(y->second)<B.error){
                    y++;
                    continue;
                }else return false;
            }else throw std::logic_error("keys not equal, less, or greater");
        }
    }
    return true;
}

sparse elem_mult(const sparse&A,const sparse&B){
    if(A.M!=B.M and A.N!=B.N)
        throw(std::invalid_argument("Matrices not same size"));
    sparse result(A.M,A.N);
    auto x=A.MAP.cbegin(),
    Af=A.MAP.cend();
    auto y=B.MAP.cbegin(),
    Bf=B.MAP.cend();

    while(x!=Af and y!=Bf){
        if(x->first==y->first){
            result.MAP[x->first]=x->second*y->second;
            x++;
            y++;
        }else{//x->first!=y->first
            if(x->first<y->first){
                //0 so do nothing
                x++;
            }else if(x->first>y->first){
                //0, so do nothing
                y++;
            }else throw std::logic_error("keys not equal, less, or greater");
        }
    }
    return result;
}
sparse elem_div(const sparse&A,const sparse&B){
    if(A.M!=B.M and A.N!=B.N)
        throw(std::invalid_argument("Matrices not same size"));
    sparse result(A.M,A.N);
    auto x=A.MAP.cbegin(),
    Af=A.MAP.cend();
    auto y=B.MAP.cbegin(),
    Bf=B.MAP.cend();

    while(x!=Af and y!=Bf){
        if(x->first==y->first){
            result.MAP[x->first]=x->second/y->second;
            x++;
            y++;
        }else{//x->first!=y->first
            if(x->first<y->first){
                result.MAP[x->first]=std::numeric_limits<Double>::quiet_NaN();
                x++;
            }else if(x->first>y->first){
                //0, so do nothing
                y++;
            }else throw std::logic_error("keys not equal, less, or greater");
        }
    }
    return result;

}

//diagonal sparse matrix
sparse diagonal(const std::size_t N,const Double d,const Double e){
    if(0==e){if(0==d) return sparse(N,N);}
    else if(std::abs(d)<e) return sparse(N,N,e);
    return sparse(Dvector(N,d),e);
}
//transmutation matrix
sparse transmutation(const tvector&I){
    sparse T(I.size(),*std::max_element(I.begin(),I.end()));
    for(size_t i=0;i<I.size();i++) T.set_entry(i,I[i],1);
    return T;
}
sparse transmutation(const tvector&I,std::size_t t){
    if(t<*std::max_element(I.begin(),I.end()))
        throw std::invalid_argument("max value larger than size");
    sparse T(I.size(),t);
    for(size_t i=0;i<I.size();i++) T.set_entry(i,I[i],1);
    return T;
}
