//
//  sparse.h
//  
//
//  Created by Johnathan Gross
//

#include "type.h"
#include <array>
#include <map>

#ifndef sparse_h
#define sparse_h

#ifndef SPARSE_JLG
#define SPARSE_JLG
#endif

typedef std::array<std::size_t,2> sparsekey;
typedef std::map<sparsekey,Double> sparsemap;

class sparse{
private:
    std::size_t M,N;
    Double error;
    sparsemap MAP;
public:
    //constructors
    sparse(const std::size_t M_=0,const std::size_t N_=0,const Double e=0):
    	M(M_),N(N_),error(std::max(Double(0),e)),MAP({}){}
    sparse(const Dmatrix&,const Double=0);
    sparse(const Dvector&,const Double=0);
    sparse(const sparse&A,const Double e):
    M(A.M),N(A.N),error(std::max(Double(0),e)),MAP(A.MAP){clean();}

    //value access
    inline std::size_t m()const noexcept{return M;}
    inline std::size_t n()const noexcept{return N;}
    inline Double err()const noexcept{return error;}
    inline bool empty()const noexcept{return MAP.empty();}
    inline auto size()const noexcept{return MAP.size();}
    inline auto max_size()const noexcept{return MAP.max_size();}
    inline const sparsemap& map()noexcept{return MAP;}
    inline auto begin()noexcept{return MAP.begin();}
    inline auto cbegin()const noexcept{return MAP.cbegin();}
    inline auto rbegin()noexcept{return MAP.rbegin();}
    inline auto crbegin()const noexcept{return MAP.crbegin();}
    inline auto end()noexcept{return MAP.end();}
    inline auto cend()const noexcept{return MAP.cend();}
    inline auto rend()noexcept{return MAP.rend();}
    inline auto crend()const noexcept{return MAP.crend();}
    inline bool is_square()const noexcept{return M==N;}
    bool is_diag()const noexcept;
    bool is_triu()const noexcept;
    bool is_tril()const noexcept;
    bool symmetric()const noexcept;
    bool has_nan()const noexcept;
    bool has_inf()const noexcept;
    bool finite()const noexcept;

    //"element access", only gives value, does not allow access
    inline Double operator[](const sparsekey k)const noexcept{
        if(MAP.count(k)) return MAP.at(k);
        else return 0;
    }
    Double at(const sparsekey k)const;
    inline Double at(const std::size_t i,std::size_t j)const{return at({i,j});}
    inline auto count(sparsekey k)const noexcept{return MAP.count(k);}
    inline auto count(const std::size_t i,const std::size_t j)const noexcept{return count({i,j});}
    inline Double operator()(const sparsekey k)const{return at(k);}
    inline Double operator()(const std::size_t i,const std::size_t j)const{return at({i,j});}

    //max/min
    sparsemap::const_iterator max()const;
    inline Double max(int)const{return max()->second;}
    sparsemap::const_iterator min()const;
    inline Double min(int)const{return min()->second;}
    sparsemap::const_iterator biggest()const;
    inline Double biggest(int)const{return biggest()->second;}
    sparsemap::const_iterator smallest()const;
    inline Double smallest(int)const{return smallest()->second;}

    //functions that manipulate the content
    void set_error(const Double,const bool=false);//true for clean
    void set_entry(const std::size_t,const std::size_t,const Double);
    void add_to_entry(const std::size_t,const std::size_t,const Double);
    void resize(const std::size_t,const std::size_t);
    void clean();
    //removes any elements smaller than error
    inline void zero(){MAP.clear();}
    //sets all values to 0

    //functions that return modifications of the sparse matrix
    sparse transpose()const;
    //returns transpose
    Dmatrix matrix()const;
    //returns matrix with zeros filled in
    sparse sub(const std::size_t,//ri
               const std::size_t,//rf
               const std::size_t,//ci
               const std::size_t)const;//cf
    //returns submatrix [ri,rf)[ci,cf)
    sparse sub(const std::size_t r,
               const std::size_t c)const{
    //overload with initial being 0
        return sub(0,r,0,c);
    }

    //implicit conversion to Dmatrix
    operator Dmatrix()const{return matrix();}

    //operator overloads
    inline friend sparse operator+(const sparse&A){return A;}
    friend sparse operator-(const sparse&);
    friend sparse operator+(const sparse&,const sparse&);
    friend sparse operator-(const sparse&,const sparse&);
    friend sparse operator*(const sparse&,const Double&);
    inline friend sparse operator*(const Double&a,const sparse&S){return S*a;}
    friend sparse operator/(const sparse&,const Double&);
    friend sparse operator/(const sparse&,const sparse&);
    //multiplication by transpose of second matrix
    inline friend sparse operator*(const sparse&A,const sparse&B){
        return A/B.transpose();
    }
    friend sparse operator*(const Dvector&,const sparse&);
    //left multiplication by diagonal matrix
    friend sparse operator*(const sparse&,const Dvector&);
    //right multiplication by diagonal matrix

    friend sparse& operator+=(sparse&,const sparse&);
    friend sparse& operator-=(sparse&,const sparse&);
    inline friend sparse& operator*=(sparse&A,const sparse&B){return A=A*B;}
    friend sparse& operator*=(sparse&,const Double&);
    inline friend sparse& operator/=(sparse&A,const sparse&B){return A=A/B;}
    friend sparse& operator/=(sparse&,const Double&);
    friend sparse& operator*=(sparse&,const Dvector&);
    //left multiplication by diagonal matrix
    friend sparse& operator/=(sparse&,const Dvector&);
    //right multiplication by diagonal matrix
    friend bool operator==(const sparse&,const sparse&);
    inline friend bool operator!=(const sparse&A,const sparse&B){
        return !(A==B);
    }

    //element-wise multiplication/division
    friend sparse elem_mult(const sparse&,const sparse&);
    friend sparse elem_div(const sparse&,const sparse&);
};

sparse diagonal(const std::size_t,const Double=1,const Double=0);
//creates a sparse matrix that is (second arg)*identity matrix, with error third argument
sparse transmutation(const tvector&);
sparse transmutation(const tvector&,const std::size_t);
//creates a sparse matrix that rearranges elements according to tvector

#endif /* sparse_h */
