// Conversion between C++ Standard Template Library and GNU Scientific Library
// vectors, matrices (vector<vector>), complex, and sparse matrices
// largely irrelevant since switch to Armadillo for eigenvalues/eigenvectors
//  stl2gsl.h
//  
//
//  Created by Johnathan Gross
//
//

#ifndef stl2gsl_hpp
#define stl2gsl_hpp

#ifdef SPARSE_JLG
#include "sparse/sparse.h"
#include <gsl/gsl_spmatrix.h>
#else
#include "type.h"
#endif
#ifdef _complex_JLG
#include <gsl/gsl_complex.h>
#endif
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

//#ifdef _long_JLG
//#warning "Uses double precision in gsl, not long double."
//#endif

#ifdef _complex_JLG
//convert between STL's complex<> and GSL's gsl_complex
#ifndef _long_JLG
//double
void s2gc(const cdouble&,//sc
          gsl_complex&);//gc
void g2sc(const gsl_complex&,//gc
          cdouble&);//sc
#else
//long double
void s2gc(const cldouble&,//sc
          gsl_complex_long_double&);//gc
void g2sc(const gsl_complex_long_double&,//gc
          cldouble&);//sc
#endif //_long_JLG
#endif

//convert between STL's vector<> and GSL's gsl_vector
//int
void s2gv(const ivector&,//sv
          gsl_vector_int&);//gv
void g2sv(const gsl_vector_int&,//gv
          ivector&);//sv
#ifndef _long_JLG
//double
void s2gv(const dvector&,//sv
          gsl_vector&);//gv
void g2sv(const gsl_vector&,//gv
          dvector&);//sv
#else
//long double
void s2gv(const ldvector&,//sv
          gsl_vector_long_double&);//gv
void g2sv(const gsl_vector_long_double&,//gv
          ldvector&);//sv
#endif
#ifdef _complex_JLG
#ifndef _long_JLG
//complex double
void s2gv(const cdvector&,//sv
          gsl_vector_complex&);//gv
void g2sv(const gsl_vector_complex&,//gv
          cdvector&);//sv
#else
//complex long double
void s2gv(const cldvector&,//sv
          gsl_vector_complex_long_double&);//gv
void g2sv(const gsl_vector_complex_long_double&,//gv
          cldvector&);//sv
#endif
#endif

//convert between STL's vector<vector<>> and GSL's gsl_matrix
//and between my sparse and gsl_matrix
//int
void s2gm(const imatrix&,//sm
          gsl_matrix_int&);//gm
void g2sm(const gsl_matrix_int&,//gm
          imatrix&);//sm
#ifndef _long_JLG
//double
void s2gm(const dmatrix&,//sm
          gsl_matrix&);//gm
void g2sm(const gsl_matrix&,//gm
          dmatrix&);//sm
#else
//long double
void s2gm(const ldmatrix&,//sm
          gsl_matrix_long_double&);//gm
void g2sm(const gsl_matrix_long_double&,//gm
          ldmatrix&);//sm
#endif
#ifdef _complex_JLG
#ifndef _long_JLG
//complex double
void s2gm(const cdmatrix&,//sm
          gsl_matrix_complex&);//gm
void g2sm(const gsl_matrix_complex&,//gm
          cdmatrix&);//sm
#else
//complex long double
void s2gm(const cldmatrix&,//sm
          gsl_matrix_complex_long_double&);//gm
void g2sm(const gsl_matrix_complex_long_double&,//gm
          cldmatrix&);//sm
#endif
#endif

//convert between JLG's sparse and GSL's gsl_matrix
#ifdef SPARSE_JLG
#ifndef _long_JLG
void sp2gpm(const sparse&,//sp
            gsl_spmatrix&);//gp
void gp2spm(const gsl_spmatrix&,//gp
            sparse&);//sp
inline void sp2gm(const sparse&sp,gsl_matrix&g){
    gsl_spmatrix gp;
    sp2gpm(sp,gp);
    g=*gsl_matrix_alloc(gp.size1,gp.size2);
    gsl_spmatrix_sp2d(&g,&gp);
}
inline void s2gpm(const dmatrix&s,gsl_spmatrix&gp){
    sp2gpm(sparse(s),gp);
}
inline void gp2sm(const gsl_spmatrix&gp,dmatrix&s){
    sparse sp;
    gp2spm(gp,sp);
    s=sp.matrix();
}
inline void g2spm(const gsl_matrix&g,sparse&sp){
    gsl_spmatrix gp=*gsl_spmatrix_alloc(g.size1,g.size2);
    gsl_spmatrix_d2sp(&gp,&g);
    gp2spm(gp,sp);
}
#else
void sp2gpm(const sparse&,//sp
            gsl_spmatrix_long_double&);//gp
void gp2spm(const gsl_spmatrix_long_double&,//gp
            sparse&);//sp
inline void sp2gm(const sparse&sp,gsl_matrix_long_double&g){
    gsl_spmatrix_long_double gp;
    sp2gpm(sp,gp);
    g=*gsl_matrix_long_double_alloc(gp.size1,gp.size2);
    gsl_spmatrix_long_double_sp2d(&g,&gp);
}
inline void s2gpm(const ldmatrix&s,gsl_spmatrix_long_double&gp){
    sp2gpm(sparse(s),gp);
}
inline void gp2sm(const gsl_spmatrix_long_double&gp,ldmatrix&s){
    sparse sp;
    gp2spm(gp,sp);
    s=sp.matrix();
}
inline void g2spm(const gsl_matrix_long_double&g,sparse&sp){
    gsl_spmatrix_long_double gp=*gsl_spmatrix_long_double_alloc(g.size1,g.size2);
    gsl_spmatrix_long_double_d2sp(&gp,&g);
    gp2spm(gp,sp);
}
#endif //_long_JLG
#endif //SPARSE_JLG

#endif /* stl2gsl_hpp */
