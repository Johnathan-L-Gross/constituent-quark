//
//  stl2gsl.cpp
//  
//
//  Created by Johnathan Gross
//
//

#include "stl2gsl.h"
#ifdef _complex_JLG
#include <gsl/gsl_complex_math.h>
#endif
#ifdef TESTING_STATES_JLG
#include <iostream>
#endif

//******************************************************************************
//Complex
#ifdef _complex_JLG
#ifndef _long_JLG
//double
void s2gc(const cdouble&sc,gsl_complex&gc){
    GSL_SET_COMPLEX(&gc,sc.real(),sc.imag());
}
void g2sc(const gsl_complex&gc,cdouble&sc){
    sc=cdouble(GSL_REAL(gc),GSL_IMAG(gc));
}
#else
//long double
void s2gc(const cldouble&sc,gsl_complex_long_double&gc){
    GSL_SET_COMPLEX(&gc,sc.real(),sc.imag());
}
void g2sc(const gsl_complex_long_double&gc,cldouble&sc){
    sc=cldouble(GSL_REAL(gc),GSL_IMAG(gc));
}
#endif //_long_JLG
#endif //_complex_JLG

//******************************************************************************
//Vector
//int
void s2gv(const ivector&sv,gsl_vector_int&gv){
    gv=*gsl_vector_int_calloc(sv.size());
#pragma omp parallel for
    for(std::size_t i=0;i<sv.size();i++) gsl_vector_int_set(&gv,i,sv[i]);
}
void g2sv(const gsl_vector_int&gv,ivector&sv){
    sv=ivector(gsl_vector_int_const_ptr(&gv,0),gsl_vector_int_const_ptr(&gv,0)+gv.size);
}
#ifndef _long_JLG
//double
void s2gv(const dvector&sv,gsl_vector&gv){
        gv=*gsl_vector_calloc(sv.size());
    #pragma omp parallel for
        for(std::size_t i=0;i<sv.size();i++) gsl_vector_set(&gv,i,sv[i]);
}
void g2sv(const gsl_vector&gv,dvector&sv){
    sv=dvector(gsl_vector_const_ptr(&gv,0),gsl_vector_const_ptr(&gv,0)+gv.size);
}
#else
//long double
void s2gv(const ldvector&sv,gsl_vector_long_double&gv){
        gv=*gsl_vector_long_double_calloc(sv.size());
    #pragma omp parallel for
        for(std::size_t i=0;i<sv.size();i++) gsl_vector_long_double_set(&gv,i,sv[i]);
}
void g2sv(const gsl_vector_long_double&gv,ldvector&sv){
    sv=ldvector(gsl_vector_long_double_const_ptr(&gv,0),
                gsl_vector_long_double_const_ptr(&gv,0)+gv.size);
}
#endif //_long_JLG
#ifdef _complex_JLG
#ifndef _long_JLG
//complex double
void s2gv(const cdvector&sv,gsl_vector_complex&gv){
    gv=*gsl_vector_complex_calloc(sv.size());
#pragma omp parallel for
    for(size_t i=0;i<sv.size();i++){
        gsl_complex gc;
        s2gc(sv[i],gc);
        gsl_vector_complex_set(&gv,i,gc);
    }
}
void g2sv(const gsl_vector_complex&gv,cdvector&sv){
    sv.resize(gv.size,0);
#pragma omp parallel for
    for(std::size_t i=0;i<gv.size;i++)
        g2sc(gsl_vector_complex_get(&gv,i),sv[i]);
}
#else
//complex long double
void s2gv(const cldvector&sv,gsl_vector_complex_long_double&gv){
    gv=*gsl_vector_complex_long_double_calloc(sv.size());
#pragma omp parallel for
    for(size_t i=0;i<sv.size();i++){
        gsl_complex_long_double gc;
        s2gc(sv[i],gc);
        gsl_vector_complex_long_double_set(&gv,i,gc);
    }
}
void g2sv(const gsl_vector_complex_long_double&gv,cldvector&sv){
    sv.resize(gv.size,0);
#pragma omp parallel for
    for(std::size_t i=0;i<gv.size;i++)
        g2sc(gsl_vector_complex_long_double_get(&gv,i),sv[i]);
}
#endif //_long_JLG
#endif //_complex_JLG

//******************************************************************************
//Matrix
//int
void s2gm(const imatrix&sm,gsl_matrix_int&gm){
    std::size_t s=0;
#pragma omp parallel for reduction(max:s)
    for(size_t i=0;i<sm.size();i++) s=std::max(s,sm[i].size());
    gm=*gsl_matrix_int_calloc(sm.size(),s);
#pragma omp parallel for
    for(std::size_t i=0;i<sm.size();i++)for(std::size_t j=0;j<sm[i].size();j++)
        gsl_matrix_int_set(&gm,i,j,sm[i][j]);
}
void g2sm(const gsl_matrix_int&gm,imatrix&sm){
    sm.resize(gm.size1,ivector(gm.size2,0));
#pragma omp parallel for
    for(std::size_t i=0;i<sm.size();i++)//for(std::size_t j=0;j<sm[0].size();j++)
        sm[i]=ivector(gsl_matrix_int_const_ptr(&gm,i,0),
                      gsl_matrix_int_const_ptr(&gm,i,0)+gm.size2);
}
#ifndef _long_JLG
//double
void s2gm(const dmatrix&sm,gsl_matrix&gm){
    std::size_t s=0;
#pragma omp parallel for reduction(max:s)
    for(size_t i=0;i<sm.size();i++) s=std::max(s,sm[i].size());
    gm=*gsl_matrix_calloc(sm.size(),s);
#pragma omp parallel for
    for(std::size_t i=0;i<sm.size();i++)for(std::size_t j=0;j<sm[i].size();j++)
        gsl_matrix_set(&gm,i,j,sm[i][j]);
}
void g2sm(const gsl_matrix&gm,dmatrix&sm){
    sm.resize(gm.size1,dvector(gm.size2,0));
#pragma omp parallel for
    for(std::size_t i=0;i<sm.size();i++)
        sm[i]=dvector(gsl_matrix_const_ptr(&gm,i,0),
                      gsl_matrix_const_ptr(&gm,i,0)+gm.size2);
}
#else
//long double
void s2gm(const ldmatrix&sm,gsl_matrix_long_double&gm){
    std::size_t s=0;
#pragma omp parallel for reduction(max:s)
    for(size_t i=0;i<sm.size();i++) s=std::max(s,sm[i].size());
    gm=*gsl_matrix_long_double_calloc(sm.size(),s);
#pragma omp parallel for
    for(std::size_t i=0;i<sm.size();i++)for(std::size_t j=0;j<sm[i].size();j++)
        gsl_matrix_long_double_set(&gm,i,j,sm[i][j]);
}
void g2sm(const gsl_matrix_long_double&gm,ldmatrix&sm){
    sm.resize(gm.size1,ldvector(gm.size2,0));
#pragma omp parallel for
    for(std::size_t i=0;i<sm.size();i++)
        sm[i]=ldvector(gsl_matrix_long_double_const_ptr(&gm,i,0),
                       gsl_matrix_long_double_const_ptr(&gm,i,0)+gm.size2);
}
#endif
#ifdef _complex_JLG
#ifndef _long_JLG
//complex double
void s2gm(const cdmatrix&sm,gsl_matrix_complex&gm){
    std::size_t s=0;
#pragma omp parallel for reduction(max:s)
    for(size_t i=0;i<sm.size();i++) s=std::max(s,sm[i].size());
    gm=*gsl_matrix_complex_calloc(sm.size(),s);
#pragma omp parallel for collapse(2)
    for(std::size_t i=0;i<sm.size();i++)for(std::size_t j=0;j<s;j++){
        gsl_complex gc;
        s2gc(sm[i][j],gc);
        gsl_matrix_complex_set(&gm,i,j,gc);
    }
}
void g2sm(const gsl_matrix_complex&gm,cdmatrix&sm){
    sm.resize(gm.size1,cdvector(gm.size2,0));
#pragma omp parallel for collapse(2)
    for(std::size_t i=0;i<sm.size();i++) for(std::size_t j=0;j<gm.size2;j++)
        g2sc(gsl_matrix_complex_get(&gm,i,j),sm[i][j]);
}
#else
//complex long double
void s2gm(const cldmatrix&sm,gsl_matrix_complex_long_double&gm){
    std::size_t s=0;
#pragma omp parallel for reduction(max:s)
    for(size_t i=0;i<sm.size();i++) s=std::max(s,sm[i].size());
    gm=*gsl_matrix_complex_long_double_calloc(sm.size(),s);
#pragma omp parallel for collapse(2)
    for(std::size_t i=0;i<sm.size();i++)for(std::size_t j=0;j<s;j++){
        gsl_complex_long_double gc;
        s2gc(sm[i][j],gc);
        gsl_matrix_complex_long_double_set(&gm,i,j,gc);
    }
}
void g2sm(const gsl_matrix_complex_long_double&gm,cldmatrix&sm){
    sm.resize(gm.size1,cldvector(gm.size2,0));
#pragma omp parallel for collapse(2)
    for(std::size_t i=0;i<sm.size();i++) for(std::size_t j=0;j<gm.size2;j++)
        g2sc(gsl_matrix_complex_long_double_get(&gm,i,j),sm[i][j]);
}
#endif //_long_JLG
#endif //_complex_JLG

//******************************************************************************
//Matrix
#ifdef SPARSE_JLG
#ifndef _long_JLG
//double
void sp2gpm(const sparse&sp,gsl_spmatrix&gp){
    gp=*gsl_spmatrix_alloc_nzmax(sp.m(),sp.n(),sp.size(),GSL_SPMATRIX_COO);
    for(auto it=sp.cbegin();it!=sp.cend();it++)
        gsl_spmatrix_set(&gp,it->first[0],it->first[1],it->second);
}
void gp2spm(const gsl_spmatrix&gp,sparse&sp){
    sp.zero();
    sp.resize(gp.size1,gp.size2);
    std::size_t j=0;
    switch (gp.sptype) {
        case GSL_SPMATRIX_COO:
            for(std::size_t k=0;k<gp.nz;k++) sp.set_entry(gp.i[k],gp.p[k],gp.data[k]);
            break;
        case GSL_SPMATRIX_CSC:
            for(std::size_t k=0;k<gp.nz;k++){
                while(gp.p[j+1]<=k) j++;
                sp.set_entry(gp.i[k],j,gp.data[k]);
            }
            break;
        case GSL_SPMATRIX_CSR:
            for(std::size_t k=0;k<gp.nz;k++){
                while(gp.p[j+1]<=k) j++;
                sp.set_entry(j,gp.i[k],gp.data[k]);
            }
            break;
    }
    sp.clean();
}
#else
//long double
void sp2gpm(const sparse&sp,gsl_spmatrix_long_double&gp){
    gp=*gsl_spmatrix_long_double_alloc_nzmax(sp.m(),sp.n(),sp.size(),GSL_SPMATRIX_COO);
    for(auto it=sp.cbegin();it!=sp.cend();it++)
        gsl_spmatrix_long_double_set(&gp,it->first[0],it->first[1],it->second);
}
void gp2spm(const gsl_spmatrix_long_double&gp,sparse&sp){
    sp.zero();
    sp.resize(gp.size1,gp.size2);
    std::size_t j=0;
    switch (gp.sptype) {
        case GSL_SPMATRIX_COO:
            for(std::size_t k=0;k<gp.nz;k++) sp.set_entry(gp.i[k],gp.p[k],gp.data[k]);
            break;
        case GSL_SPMATRIX_CSC:
            for(std::size_t k=0;k<gp.nz;k++){
                while(gp.p[j+1]<=k) j++;
                sp.set_entry(gp.i[k],j,gp.data[k]);
            }
            break;
        case GSL_SPMATRIX_CSR:
            for(std::size_t k=0;k<gp.nz;k++){
                while(gp.p[j+1]<=k) j++;
                sp.set_entry(j,gp.i[k],gp.data[k]);
            }
            break;
    }
    sp.clean();
}
#endif //_long_JLG
#endif //SPARSE_JLG
