//
//  algebra.cpp
//  
//
//  Created by Johnathan Gross
//
//

#include "algebra.h"
#include <cstdlib>
#ifdef TESTING_STATES_JLG
#include <stdexcept>
#include <iostream>
#endif
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_gamma.h>


#ifdef W_JLG_extern
map3 w3j_map_JLG;
map6 w6j_map_JLG;
map9 w9j_map_JLG;
map2 norm_map_JLG;
map2 norml_map_JLG;
#endif //W_JLG_extern

//******************************************************************************
//******************************************************************************
#ifdef W_JLG_map
double W3J(const unsigned int j1,
           const unsigned int j2,
           const unsigned int j3,
           const int m1,
           const int m2,
           const bool th){
    //careful, parameters are 2x what their value in mathematical definition
    //1/2->1, 5/2->5, 3->6

    //std::cerr<<"W3 "<<j1<<j2<<j3<<m1<<m2<<m3<<std::endl;

    if(!istriangle(j1,j2,j3)){
#ifdef TESTING_STATES_JLG
        if(th){
            std::cerr<<j1<<' '<<j2<<' '<<j3<<std::endl;
            throw std::invalid_argument("j's do not satisfy triangle inequality");
        }else return 0.0;
#else
        return 0.0;
#endif
    }
    if((mod(j1,2) != mod(m1,2)) or (mod(j2,2) != mod(m2,2))){
#ifdef TESTING_STATES_JLG
        if(th){
            std::cerr<<j1<<' '<<j2<<' '<<j3<<' '
            <<m1<<' '<<m2<<' '<<-m1-m2<<std::endl;
            throw std::invalid_argument("m and j not same half/full integer type");
        }else return 0.0;
#else
        return 0.0;
#endif
    }

    //symmetry relations
    //only calculates for j's in increasing order
    //then m's in increasing magnitude order
    //then only for m3<=0
    const int m3=-m1-m2;
    if(j1>j2) return W3J(j2,j1,j3,-m2,-m1);
    else if(j2>j3) return W3J(j1,j3,j2,-m1,-m3);
    else if(j1==j2 and std::abs(m1)>std::abs(m2)) return W3J(j2,j1,j3,-m2,-m1);
    else if(j2==j3 and std::abs(m2)>std::abs(m3)) return W3J(j1,j3,j2,-m1,-m3);
    else if(m3>0) return parity((j1+j2+j3)/2) * W3J(j1,j2,j3,-m1,-m2);

    const key5 key={int(j1),int(j2),int(j3),m1,m2};
#ifndef W_JLG_extern
    static map3 w3j_map_JLG;
#endif
#pragma omp critical(w3j)
    if(0==w3j_map_JLG.count(key)) w3j_map_JLG[key]=Wigner3(j1,j2,j3,m1,m2);
    
    return w3j_map_JLG[key];
}
#endif //W_JLG_map

//******************************************************************************
//******************************************************************************
double Wigner3(const unsigned int j1,
               const unsigned int j2,
               const unsigned int j3,
               const int m1,
               const int m2){
    //careful, parameters are 2x what their value in mathematical definition
    //1/2->1, 5/2->5, 3->6
    //std::cerr<<"W3 "<<j1<<j2<<j3<<m1<<m2<<m3<<std::endl;
    return gsl_sf_coupling_3j(j1,j2,j3,
                              m1,m2,-m1-m2);
}

//******************************************************************************
//******************************************************************************
#ifdef W_JLG_map
double W6J(const unsigned int j1,
           const unsigned int j2,
           const unsigned int j3,
           const unsigned int j4,
           const unsigned int j5,
           const unsigned int j6,
           const bool th){
    //careful, parameters are 2x what their value in mathematical definition
    //1/2->1, 5/2->5, 3->6
    //std::cerr<<"W6 "<<j1<<j2<<j3<<j4<<j5<<j6<<std::endl;

    if(!(istriangle(j1,j2,j3) and
         istriangle(j1,j5,j6) and
         istriangle(j4,j2,j6) and
         istriangle(j4,j5,j3))){
#ifdef TESTING_STATES_JLG
        if(th){
            std::cerr<<j1<<' '<<j2<<' '<<j3<<'\n'<<j4<<' '<<j5<<' '<<j6<<std::endl;
            throw std::invalid_argument("j's do not satisfy triangle inequalities");
        }else return 0.0;
#else
        return 0.0;
#endif
    }

    //symmetry relations
    //first by sum of column
    //then by first and second column
    if(j1+j4>j2+j5) return W6J(j2,j1,j3,
                               j5,j4,j6);
    else if(j2+j5>j3+j6) return W6J(j1,j3,j2,
                                    j4,j6,j5);
    else if(j1>j4){
        if(j2>j5) return W6J(j4,j5,j3,
                             j1,j2,j6);
        else return W6J(j4,j2,j6,
                        j1,j5,j3);
    } else if(j2>j5) return W6J(j1,j5,j6,
                                j4,j2,j3);

    const key6 key={j1,j2,j3,j4,j5,j6};
#ifndef W_JLG_extern
    static map6 w6j_map_JLG;
#endif
#pragma omp critical(w6j)
        if(0==w6j_map_JLG.count(key))
            w6j_map_JLG[key]=Wigner6(j1,j2,j3,
                                     j4,j5,j6);

    return w6j_map_JLG[key];
}
#endif //W_JLG_map

//******************************************************************************
//******************************************************************************
double Wigner6(const unsigned int j1,
               const unsigned int j2,
               const unsigned int j3,
               const unsigned int j4,
               const unsigned int j5,
               const unsigned int j6){
    //careful, parameters are double what their value in mathematical definition
    //1/2->1, 5/2->5, 3->6
    //std::cerr<<"W3 "<<j1<<j2<<j3<<m1<<m2<<m3<<std::endl;
    return gsl_sf_coupling_6j(j1,j2,j3,
                              j4,j5,j6);
}

//******************************************************************************
//******************************************************************************
#ifdef W_JLG_map
double W9J(const unsigned int j1,
           const unsigned int j2,
           const unsigned int j3,
           const unsigned int j4,
           const unsigned int j5,
           const unsigned int j6,
           const unsigned int j7,
           const unsigned int j8,
           const unsigned int j9,
           const bool th){
    //careful, parameters are 2x what their value in mathematical definition
    //1/2->1, 5/2->5, 3->6
    //std::cerr<<"W6 "<<j1<<j2<<j3<<j4<<j5<<j6<<j7<<j8<<j9<<std::endl;

    if(!(istriangle(j1,j2,j3) and
         istriangle(j4,j5,j6) and
         istriangle(j7,j8,j9) and
         istriangle(j1,j4,j7) and
         istriangle(j2,j5,j8) and
         istriangle(j3,j6,j9))){
#ifdef TESTING_STATES_JLG
        if(th){
            std::cerr<<j1<<' '<<j2<<' '<<j3<<'\n'
            <<j4<<' '<<j5<<' '<<j6<<'\n'
            <<j7<<' '<<j8<<' '<<j9<<std::endl;
            throw std::invalid_argument("j's do not satisfy triangle inequalities");
        }else return 0.0;
#else
        return 0.0;
#endif
    }

    //symmetry relations
    //first by sum of column
    //then by sum of row
    //c1<=r1 then c2<=r2 then c3<=r3
    const unsigned int c1=j1+j4+j7;
    const unsigned int c2=j2+j5+j8;
    const unsigned int c3=j3+j6+j9;
    if(c1<=c2 and c2<=c3);//do nothing, already permuted
    else if(c2<=c3 and c3<=c1) return W9J(j2,j3,j1,
                                          j5,j6,j4,
                                          j8,j9,j7);
    else if(c3<=c1 and c1<=c2) return W9J(j3,j1,j2,
                                          j6,j4,j5,
                                          j9,j7,j8);
    else if(c1<c3 and c3<c2) return W9J(j1,j3,j2,
                                        j7,j9,j8,
                                        j4,j6,j5);
    else if(c2<c1 and c1<c3) return W9J(j2,j1,j3,
                                        j8,j7,j9,
                                        j5,j4,j6);
    else return W9J(j3,j2,j1,
                    j9,j8,j7,
                    j6,j5,j4);

    const unsigned int r1=j1+j2+j3;
    const unsigned int r2=j4+j5+j6;
    const unsigned int r3=j7+j8+j9;
    const double eps=parity((c1+c2+c3)/2,double(1));
    if(r1<=r2 and r2<=r3);//do nothing, already permuted
    else if(r2<=r3 and r3<=r1) return W9J(j4,j5,j6,
                                          j7,j8,j9,
                                          j1,j2,j3);
    else if(r3<=r1 and r1<=r2) return W9J(j7,j8,j9,
                                          j1,j2,j3,
                                          j4,j5,j6);
    else if(r1<r3 and r3<r2) return eps*W9J(j1,j2,j3,
                                            j7,j8,j9,
                                            j4,j5,j6);
    else if(r2<r1 and r1<r3) return eps*W9J(j4,j5,j6,
                                            j1,j2,j3,
                                            j7,j8,j9);
    else return eps*W9J(j7,j8,j9,
                        j4,j5,j6,
                        j1,j2,j3);

    const key9 key={j1,j2,j3,j4,j5,j6,j7,j8,j9};
#ifndef W_JLG_extern
    static map9 w9j_map_JLG;
#endif
#pragma omp critical(w6j)
    if(0==w9j_map_JLG.count(key))
        w9j_map_JLG[key]=Wigner9(j1,j2,j3,
                                 j4,j5,j6,
                                 j7,j8,j9);

    return w9j_map_JLG[key];
}
#endif //W_JLG_map

//******************************************************************************
//******************************************************************************
double Wigner9(const unsigned int j1,
               const unsigned int j2,
               const unsigned int j3,
               const unsigned int j4,
               const unsigned int j5,
               const unsigned int j6,
               const unsigned int j7,
               const unsigned int j8,
               const unsigned int j9){
    //careful, parameters are double what their value in mathematical definition
    //1/2->1, 5/2->5, 3->6
    //std::cerr<<"W3 "<<j1<<j2<<j3<<m1<<m2<<m3<<std::endl;
    return gsl_sf_coupling_9j(j1,j2,j3,
                              j4,j5,j6,
                              j7,j8,j9);
}

//******************************************************************************
//******************************************************************************
#ifdef W_JLG_map
Double HN_l(const unsigned int n,
            const unsigned int l){
    //return (M_LN2+lgamma(n+1)-lgamma(n+l+1.5))/2;
    //logarithm of 3D HO normalization constant
    const key2 k={{n,l}};
#ifndef W_JLG_extern
    static map2 norml_map_JLG;
#endif
#pragma omp critical(hnl)
        if(0==norml_map_JLG.count(k))
            norml_map_JLG[k]=harm_norm_l(n,l);

    return norml_map_JLG[k];
}
#endif //W_JLG_map

//******************************************************************************
//******************************************************************************
Double harm_norm_l(const unsigned int n,
                   const unsigned int l){
    return (M_LN2-gsl_sf_lnfact(n)-gsl_sf_lngamma(n+l+1.5))/2.;
}
