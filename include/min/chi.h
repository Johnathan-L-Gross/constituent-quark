//
//  chi.h
//
//
//  Created by Johnathan Gross
//
//

#ifndef ____chi_JLG_
#define ____chi_JLG_

#include "type.h"
#include <functional>

typedef std::function<void(const Dvector&,Dvector&,Dvector&)> cfunc3;
typedef std::function<void(const Dvector&,Dvector&)> cfunc2;

Double chisquared(const Dvector&,//model
                  const Dvector&,//data
                  const Dvector&);//errors
//Calculates the chi^2 of model to a set of data values with errors
Double chisquared(const Dvector&,//model
                  const Dvector&,//data
                  const Double);//error
//overload with only one error

Double deltaC(const Dvector&,//out
              const Dvector&,//data
              const Dvector&);//err;
/*
 Calculates the offset to be added to out that minimizes chi^2
 when calculated with data and err
 Result: value of offset
 */
void newC(Dvector&,//o
          Double&,//C
          const Dvector&,//d
          const Dvector&);//e
/*
 Adds deltaC(o,d,e) to values of o and C
 Effect: modifies o and C
 */

inline cfunc3 c2to3(const cfunc2&f){
    return [&f](const Dvector&a,Dvector&b[[maybe_unused]],Dvector&c){f(a,c);};
}
inline cfunc2 c3to2(const cfunc3&f){
    return [&f](const Dvector&a,Dvector&c){
        Dvector b{};
        f(a,b,c);
    };
}

class chimin {
private:
    Dvector data,//data to be fit to
    sig,//errors to weigh fit
    eparam,//set of parameters to fit with
    iparam,//set of parameters stored between iterations of func
    vals;//output values to fit against data
    std::size_t nd,//size of data
    ne,//size of eparam
    neo,//size of parameters not constant offset
    ni,//size of iparam
    df;//degrees of freedom
    Double oerror,//chi^2 precision
    ierror,//internal precision
    chisq;//value of chi^2
    bool off,//last parameter is a constant offset
    rel;//relative or absolute precision
    cfunc3 func;//output
    /*function that uses fixed external parameters to produce output values
     Effect:overwrites internal parameters and output
     */

     chimin()=delete;//no default constructor

public:
    chimin(const cfunc3&,//func
           const Dvector&,//data
           const Dvector&,//sig
           const Dvector&,//eparam
           const Dvector&,//iparam
           const Dvector&,//vals
           const Double=0,//oerror
           const Double=0,//ierror
           const bool=false,//off
           const bool=false);//rel
    /*
     Constructor
     */
    chimin(const cfunc3&f,
           const Dvector&d,
           const Dvector&s,
           const Dvector&e,
           const Dvector&i,
           const Double oer=0,
           const Double ier=0,
           const bool of=false,
           const bool r=false):
    chimin(f,d,s,e,i,{},oer,ier,of,r){}
    //overload with no values given
    
    chimin(const cfunc2&,//func
           const Dvector&,//data
           const Dvector&,//sig
           const Dvector&,//eparam
           const Dvector&,//iparam
           const Dvector&,//vals
           const Double=0,//oerror
           const Double=0,//ierror
           const bool=false,//off
           const bool=false);//rel
    //overload for cfunc2
    chimin(const cfunc2&f,
           const Dvector&d,
           const Dvector&s,
           const Dvector&e,
           const Dvector&v,
           const Double oer=0,
           const Double ier=0,
           const bool of=false,
           const bool r=false):
    chimin(f,d,s,e,{},v,oer,ier,of,r){}
    //overload for cfunc2 with no iparam given
    chimin(const cfunc2&f,
           const Dvector&d,
           const Dvector&s,
           const Dvector&e,
           const Double oer=0,
           const Double ier=0,
           const bool of=false,
           const bool r=false):
    chimin(f,d,s,e,{},{},oer,ier,of,r){}
    //overload for cfunc2 with no iparam or values given

    Double eval(const Dvector&,//e
                Dvector&,//i
                Dvector&);//o
    /*
     Evaluates func at e
     Effect:modifies i with ne internal parameters,returns outputs in o
     Result:chi^2 value of o
     */
    inline Double eval(const Dvector&e,Dvector&o){return eval(e,iparam,o);}
    //overload modifying iparam
    inline Double eval(Dvector&o){return eval(eparam,iparam,o);}
    //overload using eparam and modifying iparam
    inline Double eval(){
    //overload that uses eparam, modifies iparam, and stores chi^2 and output
        chisq=eval(eparam,iparam,vals);
        return chisq;
    }

    void minimize(Dvector&,//initial guess
                  Dmatrix&,//initial directions
                  bool=false);//display flag
    /*
     Routine to minimize chi^2 with respect to parameters guess,
     along directions
     Effect:modifies eparam, iparam,vals, and chisq
     */
    inline void minimize(Dmatrix&mat,bool dis=false){minimize(eparam,mat,dis);}
    //overload that uses eparam as guess
    inline void minimize(Dvector&p,bool dis=false){
        Dmatrix M(0);
        minimize(p,M,dis);
    }
    inline void minimize(bool dis=false){
        Dmatrix M(0);
        minimize(eparam,M,dis);
    }

    //inline void sd(const Dvector&x){data=x;}
    //inline void ss(const Dvector&x){sig=x;}
    //inline void se(const Dvector&x){eparam=x;}
    //inline void si(const Dvector&x){iparam=x;}
    inline const Dvector& gd()const{return data;}
    inline const Dvector& gs()const{return sig;}
    inline const Dvector& ge()const{return eparam;}
    inline const Dvector& gi()const{return iparam;}
    inline const Dvector& gv()const{return vals;}
    inline Double gd(std::size_t i)const{return data[i];}
    inline Double gs(std::size_t i)const{return sig[i];}
    inline Double ge(std::size_t i)const{return eparam[i];}
    inline Double gi(std::size_t i)const{return iparam[i];}
    inline Double gv(std::size_t i)const{return vals[i];}
    inline std::size_t gnd()const{return nd;}
    inline std::size_t gne()const{return ne;}
    inline std::size_t gneo()const{return neo;}
    inline std::size_t gni()const{return ni;}
    inline std::size_t gdf()const{return df;}
    inline Double gierror()const{return ierror;}
    inline Double goerror()const{return oerror;}
    inline Double gchi()const{return chisq;}
    inline bool grel()const{return rel;}
    inline bool goff()const{return off;}
    inline const cfunc3& gfunc()const{return func;}
};
#endif //____chi_JLG_
