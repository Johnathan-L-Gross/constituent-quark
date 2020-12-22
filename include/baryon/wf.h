//
//  wf.h
//  
//
//  Created by Johnathan Gross
//
//

#ifndef ____wf__
#define ____wf__

#include "sparse/sparse.h"
#include <iostream>

typedef std::array<Double,3> marray;

class wf1;
class wf2;

class wf1{//just one set of states
    friend class wf2;
    //friend void swap(wf1&,wf1&);
protected:
    imatrix states;
    int Pf,J,P,M,Nmin,Nmax;
    marray mass;
    bool sym;
    size_t size;
    sparse mosh1,mosh2;//231 and 312

    //funcitons to change the state of the system
    inline void sstates(imatrix&s){states=s;}
    inline void smosh1(sparse&m){mosh1=m;}
    inline void smosh2(sparse&m){mosh2=m;}
public:
    wf1(const imatrix&,//states
        const marray&);//mass
    //given statelist and masses
    wf1(const int,//Pf
        const int,//J
        const int,//P
        const int,//M
        const int,//Nmin
        const int,//Nmax 6
        const marray&,//mass
        const bool=true);//sym
    //given all q# and masses
    wf1(const int Pf_,
        const int J_,
        const int P_,
        const int M_,
        const int Nmax_,//5
        const marray&mass_,
        const bool sym_=true):
    wf1(Pf_,J_,P_,M_,0,Nmax_,mass_,sym_){}
    //no min
    wf1(const int Pf_,
        const int J_,
        const int P_,
        const int Nmax_,//4
        const marray&mass_,
        const bool sym_=true):
    wf1(Pf_,J_,P_,J_,Nmax_,mass_,sym_){}
    //no M
    
    wf1()=default;//default constructor
    //wf1(const wf1&)=default;//copy constructor
    //wf1(wf1&&)=default;//move constructor
    //~wf1()=default;//destructor
    //wf1& operator=(wf1);//copy assignment

    /*wf1(const imatrix&,//states
        const marray&,//mass
        const Dmatrix&,//mosh1
        const Dmatrix&);//mosh2*/
    //constructor not to be used by anything but wf2

    inline const ivector& operator[](size_t i)const{
        return states[i];
    }
    inline const ivector& at(size_t i)const{return states.at(i);}

    //return const references
    inline const imatrix&gstates()const{return states;}
    inline const sparse&gmosh1()const{return mosh1;}
    inline const sparse&gmosh2()const{return mosh2;}
    inline const marray&gmass()const{return mass;}
    inline const Double&gmass(std::size_t i)const{return mass[i];}

    //return values
    inline size_t gsize() const{return size;}
    inline int gNmin() const{return Nmin;}
    inline int gNmax() const{return Nmax;}
    inline bool gsym() const{return sym;}

    //return state values
#ifdef TESTING_STATES_JLG
    inline int gPf(size_t n=0)const{return states.at(n)[0];}
    inline int gJ(size_t n=0)const{return states.at(n)[1];}
    inline int gP(size_t n=0)const{return states.at(n)[2];}
    inline int gM(size_t n=0)const{return states.at(n)[3];}
    inline int gN()const{return states.at(size-1)[4];}
    inline int gN(size_t n)const{return states.at(n)[4];}
    inline int gL(size_t n)const{return states.at(n)[5];}
    inline int gnp(size_t n)const{return states.at(n)[6];}
    inline int glp(size_t n)const{return states.at(n)[7];}
    inline int gnl(size_t n)const{return states.at(n)[8];}
    inline int gll(size_t n)const{return states.at(n)[9];}
    inline int gS(size_t n)const{return states.at(n)[10];}
    inline int gPs(size_t n)const{return states.at(n)[11];}
#else
    inline int gPf(size_t n=0)const{return states[n][0];}
    inline int gJ(size_t n=0)const{return states[n][1];}
    inline int gP(size_t n=0)const{return states[n][2];}
    inline int gM(size_t n=0)const{return states[n][3];}
    inline int gN()const{return states[size-1][4];}
    inline int gN(size_t n)const{return states[n][4];}
    inline int gL(size_t n)const{return states[n][5];}
    inline int gnp(size_t n)const{return states[n][6];}
    inline int glp(size_t n)const{return states[n][7];}
    inline int gnl(size_t n)const{return states[n][8];}
    inline int gll(size_t n)const{return states[n][9];}
    inline int gS(size_t n)const{return states[n][10];}
    inline int gPs(size_t n)const{return states[n][11];}
#endif

    void display(const std::size_t s=-1,
                 std::ostream&o=std::cout)const;
    void display(std::ostream&o)const{display(-1,o);}

    //parity transformations
    sparse P12()const;
    sparse P12ss()const;
    sparse P13()const;
    sparse P13ss()const;
    sparse P23()const;
    sparse P23ss()const;
};

class wf2{//just two sets of states
    //B should be bigger than A
    //friend void swap(wf2&,wf2&);
protected:
    wf1 wfa,wfb;
    bool intern=false;
    tvector trans;
    sparse mosh1,mosh2;

public:
    wf2(const wf1&,//wfa
        const wf1&);//wfb
    wf2(const imatrix&s1,
        const imatrix&s2,
        const marray&mass):
    wf2(wf1(s1,mass),wf1(s2,mass)){}
    wf2(const int Pf1,
        const int J1,
        const int P1,
        const int M1,
        const int Nmin1,
        const int Nmax1,
        const int Pf2,
        const int J2,
        const int P2,
        const int M2,
        const int Nmin2,
        const int Nmax2,//12
        const marray&mass,
        const bool sym1=true,
        const bool sym2=true):
    wf2(wf1(Pf1,J1,P1,M1,Nmin1,Nmax1,mass,sym1),
        wf1(Pf2,J2,P2,M2,Nmin2,Nmax2,mass,sym2)){}
    wf2(const int Pf1,
        const int J1,
        const int P1,
        const int M1,
        const int Nmax1,
        const int Pf2,
        const int J2,
        const int P2,
        const int M2,
        const int Nmax2,//10
        const marray&mass,
        const bool sym1=true,
        const bool sym2=true):
    wf2(Pf1,J1,P1,M1,0,Nmax1,Pf2,J2,P2,M2,0,Nmax2,mass,sym1,sym2){}
    //no min
    wf2(const int Pf1,
        const int J1,
        const int P1,
        const int Nmax1,
        const int Pf2,
        const int J2,
        const int P2,
        const int Nmax2,//8
        const marray&mass,
        const bool sym1=true,
        const bool sym2=true):
    wf2(Pf1,J1,P1,J1,Nmax1,Pf2,J2,P2,J2,Nmax2,mass,sym1,sym2){}
    //no M
    wf2(const int Pf,
        const int J,
        const int P,
        const int M,
        const int Nmin,
        const int Nmax,//6
        const marray&mass,
        const bool sym1=true,
        const bool sym2=false):
    wf2(Pf,J,P,M,Nmin,Nmax,Pf,J,P,M,Nmin,Nmax,mass,sym1,sym2){}
    //only1
    wf2(const int Pf,
        const int J,
        const int P,
        const int M,
        const int Nmax,//5
        const marray&mass,
        const bool sym1=true,
        const bool sym2=false):
    wf2(Pf,J,P,M,0,Nmax,mass,sym1,sym2){}
    //no min, only 1
    wf2(const int Pf,
        const int J,
        const int P,
        const int Nmax,//4
        const marray&mass,
        const bool sym1=true,
        const bool sym2=false):
    wf2(Pf,J,P,J,Nmax,mass,sym1,sym2){}
    //no M, only 1
    
    wf2()=default;//default constructor
    //wf2(const wf2&);//copy constructor
    //wf2(wf2&&);//move constructor
    //~wf2();//destructor
    //wf2& operator=(wf2);//copy assignment

    inline const wf1&w1() const{return wfa;}
    inline const wf1&w2() const{return wfb;}
    inline const imatrix&states1() const{return wfa.gstates();}
    inline const imatrix&states2() const{return wfb.gstates();}
    inline const tvector&gtrans() const{return trans;}
    inline const sparse&gmosh1() const{return mosh1;}
    inline const sparse&gmosh2() const{return mosh2;}
    inline size_t size1() const{return wfa.size;}
    inline size_t size2() const{return wfb.size;}
    inline size_t tsize() const{return trans.size();}
    inline bool gintern() const{return intern;}
    inline const marray&gmass() const{return wfa.gmass();}
    inline const Double&gmass(std::size_t i)const{return wfa.gmass(i);}
    
    inline void display1(const std::size_t s=-1,
                         std::ostream&o=std::cout)const{wfa.display(s,o);}
    inline void display1(std::ostream&o)const{wfa.display(-1,o);}
    inline void display2(const std::size_t s=-1,
                         std::ostream&o=std::cout)const{wfb.display(s,o);}
    inline void display2(std::ostream&o)const{wfb.display(-1,o);}
};

#endif /* defined(____wf_b__) */
