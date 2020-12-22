//
//  quark.cpp
//  
//
//  Created by Johnathan Gross
//

#include "quark.h"
#include "functions.h"
#include "baryon/bar_calc.h"

hfunc strong(const std::vector<std::pair<Double,Double>>&a,
             const Double b,
             const Double sigma0,
             const Double s,
             const Double eCoul,
             const Double eCont,
             const Double eTen,
             const wf2&wf,
             const int3&par,
             const Double error){
    const bool eeq=(eCont==eTen),
    calcCoul=(-0.5!=eCoul),
    calcCont=(-0.5!=eCont),
    mneq=(wf.gmass(0)!=wf.gmass(1));
    const bool calcTen=(!eeq and -0.5!=eTen);
    const std::size_t s1=wf.size1(),
    s2=wf.size2();
    const Double s12=sigmaij2(sigma0,s,wf.gmass()),
    s23=sigmaij2(sigma0,s,mrot(1,wf.gmass())),
    s31=sigmaij2(sigma0,s,mrot(2,wf.gmass()));
    return[=,&a,&wf,&par](sparse&h,const Double omega){
        const std::size_t s2b=2*(int(calcCoul)+int(calcCont)+int(calcTen)),
        s3b=int(calcCoul)+int(calcCont)+int(calcTen);
        //std::cerr<<eeq<<calcCoul<<calcCont<<calcTen<<mneq<<' '<<s2b<<' '<<s3b<<std::endl;
        std::vector<sparse> ham0{sparse(s1,s1,error/10)},
        ham1{sparse(s1,s1,error/10)},
        ham2(7,sparse(s2,s2,error/10)),
        ham2b(s2b,sparse(s2,s2,error/10)),
        ham3(4,sparse(s2,s2,error/10)),
        ham3b(s3b,sparse(s2,s2,error/10));
        
#pragma omp parallel sections
        {
            //3 kinetic energy, 1D integral
#pragma omp section
            hamcalc(ham0,wf.states1(),par,omega,wf.gmass(),
                    //functions
                    {relkin},
                    //k
                    {0b01},
                    //tens
                    {0},
                    //spin
                    {3},
                    //rot
                    {3},
                    true,error/10);
            //3body confinement, 3D integration
#pragma omp section
            hamcalc(ham1,wf.states1(),par,omega,wf.gmass(),
                    //functions
                    {conf(b)},
                    //k
                    {0b00},
                    //tens
                    {0},
                    //spin
                    {3},
                    //rot
                    {3},
                    true,error/10);
            //12 suppressed terms and 23 terms
#pragma omp section
            {
                std::vector<func1>
                f2{relkin,//0 23
                    fullcoulomb(a,s12),//1 12
                    fullcoulomb(a,s23),//2 23
                    fullcontact(a,s12),//3 12
                    fullcontact(a,s23),//4 23
                    fulltensor(a,s12),//5 12
                    fulltensor(a,s23)};//6 23
                ivector
                k2(7,0b10),
                t2{0,0,0,1,1,2,2},
                sp2(7,12),
                r2{1,3,1,3,1,3,1};
                k2[0]=0b01;
                
                hamcalc(ham2,wf.states2(),par,omega,wf.gmass(),
                        f2,k2,t2,sp2,r2,true);
            }
            //12 and 23 betas and deltas
#pragma omp section
            {
                std::vector<func1> f2b(s2b);
                ivector k2b(s2b,0b11),
                t2b(s2b,0),
                sp2b(s2b,12),
                r2b(s2b,3);
                
                if(calcCoul){
                    f2b[0]=beta(eCoul);
                    f2b[1]=beta(eCoul);
                    r2b[1]=1;
                }
                if(calcCont){
                    f2b[2*int(calcCoul)]=delta(eCont);
                    f2b[2*int(calcCoul)+1]=delta(eCont);
                    r2b[2*int(calcCoul)+1]=1;
                }
                if(calcTen){
                    f2b[s2b-2]=delta(eTen);
                    f2b[s2b-1]=delta(eTen);
                    r2b[s2b-1]=1;
                }
                
                hamcalc(ham2b,wf.states2(),par,omega,wf.gmass(),
                        f2b,k2b,t2b,sp2b,r2b,true);
            }
            //31 terms
#pragma omp section
            if(mneq){
                std::vector<func1>
                f3{relkin,//0
                    fullcoulomb(a,s31),//1
                    fullcontact(a,s31),//2
                    fulltensor(a,s31)};//3
                ivector k3(4,0b10),
                t3{0,0,1,2},
                sp3(4,12),
                r3(4,2);
                k3[0]=0b01;
                
                hamcalc(ham3,wf.states2(),par,omega,wf.gmass(),
                        f3,k3,t3,sp3,r3,true);
            }
            //31 betas and deltas
#pragma omp section
            if(mneq){
                std::vector<func1> f3b(s3b);
                ivector k3b(s3b,0b11),
                t3b(s3b,0),
                sp3b(s3b,12),
                r3b(s3b,2);
                if(calcCoul) f3b[0]=beta(eCoul);
                if(calcCont) f3b[calcCoul]=delta(eCont);
                if(calcTen) f3b[s3b-1]=delta(eTen);
                    
                hamcalc(ham3b,wf.states2(),par,omega,wf.gmass(),
                        f3b,k3b,t3b,sp3b,r3b,true);
            }
        }//sections
        
        sparse h12(s2,s2,error/10),
        h23(s2,s2,error/10),
        h31(s2,s2,error/10);
        
        //contact and tensor terms
        if(eeq and calcCont){
            //deltaC*(Cont+Ten)*deltaC
            h12=ham2b[2*calcCoul]*(ham2[3]+ham2[5])*ham2b[2*calcCoul];
            h23=ham2b[2*calcCoul+1]*(ham2[4]+ham2[6])*ham2b[2*calcCoul+1];
            if(mneq) h31=ham3b[calcCoul]*(ham3[2]+ham3[3])*ham3b[calcCoul];
        }else if(calcCont and calcTen){
            //deltaC*Cont*deltaC+deltaT*Ten*deltaT
            h12=ham2b[2*calcCoul]*ham2[3]*ham2b[2*calcCoul]
            +ham2b[s2b-2]*ham2[5]*ham2b[s2b-2];
            h23=ham2b[2*calcCoul+1]*ham2[4]*ham2b[2*calcCoul+1]
            +ham2b[s2b-1]*ham2[6]*ham2b[s2b-1];
            if(mneq)
                h31=ham3[calcCoul]*ham3[2]*ham3b[calcCoul]
                +ham3b[s3b-1]*ham3[3]*ham3b[int(calcCoul)+1];
        }else if(calcCont){
            //deltaC*Cont*deltaC+Ten
            h12=ham2b[2*calcCoul]*ham2[3]*ham2b[2*calcCoul]+ham2[5];
            h23=ham2b[2*calcCoul+1]*ham2[4]*ham2b[2*calcCoul+1]+ham2[6];
            if(mneq) h31=ham3b[calcCoul]*ham3[2]*ham3b[calcCoul]+ham3[3];
        }else if(calcTen){
            //Cont+deltaT*Ten*deltaT
            h12=ham2[3]+ham2b[s2b-2]*ham2[5]*ham2b[s2b-2];
            h23=ham2[4]+ham2b[s2b-1]*ham2[6]*ham2b[s2b-1];
            if(mneq) h31=ham3[2]+ham3b[s3b-1]*ham3[3]*ham3b[s3b-1];
        }else{
            //Cont+Ten
            h12=ham2[3]+ham2[5];
            h23=ham2[4]+ham2[6];
            if(mneq) h31=ham3[2]+ham3[3];
        }

        sparse tmat=transmutation(wf.gtrans(),s2);
        sparse mosh1(wf.w2().gmosh1()),
        tmosh2(tmat*wf.w2().gmosh2());
                
        //3 kinetic energy
        h=ham0[0]
        //3body confinement
        +ham1[0]
        +tmat*(
               //12 Coulomb
               (calcCoul?ham2b[0]*ham2[1]*ham2b[0]:ham2[1])
               //12 contact and 12 tensor
               +h12
               +(mneq?1:2)
               *mosh1*(
                       //1 kinetic energy
                       ham2[0]
                       //23 Coulomb
                       +(calcCoul?ham2b[1]*ham2[2]*ham2b[1]:ham2[2])
                       //23 contact and 23 tensor
                       +h23
                       )/mosh1
               )/tmat;
        if(mneq) h+=
            tmosh2*(
                    //2 kinetic energy
                    ham3[0]
                    //31 Coulomb
                    +(calcCoul?ham3b[0]*ham3[1]*ham3b[0]:ham3[1])
                    //31 contact and 31 tensor
                    +h31
                    )/tmosh2;
        h=(h+h.transpose())/2;
        h.set_error(error,true);
    };
}

hfunc electric(const Double aEM,
               const Double gEM,
               const marray&ch,
               const wf2&wf,
               const int1&par,
               const Double error){
    const bool mneq=(wf.gmass(0)!=wf.gmass(1) or ch[0]!=ch[1]);
    const std::size_t s1=wf.size1(),
    s2=wf.size2();
    return [=,&ch,&wf,&par](sparse&h,const Double omega){
        std::vector<sparse> ham1{sparse(s1,s1,error/10)},
        ham2{sparse(s2,s2,error/10)},
        ham3{sparse(s2,s2,error/10)};
        
#pragma omp parallel sections
        {
#pragma omp section
            hamcalc(ham1,wf.states1(),par,omega,wf.gmass(),
                    {EMsmearcoulomb(aEM,gEM,ch)},
                    {0b01},
                    {0},
                    {12},
                    {3},
                    true);
#pragma omp section
            hamcalc(ham2,wf.states2(),par,omega,wf.gmass(),
                    {EMsmearcoulomb(aEM,gEM,mrot(1,ch))},
                    {0b10},
                    {0},
                    {12},
                    {1},
                    true);
#pragma omp section
            if(mneq){
                hamcalc(ham3,wf.states2(),par,omega,wf.gmass(),
                        {EMsmearcoulomb(aEM,gEM,mrot(2,ch))},
                        {0b10},
                        {0},
                        {12},
                        {2},
                        true);
            }
        }
        
        sparse tmat=transmutation(wf.gtrans(),s2),
        tmosh1=tmat*wf.w2().gmosh1();
        sparse tmosh2=(mneq?tmat*wf.w2().gmosh2():sparse());
        h=ham1[0]+(mneq?1:2)*tmosh1*ham2[0]/tmosh1;
        if(mneq) h+=tmosh2*ham3[0]/tmosh2;
        h=(h+h.transpose())/2;
        h.set_error(error,true);
    };
}

hfunc electromagnetic(const Double aEM,
                      const Double gEM,
                      const Double eEM,
                      const marray&ch,
                      const wf2&wf,
                      const int1&par,
                      const Double error){
    const bool mneq=(wf.gmass(0)!=wf.gmass(1) or ch[0]!=ch[1]),
    calcEM=(-0.5!=eEM);
    const std::size_t s1=wf.size1(),
    s2=wf.size2();
    return [=,&ch,&wf,&par](sparse&h,const Double omega){
        std::vector<sparse> ham1{sparse(s1,s1,error/10)},
        ham2(3,sparse(s2,s2,error/10)),
        ham2b(2,sparse(s2,s2,error/10)),
        ham3(2,sparse(s2,s2,error/10)),
        ham3b(1,sparse(s2,s2,error/10));
        
#pragma omp parallel sections
        {
#pragma omp section
            hamcalc(ham1,wf.states1(),par,omega,wf.gmass(),
                    {EMsmearcoulomb(aEM,gEM,ch)},
                    {0b01},
                    {0},
                    {12},
                    {3},
                    true);
#pragma omp section
            hamcalc(ham2,wf.states2(),par,omega,wf.gmass(),
                    {EMsmearcoulomb(aEM,gEM,mrot(1,ch)),//0
                    EMsmearcontact(aEM,gEM,ch),//1
                    EMsmearcontact(aEM,gEM,mrot(1,ch))},//2
                    {0b10,0b10,0b10},
                    {0,1,1},
                    {12,12,12},
                    {1,3,1},
                    true);
#pragma omp section
            if(calcEM){
                hamcalc(ham2b,wf.states2(),par,omega,wf.gmass(),
                        {delta(eEM),delta(eEM)},
                        {0b11,0b11},
                        {0,0},
                        {12,12},
                        {3,1},
                        true);
            }
#pragma omp section
            if(mneq){
                hamcalc(ham3,wf.states2(),par,omega,wf.gmass(),
                        {EMsmearcoulomb(aEM,gEM,mrot(2,ch)),
                        EMsmearcontact(aEM,gEM,mrot(2,ch))},
                        {0b10,0b10},{0,1},{12,12},{2,2},true);
            }
#pragma omp section
            if(mneq and calcEM)
                hamcalc(ham3b,wf.states2(),par,omega,wf.gmass(),
                        {delta(eEM)},{0b11},{0},{12},{2},true);
        }
        
        sparse tmat=transmutation(wf.gtrans(),s2);
        sparse mosh1(wf.w2().gmosh1()),
        tmosh2=(mneq?tmat*wf.w2().gmosh2():sparse());
        h=
        ham1[0]
        //12 Coulomb
        +tmat*(
               (calcEM?ham2b[0]*ham2[1]*ham2b[0]:ham2[1])
               //12 contact
               +(mneq?1:2)
               *mosh1*(
                       ham2[0]
                       //23 Coulomb
                       +(calcEM?ham2b[1]*ham2[2]*ham2[1]:ham2[2])
                       //23 contact
                       )/mosh1
               )/tmat;
        if(mneq) h+=
            tmosh2*(
                    ham3[0]
                    //31 Coulomb
                    +(calcEM?ham3[0]*ham3[1]*ham3[0]:ham3[1])
                    //31 contact
                    )/tmosh2;
        h=(h+h.transpose())/2;
        h.set_error(error,true);
    };
}
