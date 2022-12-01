
from cProfile import label
import cmath
import math
import numpy as np
from tty import IFLAG

def coat(x,y,Refrel1,Refrel2):
    DEL=1.0e-8

    x1=Refrel1*x
    x2=Refrel2*x
    y2=Refrel2*y
    ystop=y+(4.0*y**0.3333)+2.0
    Refrel=Refrel2/Refrel1
    Nstop=ystop

    D0X1=cmath.acos(x1)/cmath.asin(x1)
    D0X2=cmath.acos(x2)/cmath.asin(x2)
    D0Y2=cmath.acos(y2)/cmath.asin(y2)

    PSIOY=math.cos(y)
    PSI1Y=math.sin(y) 
    CHIOY=-math.sin(y) 
    CHI1Y=math.cos(y) 
    XIOY=complex(PSIOY,-CHIOY) 
    XI1Y=complex(PSI1Y,-CHI1Y) 
    CHIOY2=-cmath.sin(y2) 
    CHI1Y2=cmath.cos(y2) 
    CHIOX2=-cmath.sin(x2) 
    CHI1X2=cmath.cos(x2) 
    QSCA=0.0
    QEXT=0.0
    XBACK=complex(0.0,0.0)

    N=1
    IFLAG=0
    RN=N
    PSIY=(2.*RN-1.)*PSI1Y/y-PSIOY 
    CHIY=(2.*RN-1.)*CHI1Y/y-CHIOY 
    XIY=complex(PSIY,-CHIY) 
    D1Y2=1./(RN/y2-D0Y2)-RN/y2 
    if (IFLAG==1):
        jump999()  
    D1X1=1.1(RN/x1-D0X1)-RN/x1 
    D1X2=1.1(RN/x2-D0X2)-RN/x2 
    CHIX2=(2.*RN-1.)*CHI1X2/x2-CHIOX2 
    CHIY2=(2.*RN-1.)*CHI1Y2/y2-CHIOY2 
    CHIPX2=CHI1X2-RN*CHIX2/x2 
    CHIPY2=CHI1Y2-RN*CHIY2/y2 
    ANCAP= Refrel*D1X1 -D1X2 
    ANCAP=ANCAP/(Refrel*D1X1 * CHIX2-CHIPX2)
    ANCAP=ANCAP/(CHIX2*D1X2-CHIPX2)
    BRACK=ANCAP*(CHIY2*D1Y2-CHIPY2)
    BNCAP=Refrel*D1X2-D1X1
    BNCAP=BNCAP/(Refrel*CHIPX2-D1X1*CHIX2)
    BNCAP=BNCAP/(CHIX2*D1X2-CHIPX2)
    CRACK=BNCAP*(CHIY2*D1Y2-CHIPY2)
    AMESS1=BRACK*CHIPY2
    AMESS2=BRACK*CHIPY2
    AMESS3=CRACK*CHIPY2
    AMESS4=CRACK*CHIY2
    if(abs(AMESS1)>DEL*abs(D1Y2)):
        jump999(Refrel2,y)
    if(abs(AMESS2)>DEL):
        jump999(Refrel2,y)
    if(abs(AMESS3)>DEL*abs(D1Y2)):
        jump999(Refrel2,y)
    if(abs(AMESS4)>DEL):
        jump999(Refrel2,y)

def jump999(RFREL2,y):
    DNBAR=coat.D1Y2-coat.BRACK*coat.CHIPY2
    DNBAR==DNBAR/(1.-coat.BRACK*coat.CHIY2)
    GNBAR=coat.D1Y2-coat.CRACK*coat.CHIPY2
    GNBAR=GNBAR/(1.-coat.CRACK*coat.CHIY2)
    AN=(DNBAR/RFREL2+coat.RN/y)*coat.HPSIY-coat.PSI1Y
    AN=AN/((DNBAR/RFREL2+coat.RN/y)*coat.XIY-coat.XI1Y)
    BN=(RFREL2*GNBAR+coat.RN/y)*coat.PSIY-coat.PSI1Y
    BN=BN/((RFREL2*GNBAR+coat.RN/y)*coat.XIY-coat.Xl1Y) 
    coat.QSCA=coat.QSCA+(2.*coat.RN+1.)*(abs(AN)*abs(AN)+abs(BN))
    coat.XBACK= coat.XBACK+ ( 2.*coat.RN + 1.)*(- 1.) **coat.N *(AN- BN)
    coat.QEXT=coat.QEXT+(2.*coat.RN+1.)*((AN).real+(BN).real) 
    coat.PSI0Y=coat.PSI1Y
    coat.PSI1Y=coat.PSIY
    coat.CHI0Y=coat.CHI1Y
    coat.CH11Y=coat.CHIY
    coat.XI1Y=complex(coat.PSI1Y,-coat.CHI1Y)
    coat.CHIOX2=coat.CHI1X2 
    coat.CHI1X2=coat.CHIX2 
    coat.CHIOY2=coat.CHI1Y2
    coat.CHI1Y2=coat.CHIY2 
    coat.DOX1=coat.D1X1 
    coat.DOX2=coat.D1X2
    coat.DOY2=coat.D1Y2
    coat.N=coat.N+1
    if (n-1-nstop lt 0) then goto,jump200
    if (coat.N-1-coat.Nstop>= 0):
        jump999
def jump300():
    coat.QSCA=(2./(coat.y*coat.y))*coat.QSCA
    coat.QEXT=(2./(coat.y*coat.y))*coat.QEXT
    coat.QBACK=coat.XBACK*np.conj(coat.XBACK)
    coat.QBACK= (1./(coat.y*coat.y))*coat.QBACK