TITLE  squid sodium, potassium delayed rectifier, and potassium A channels
 
UNITS {
        (molar) = (1/liter)
        (S) = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
         F = (faraday) (coulomb)
         R = (mole k) (mV-coulomb/degC)
        (mM) =  (millimolar)
}
 
NEURON {
        SUFFIX hh3
        USEION na READ nai WRITE ina
        USEION k WRITE ik
        RANGE  gnabar,gkhhbar,gkabar,ina,ikhh,ika,ik,ena,miv,hiv,htv1,htv2
        GLOBAL minf,hinf,ninf
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v   (mV)
        dt  (ms)
	nai (mM)
	celsius = 35.0 (degC)
        gnabar  = 5500.0e-6 (S/cm2)
        gkhhbar = 1000.0e-6 (S/cm2)
        gkabar  = 100.0e-6  (S/cm2)
        miv  = 44.6 (mV)
	hiv  = 66.8 (mV)
 	htv1 = 39.0 (mV)
 	htv2 = 59.0 (mV)
        ek  = -100  (mV)
        nao =  145  (mM)
        
 	
}
 
STATE {
        m <1e-4> h <1e-4> n <1e-4> p <1e-4> q <1e-4>
}
 
ASSIGNED {
        ina (mA/cm2)
        ik (mA/cm2)
        ika (mA/cm2)
        ikhh (mA/cm2)
        ena (mV)
        minf hinf ninf qinf pinf
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        ena = R*(celsius+273.15)/F*log(nao/nai)
        ina = gnabar*m*m*m*h*(v - ena)
        ikhh = gkhhbar*n*n*n*(v - ek)      
        ika = gkabar*p*p*p*q*(v - ek)      
        ik = ika + ikhh
}
 
UNITSOFF
 
INITIAL {
        m = boltz(v,-miv,6.0)
        h = boltz(v,-hiv,-7.8)
        n = boltz(v,-35,12.0)
        p = boltz(v,-42,4.0)
        q = boltz(v,-63,-4.0)
}

DERIVATIVE states {  :Computes state variables m, h, and n 
LOCAL minf,hinf,ninf,pinf,qinf,mtau,htau,ntau,ptau,qtau
        minf = boltz(v,-miv,6.0)
        hinf = boltz(v,-hiv,-7.8)
        ninf = boltz(v,-35,12.0)
        pinf = boltz(v,-42,4.0)
        qinf = boltz(v,-63,-4.0)
        mtau = boltz(v,-45.0,-1.5) - boltz(v,-65.0,-0.5) +0.04
        htau = 56.0*boltz(v,-htv1,-4.5) - 56.0*boltz(v,-htv2,-2.0) +1.0
        ntau = 10.0
        ptau = 5.5*exp(-(v+42)*(v+42)/100)+4.0
        qtau = 50.0
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
        p' = (pinf-p)/ptau
        q' = (qinf-q)/qtau
}
 
 
 
FUNCTION boltz(x,y,z) {
                boltz = 1/(1 + exp(-(x - y)/z))
}
 
UNITSON

