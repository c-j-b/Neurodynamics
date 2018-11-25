TITLE Calcium dependent potassium channel (small conductance SK)
 

UNITS {
        (molar) = (1/liter)
        (pA) =  (picoamp)
	(mV) =	(millivolt)
        (S)  =  (siemens)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


INDEPENDENT {v FROM -100 TO 50 WITH 50 (mV)}

NEURON {
	SUFFIX kca
	USEION ca READ cai
	USEION k WRITE ik
	RANGE  gkbar,km,oinf,n
 
}


PARAMETER {
        dt  (ms)
        cai (mM)
        celsius = 35   (degC)
        gkbar = 800e-6 (S/cm2)
        ek = -100      (mV)
        km = 0.00019   (mM)
        n  = 4.0       (1)
        
        
}

ASSIGNED { 
           ik		(mA/cm2)
           oinf           
}


BREAKPOINT {
        oinf = 1/(1 + pow(km/cai,n))
	ik = oinf*gkbar*(v - ek)
}

