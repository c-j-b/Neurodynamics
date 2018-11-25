TITLE NMDA receptor as a distributed mechanism as in Canavier and Landry, 2006

COMMENT
Landry did not multiply the calcium current by z squared,
hence the effective permeability ratio is not 10.6 but 2.65.
(To make this file equivalent to Komendantov's nmda.mod 
(see Komendantov et al. 2004 ModelDB entry),
the 91st line here should be ica = 4*10.6*power*numerca/denom2).
Also Landry used the somatic calcium concentration
to drive the dendritic calcium component of the NMDA current, while
in Komendantov et al. 2004, a constant calcium concentration cai
is used in the dendrites instead.
Other than there two discrepancies, the description of the current
is identical despite the different formulations.
ENDCOMMENT

UNITS {
        (pA) = (picoamp)
        (molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX nmda
	USEION ca WRITE ica
	USEION na READ nai WRITE ina
	USEION k WRITE ik
	RANGE  ica,ina,ik,inmda,Pbar,mg,km,nai,pr
        GLOBAL pinf
        POINTER caisoma  :will take calcium concentration from the soma
       
}

UNITS {
	:FARADAY = 96520 (coul)
	:R = 8.3134 (joule/degC)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
        v (mV)
	celsius= 35  	(degC)
	Pbar = 0.0	(cm/s)	: Maximum Permeability PNMDA in Laundry
	cao = 2.0	(mM)
	lamdaca = 0.3 
        lamda = 0.75
        pr = 0.0225
	nao = 145	(mM)
	ki =  140	(mM)
	ko = 2.5	(mM)
        dt (ms)
        q=9 (mV)
        km=50.7 (mM)
        mg = 1.2 (mM)
}

STATE {
         p <1e-4>
}

ASSIGNED { 
	   nai          (mM)
           ica		(mA/cm2)
           ina		(mA/cm2)
           ik		(mA/cm2)
           inmda 	(mA/cm2)
           pinf
           caisoma      (mM)
}

LOCAL arg, power, numerna, numerk, numerca, denom, denom2

BREAKPOINT {
     SOLVE states METHOD cnexp
	pinf = pr + (1.0 - pr)/(1 + (mg/km)*exp(-v/q))
        arg = v*FARADAY /((1000)*R*(celsius+273.15))
	power =  Pbar*(0.000001)*v*p*FARADAY*FARADAY/(R*(celsius+273.15))
	numerna = lamda*nai - lamda*nao*exp(-arg)
	numerk = lamda*ki - lamda*ko*exp(-arg)
	numerca = caisoma - lamdaca*cao*exp(-2*arg)
	denom = 1 - exp(-arg)
	denom2 = 1 - exp(-2*arg)
	
        ina = power*numerna/denom
	ik = power*numerk/denom
	ica = 10.6*power*numerca/denom2
	inmda = ina + ik + ica
}
UNITSOFF
DERIVATIVE states {
	pinf = pr + (1.0 - pr)/(1 + (mg/km)*exp(-v/q))
	p' = pinf - p 
	}
UNITSON
