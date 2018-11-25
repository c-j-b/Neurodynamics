TITLE Calcium pump
 

UNITS {
        (molar) = (1/liter)
        (pA) =  (picoamp)
	(mV) =	(millivolt)
        (uS) =  (micromho)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


INDEPENDENT {v FROM -100 TO 50 WITH 50 (mV)}

NEURON {
	SUFFIX capump
	USEION ca READ cai WRITE ica
	RANGE  icapump,icapumpmax,km
 
}


PARAMETER {
        dt    (ms)
        cai   (mM)
        celsius = 35  (degC)
        icapumpmax  = 0.0312  (mA/cm2)
        km = 0.000500         (mM)
                 
}

ASSIGNED { 
          ica     (mA/cm2)
          icapump (mA/cm2)
}


BREAKPOINT {
        icapump = icapumpmax*(1/(1 + km/cai))
	ica = icapump
}
