TITLE Calcium ion accumulation without diffusion and buffering

NEURON {
	SUFFIX cabalan
	USEION ca READ cai, ica WRITE cai
	RANGE  cainit, fCa
}

UNITS {
	(molar) = (1/liter)
	(mM) =  (millimolar)
	(um) =  (micron)
	(mA) =  (milliamp)
	FARADAY = (faraday) (coulomb) 
	PI = (pi) (1)
}

PARAMETER {
         fCa = 0.005  (1)
         cainit = 0.000250 (mM)
}

ASSIGNED {
	diam  (um)
	ica   (mA/cm2)
		
}

STATE {
	cai (mM) <1e-10>
}

BREAKPOINT {
	SOLVE state METHOD sparse
}

INITIAL{
        cai=cainit
}

KINETIC state {
	COMPARTMENT PI*diam*diam/4 {cai}
	 ~ cai << (-fCa*ica*PI*diam*(1e4)/(2*FARADAY))
}
