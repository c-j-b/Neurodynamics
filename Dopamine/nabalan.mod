TITLE Sodium ion accumulation without diffusion

NEURON {
	SUFFIX nabalan
	USEION na READ ina WRITE nai
	RANGE nainit,f
}

UNITS {
	(mM) = (milli/liter)
	(um) = (micron)
	FARADAY = (faraday) (coulomb)
	PI = (pi) (1)
}

PARAMETER {
        nainit = 3.54 (mM) 
        f = 1.00
}

ASSIGNED {
	ina (milliamp/cm2)
	diam (um)
}

STATE {
	nai (mM) <1e-4>
}

BREAKPOINT {
	SOLVE state METHOD sparse
}

INITIAL{
        nai=nainit
}

KINETIC state {
	COMPARTMENT PI*diam*diam/4 {nai}
	~ nai << (-f*ina*PI*diam*(1e4)/(FARADAY))
}

