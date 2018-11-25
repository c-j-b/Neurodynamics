TITLE AMPA receptor as distributed mechanism with random synaptic dynamics
 

UNITS {
        (pA) = (picoamp)
        (molar) = (1/liter)
	(mV) =	(millivolt)
        (S)  =  (siemens)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
         F = (faraday) (coulomb)
         R = (k-mole)  (joule/degC)
}


INDEPENDENT {v FROM -100 TO 50 WITH 50 (mV)}

NEURON {
	SUFFIX ampa
	USEION na READ nai,ena  WRITE ina
	USEION k  WRITE ik
	RANGE  iampa,ina,ik,gampa,gampak,nai,ratio
        POINTER ampasyn :take data from data file
 
}


PARAMETER {
        dt   (ms)
        ena  (mV)
        nai  (mM)
        celsius = 35  (degC)
        ratio = 1     (1) :ratio = 2 means doubled components
                          : of gampa, that is AMPA/NMDA ratio
                          : (peak AMPA EPSC/peak NMDA EPSC)
                          : is doubled	
        gampa  =  2.68e-6 (S/cm2)
        gampak =  3.37e-6 (S/cm2) 
        ek =  -100 (mV)
        nao = 145 (mM)
        
}

ASSIGNED { 
           ina	    (mA/cm2)
           ik       (mA/cm2)
           iampa    (mA/cm2)
           ampasyn  (1)
}


BREAKPOINT {
        ena = (1000)*R*(celsius+273.15)/F*log(nao/nai)
	ina = ratio*ampasyn*gampa*(v-ena)
	ik = ratio*ampasyn*gampak*(v-ek)
        iampa = ina + ik 
}
