TITLE AMPA receptor as a distributed mechanism
 

UNITS {
        (pA) = (picoamp)
        (molar) = (1/liter)
	(mV) =	(millivolt)
        (S) = (siemens)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
        F = (faraday) (coulomb)
        R = (k-mole) (joule/degC)
}


INDEPENDENT {v FROM -100 TO 50 WITH 50 (mV)}

NEURON {
	SUFFIX ampa
	USEION na READ nai,ena  WRITE ina
	USEION k  WRITE ik
	RANGE  iampa,ina,ik,gampa,gampak,nai,ratio
 
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
        gampa =  0.0e-6 (S/cm2) :2.68e-6*3.626/2.2237
        gampak = 0.0e-6 (S/cm2) :3.37e-6*3.626/2.2237
        ek  = -100 (mV)
        nao =  145 (mM) 
        

}

ASSIGNED { 
           ina	  (mA/cm2)
           ik	  (mA/cm2)
           iampa  (mA/cm2)
}


BREAKPOINT {
        ena = (1000)*R*(celsius+273.15)/F*log(nao/nai)
	ina = ratio*gampa*(v-ena)
	ik = ratio*gampak*(v-ek)
        iampa= ina + ik 
}
