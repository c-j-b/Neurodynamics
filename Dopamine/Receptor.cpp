
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "time.h"


typedef struct {
	float X1;
	float X2;
	float X3;
}*UserData;

main()
{
	FILE *rp, *ap, *gaba_eventp, *glu_eventp;
	srand(time(NULL));
	int counter = 0, inc = 0;
	double alpha1, beta1, tau1, alpha2, beta2, tau2, rinf1, rinf2, meanIEI; 
	double t = 0.0, endtime = 90000.0, scale = 0.1;
	double cut = 1000.0;
	double duration;
double a,b,c;
	int check2;
	float glu_stimtime;
	
	meanIEI = 2.2237;
	alpha1 = 1.1;
	beta1 = 0.19;
	alpha2 = 0.072;
	beta2 = 0.0066;
	tau1 = alpha1 + beta1;
	tau2 = alpha2 + beta2;
	rinf1 = alpha1/tau1;
	rinf2 = alpha2/tau2;
	duration = 1.0;
	
	UserData data;
	data = (UserData)malloc(sizeof *data);
	
	
	
	glu_eventp = fopen("glu_event.data","w");
	glu_stimtime = 0;
	do
	{ 
               a = rand();
               a=1.0*a;
		b=1.0*(1.0*RAND_MAX+1.0);
		c = -meanIEI*log(a/b);
		glu_stimtime = glu_stimtime + c;
		fprintf(glu_eventp, "%f\n", glu_stimtime);
		
	}
	while(glu_stimtime < endtime);
	fclose(glu_eventp);
	
	
	rp = fopen("AMPA.data", "w");
	ap = fopen("NMDA.data", "w");
	
	while(t < endtime)
	{
		t += scale;
		data->X1 = 0.0;
		data->X2 = 0.0;
		data->X3 = 0.0;
		glu_eventp = fopen("glu_event.data","r");
		rewind(glu_eventp);
		check2 = fscanf(glu_eventp, "%f\n", &glu_stimtime);
		
		while(check2 != EOF)
		{
			if(t - glu_stimtime > 0 && t - glu_stimtime <= duration)
			{
				data->X1 = data->X1 + rinf1*(1.0 - exp(-tau1*(t-glu_stimtime)));
				data->X2 = data->X2 + rinf2*(1.0 - exp(-tau2*(t-glu_stimtime)));
			}
			
			if (t - glu_stimtime > duration && t - glu_stimtime < cut)
			{
				data->X1 = data->X1 + rinf1*(1.0 - exp(-tau1*duration))*exp(-beta1*(t-(glu_stimtime+duration)));
				data->X2 = data->X2 + rinf2*(1.0 - exp(-tau2*duration))*exp(-beta2*(t-(glu_stimtime+duration)));
			}
			
			if(glu_stimtime > t)
			{
				fclose(glu_eventp);
				break;
			}
			
			check2 = fscanf(glu_eventp, "%f\n", &glu_stimtime);
			if(check2==EOF) 
			{	
				glu_stimtime = endtime + 1.0;
			}		
		}
		fprintf(rp, "%f\n", data->X1);
		fprintf(ap, "%f\n", data->X2);
	}
	fclose(rp);
	fclose(ap);

	return(0);
}

