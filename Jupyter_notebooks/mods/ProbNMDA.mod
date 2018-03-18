TITLE NMDA receptor with presynaptic short-term plasticity 


COMMENT
NMDA receptor conductance using a dual-exponential profile
Presynaptic short-term plasticity based on Fuhrmann et al, 2002
Implemented by Srikanth Ramaswamy, Blue Brain Project, March 2009
ENDCOMMENT

NEURON {
    THREADSAFE
	POINT_PROCESS ProbNMDA
	RANGE tau_r, tau_d
	RANGE Use, u, Dep, Fac, u0
	RANGE i, g, e
	RANGE mg,mggate
    NONSPECIFIC_CURRENT i
}

PARAMETER {
	tau_r = 0.29	 	(ms) : dual-exponential conductance profile
	tau_d = 43		(ms) : IMPORTANT: tau_r < tau_d
	Use = 0.67 		(1)      : Utilization of synaptic efficacy (just initial values! Use,Dep and Fac are overwritten by BlueBuilder assigned values) 	
	Dep = 800 		(ms) 	 : relaxation time constant from depression
	Fac = 3 		(ms)     :  relaxation time constant from facilitation
    e  = 0     (mV)  : NMDA reversal potential	
    gmax = 0.001     (us) : weight conversion factor (from nS to uS)
	mg = 1			(mM)  : initial concentration of mg2+
	mggate       
    u0 = 0      :initial value of u, which is the running value of Use
}

COMMENT
The Verbatim block is needed to generate random nos. from a uniform distribution between 0 and 1 
for comparison with Pr to decide whether to activate the synapse or not
ENDCOMMENT
   
VERBATIM
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

float ranfNMDA(){
        //double MAX = (double)RAND_MAX;
        double r = (rand() / (double) RAND_MAX);
        return r;
}

void SetSeedNowNMDA(){
#ifdef SYN_DEBUG
srand(time(NULL));
#else
srand(888);
#endif
return;
}
ENDVERBATIM
   


ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
    factor
}

STATE {
	A  : state variable to construct the dual-exponential profile - decays with conductance tau_r
	B  : state variable to construct the dual-exponential profile - decays with conductance tau_d
}

INITIAL{
	LOCAL tp
	A = 0
	B = 0
	tp = (tau_r*tau_d)/(tau_d-tau_r)*log(tau_d/tau_r) :time to peak of the conductance
	factor = -exp(-tp/tau_r)+exp(-tp/tau_d) :Normalization factor - so that when t = tp, gsyn = gpeak
	factor = 1/factor
    SetSeedNowNMDA()
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	mggate = 1 / (1 + exp(0.062 (/mV) * -(v)) * (mg / 3.57 (mM))) :mggate kinetics - Jahr & Stevens 1990
	g = gmax*(B - A)*mggate
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau_r
	B' = -B/tau_d
}

NET_RECEIVE (weight, Pv, Pr, u, tsyn (ms)){
weight = weight*0.71 
:printf("weight NMDA = %g \n ",weight)
:the NETCON.weight = gsyn (per synaptic contact) * scaling factor * 0.71, as gNMDA/gAMPA = 0.71 from Chaelon et al. 2003, and Markram et al. 1997
    INITIAL{
		Pv=1
		u=u0
		tsyn=t
	    }
        : calc u at event-
    	if (Fac > 0) {
	      	u = u*exp(-(t - tsyn)/Fac) :update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.
	   } else {
		  u = Use  
	   } 
	   if(Fac > 0){
		  u = u + Use*(1-u) :update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.
	   }	

        
            Pv   = 1 - (1-Pv) * exp(-(t-tsyn)/Dep) :Probability Pv for a vesicle to be available for release, analogous to the pool of synaptic
                                                 :resources available for release in the deterministic model. Eq. 3 in Fuhrmann et al.
            Pr  = u * Pv                         :Pr is calculated as Pv * u (running value of Use)
            Pv   = Pv - u * Pv                   :update Pv as per Eq. 3 in Fuhrmann et al.
            :printf("Pv = %g\n", Pv)
            :printf("Pr = %g\n", Pr)
                tsyn = t
                if (ranfNMDA() < Pr){
	   	        A = A + weight*factor
	            B = B + weight*factor
                }
}    
		
