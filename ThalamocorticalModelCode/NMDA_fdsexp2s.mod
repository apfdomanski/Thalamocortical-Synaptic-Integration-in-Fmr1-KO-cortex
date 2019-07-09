COMMENT
NMDA channel
This is an adapted version of Exp2Syn.
Adapted by Kevin M Biddell similar to as described by wolf et al 2006
4/21/07
verified 3/29/2012 KMB 
kevin.biddell@gmail.com

Voltage dependence of Mg2+ block:
	Jahr & Stevens 1990. J Neurosci 10: 1830.
	Jahr & Stevens 1990. J Neurosci 10: 3178.

Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 -> 0 then we have a alphasynapse.
and if tau1 -> 0 then we have just single exponential decay.

The factor is evaluated in the
initial block such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT

NEURON {
	POINT_PROCESS NMDA_FDSExp2Syn
	RANGE tau1, tau2, g, e, i, mg, f, tau_F, d1, tau_D1, d2, tau_D2
	NONSPECIFIC_CURRENT i
	GLOBAL total,vmin, vmax
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}	

PARAMETER {
	e	= 0    (mV)	: reversal potential
	tau1	= 2.23  (ms)<1e-9,1e9>	: changed by KMB 5/07/07 from 2.82 to 20.77% less (wolf et 2006)
	tau2	= 150   (ms)<1e-9,1e9> 	: changed by KMB 4/23/07 from 160 to 52.7% less (wolf et al 2006)
	mg	= 1.2    (mM)	: external magnesium concentration
	vmin = -120	(mV)
	vmax = 100	(mV)
	f = 0.917 (1) < 0, 1e9 >    : facilitation
        tau_F = 94 (ms) < 1e-9, 1e9 >
        d1 = 0.416 (1) < 0, 1 >     : fast depression
        tau_D1 = 380 (ms) < 1e-9, 1e9 >
        d2 = 0.975 (1) < 0, 1 >     : slow depression
        tau_D2 = 9200 (ms) < 1e-9, 1e9 >
}

ASSIGNED {
	v (mV)
	i (nA)
	g (umho)
	factor
	total (umho)
}

STATE {
	A (umho)
	B (umho)
	Fraction_active		: fraction free of Mg2+ block
}

INITIAL {
	LOCAL tp
	rates(v)
	total = 0
	if (tau1/tau2 > 0.9999) {
		tau1 = 0.9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
}

BREAKPOINT {
	rates(v)
	SOLVE state METHOD cnexp
        g = (B-A) 	
        i = g*Fraction_active*(v - e) : N.B. C deals with fraction free of Mg++ block	
}


DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}
	
PROCEDURE rates(v(mV)) {
	TABLE Fraction_active
	DEPEND mg
	FROM vmin TO vmax WITH 200

	: from Jahr & Stevens

	Fraction_active = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}


NET_RECEIVE(weight (umho), F, D1, D2, tsyn (ms)) {
INITIAL {
: these are in NET_RECEIVE to be per-stream
        F = 1
        D1 = 1
        D2 = 1
        tsyn = t
: this header will appear once per stream
: printf("t\t t-tsyn\t F\t D1\t D2\t amp\t newF\t newD1\t newD2\n")
}

        F = 1 + (F-1)*exp(-(t - tsyn)/tau_F)
        D1 = 1 - (1-D1)*exp(-(t - tsyn)/tau_D1)
        D2 = 1 - (1-D2)*exp(-(t - tsyn)/tau_D2)
: printf("%g\t%g\t%g\t%g\t%g\t%g", t, t-tsyn, F, D1, D2, weight*F*D1*D2)
        tsyn = t

	state_discontinuity(A, A + weight*factor*F*D1*D2)
	state_discontinuity(B, B + weight*factor*F*D1*D2)
	total = total+weight*F*D1*D2

        F = F + f
        D1 = D1 * d1
        D2 = D2 * d2
: printf("\t%g\t%g\t%g\n", F, D1, D2)
}
















