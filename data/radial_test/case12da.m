function mpc = case12da
%CASE12DA

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 1;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	0	11	1	1	1;
	2	1	0.06	0.06	0	0	1	1	0	11	1	1.1	0.9;
	3	1	0.04	0.03	0	0	1	1	0	11	1	1.1	0.9;
	4	1	0.055	0.055	0	0	1	1	0	11	1	1.1	0.9;
	5	1	0.03	0.03	0	0	1	1	0	11	1	1.1	0.9;
	6	1	0.02	0.015	0	0	1	1	0	11	1	1.1	0.9;
	7	1	0.055	0.055	0	0	1	1	0	11	1	1.1	0.9;
	8	1	0.045	0.045	0	0	1	1	0	11	1	1.1	0.9;
	9	1	0.04	0.04	0	0	1	1	0	11	1	1.1	0.9;
	10	1	0.035	0.03	0	0	1	1	0	11	1	1.1	0.9;
	11	1	0.04	0.03	0	0	1	1	0	11	1	1.1	0.9;
	12	1	0.015	0.015	0	0	1	1	0	11	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0	10	-10	1	100	1	10	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.00903305785	0.00376033058	0	0	0	0	0	0	1	-360	360;
	2	3	0.00978512397	0.00408264463	0	0	0	0	0	0	1	-360	360;
	3	4	0.0173140496	0.00721487603	0	0	0	0	0	0	1	-360	360;
	4	5	0.0263471074	0.0109834711	0	0	0	0	0	0	1	-360	360;
	5	6	0.00903305785	0.00376033058	0	0	0	0	0	0	1	-360	360;
	6	7	0.00828099174	0.00344628099	0	0	0	0	0	0	1	-360	360;
	7	8	0.0363884298	0.0100413223	0	0	0	0	0	0	1	-360	360;
	8	9	0.0466280992	0.0131983471	0	0	0	0	0	0	1	-360	360;
	9	10	0.0238842975	0.00676033058	0	0	0	0	0	0	1	-360	360;
	10	11	0.0125123967	0.00353719008	0	0	0	0	0	0	1	-360	360;
	11	12	0.010231405	0.00290082645	0	0	0	0	0	0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	3	0	20	0;
];
