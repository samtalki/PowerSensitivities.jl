function mpc = case10ba
%CASE10BA

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 10;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	0	23	1	1	1;
	2	1	1.84	0.46	0	0	1	1	0	23	1	1.1	0.9;
	3	1	0.98	0.34	0	0	1	1	0	23	1	1.1	0.9;
	4	1	1.79	0.446	0	0	1	1	0	23	1	1.1	0.9;
	5	1	1.598	1.84	0	0	1	1	0	23	1	1.1	0.9;
	6	1	1.61	0.6	0	0	1	1	0	23	1	1.1	0.9;
	7	1	0.78	0.11	0	0	1	1	0	23	1	1.1	0.9;
	8	1	1.15	0.06	0	0	1	1	0	23	1	1.1	0.9;
	9	1	0.98	0.13	0	0	1	1	0	23	1	1.1	0.9;
	10	1	1.64	0.2	0	0	1	1	0	23	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0	10	-10	1	100	1	10	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.00233081285	0.00780151229	0	0	0	0	0	0	1	-360	360;
	2	3	0.000264650284	0.0114385633	0	0	0	0	0	0	1	-360	360;
	3	4	0.0141077505	0.022778828	0	0	0	0	0	0	1	-360	360;
	4	5	0.0132022684	0.0115009452	0	0	0	0	0	0	1	-360	360;
	5	6	0.0374877127	0.032657845	0	0	0	0	0	0	1	-360	360;
	6	7	0.0171134216	0.0149073724	0	0	0	0	0	0	1	-360	360;
	7	8	0.0388506616	0.0220037807	0	0	0	0	0	0	1	-360	360;
	8	9	0.0906483932	0.051342155	0	0	0	0	0	0	1	-360	360;
	9	10	0.101009452	0.0572098299	0	0	0	0	0	0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	3	0	20	0;
];
