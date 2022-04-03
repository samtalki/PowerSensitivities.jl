function mpc = case74ds
%CASE74DS  Power flow data for 74 bus distribution system from Myint & Naing
%   Please see CASEFORMAT for details on the case file format.
%
%   Data from ...
%       Myint SM, Naing SW (2015) "Network Reconfiguration For Loss Reduction
%       And Voltage Stability Improvement Of 74-Bus Radial Distribution System
%       Using Particle Swarm Optimization Algorithm", International Journal of
%       Electrical, Electronics and Data Communication, ISSN: 2320-2084
%       Volume 3, Issue 6, June 2015.

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 1;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [ %% (Pd and Qd are specified in kW & kVAr here, converted to MW & MVAr below)
	1	3	0	0	0	0	1	1	0	11	1	1	1;
	2	1	193	130	0	0	1	1	0	11	1	1.1	0.9;
	3	1	1	1	0	0	1	1	0	11	1	1.1	0.9;
	4	1	3	2	0	0	1	1	0	11	1	1.1	0.9;
	5	1	108	73	0	0	1	1	0	11	1	1.1	0.9;
	6	1	5	3	0	0	1	1	0	11	1	1.1	0.9;
	7	1	18	12	0	0	1	1	0	11	1	1.1	0.9;
	8	1	16	11	0	0	1	1	0	11	1	1.1	0.9;
	9	1	76	51	0	0	1	1	0	11	1	1.1	0.9;
	10	1	52	35	0	0	1	1	0	11	1	1.1	0.9;
	11	1	46	31	0	0	1	1	0	11	1	1.1	0.9;
	12	1	123	83	0	0	1	1	0	11	1	1.1	0.9;
	13	1	46	31	0	0	1	1	0	11	1	1.1	0.9;
	14	1	7	4	0	0	1	1	0	11	1	1.1	0.9;
	15	1	70	47	0	0	1	1	0	11	1	1.1	0.9;
	16	1	70	47	0	0	1	1	0	11	1	1.1	0.9;
	17	1	4	3	0	0	1	1	0	11	1	1.1	0.9;
	18	1	49	33	0	0	1	1	0	11	1	1.1	0.9;
	19	1	57	38	0	0	1	1	0	11	1	1.1	0.9;
	20	1	118	79	0	0	1	1	0	11	1	1.1	0.9;
	21	1	19	13	0	0	1	1	0	11	1	1.1	0.9;
	22	1	14	10	0	0	1	1	0	11	1	1.1	0.9;
	23	1	8	6	0	0	1	1	0	11	1	1.1	0.9;
	24	1	210	141	0	0	1	1	0	11	1	1.1	0.9;
	25	1	136	92	0	0	1	1	0	11	1	1.1	0.9;
	26	1	189	127	0	0	1	1	0	11	1	1.1	0.9;
	27	1	16	10	0	0	1	1	0	11	1	1.1	0.9;
	28	1	74	49	0	0	1	1	0	11	1	1.1	0.9;
	29	1	117	79	0	0	1	1	0	11	1	1.1	0.9;
	30	1	6	4	0	0	1	1	0	11	1	1.1	0.9;
	31	1	125	84	0	0	1	1	0	11	1	1.1	0.9;
	32	1	17	11	0	0	1	1	0	11	1	1.1	0.9;
	33	1	4	3	0	0	1	1	0	11	1	1.1	0.9;
	34	1	19	13	0	0	1	1	0	11	1	1.1	0.9;
	35	1	31	21	0	0	1	1	0	11	1	1.1	0.9;
	36	1	40	27	0	0	1	1	0	11	1	1.1	0.9;
	37	1	11	7	0	0	1	1	0	11	1	1.1	0.9;
	38	1	217	146	0	0	1	1	0	11	1	1.1	0.9;
	39	1	37	25	0	0	1	1	0	11	1	1.1	0.9;
	40	1	114	77	0	0	1	1	0	11	1	1.1	0.9;
	41	1	12	8	0	0	1	1	0	11	1	1.1	0.9;
	42	1	55	37	0	0	1	1	0	11	1	1.1	0.9;
	43	1	72	48	0	0	1	1	0	11	1	1.1	0.9;
	44	1	40	27	0	0	1	1	0	11	1	1.1	0.9;
	45	1	4	3	0	0	1	1	0	11	1	1.1	0.9;
	46	1	41	27	0	0	1	1	0	11	1	1.1	0.9;
	47	1	42	28	0	0	1	1	0	11	1	1.1	0.9;
	48	1	49	33	0	0	1	1	0	11	1	1.1	0.9;
	49	1	43	29	0	0	1	1	0	11	1	1.1	0.9;
	50	1	61	41	0	0	1	1	0	11	1	1.1	0.9;
	51	1	46	31	0	0	1	1	0	11	1	1.1	0.9;
	52	1	10	7	0	0	1	1	0	11	1	1.1	0.9;
	53	1	9	6	0	0	1	1	0	11	1	1.1	0.9;
	54	1	65	44	0	0	1	1	0	11	1	1.1	0.9;
	55	1	45	30	0	0	1	1	0	11	1	1.1	0.9;
	56	1	76	51	0	0	1	1	0	11	1	1.1	0.9;
	57	1	71	48	0	0	1	1	0	11	1	1.1	0.9;
	58	1	394	265	0	0	1	1	0	11	1	1.1	0.9;
	59	1	234	157	0	0	1	1	0	11	1	1.1	0.9;
	60	1	227	153	0	0	1	1	0	11	1	1.1	0.9;
	61	1	225	151	0	0	1	1	0	11	1	1.1	0.9;
	62	1	192	129	0	0	1	1	0	11	1	1.1	0.9;
	63	1	118	79	0	0	1	1	0	11	1	1.1	0.9;
	64	1	418	281	0	0	1	1	0	11	1	1.1	0.9;
	65	1	234	157	0	0	1	1	0	11	1	1.1	0.9;
	66	1	221	148	0	0	1	1	0	11	1	1.1	0.9;
	67	1	198	133	0	0	1	1	0	11	1	1.1	0.9;
	68	1	249	167	0	0	1	1	0	11	1	1.1	0.9;
	69	1	113	76	0	0	1	1	0	11	1	1.1	0.9;
	70	1	96	64	0	0	1	1	0	11	1	1.1	0.9;
	71	1	91	61	0	0	1	1	0	11	1	1.1	0.9;
	72	1	116	78	0	0	1	1	0	11	1	1.1	0.9;
	73	1	168	113	0	0	1	1	0	11	1	1.1	0.9;
	74	1	116	78	0	0	1	1	0	11	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0	10	-10	1	100	1	10	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [  %% (r and x specified in ohms here, converted to p.u. below)
	1	2	0.0578	0.0438	0	0	0	0	0	0	1	-360	360;
	2	3	0.0564	0.0428	0	0	0	0	0	0	1	-360	360;
	3	4	0.0431	0.0327	0	0	0	0	0	0	1	-360	360;
	4	5	0.0431	0.0327	0	0	0	0	0	0	1	-360	360;
	5	6	0.036	0.0273	0	0	0	0	0	0	1	-360	360;
	6	7	0.036	0.0273	0	0	0	0	0	0	1	-360	360;
	7	8	0.0289	0.0219	0	0	0	0	0	0	1	-360	360;
	8	9	0.0289	0.0219	0	0	0	0	0	0	1	-360	360;
	9	10	0.0289	0.0219	0	0	0	0	0	0	1	-360	360;
	10	11	0.0289	0.0219	0	0	0	0	0	0	1	-360	360;
	11	12	0.0289	0.0219	0	0	0	0	0	0	1	-360	360;
	12	13	0.0289	0.0219	0	0	0	0	0	0	1	-360	360;
	13	14	0.0289	0.0219	0	0	0	0	0	0	1	-360	360;
	14	15	0.0289	0.0219	0	0	0	0	0	0	1	-360	360;
	15	16	0.0465	0.0352	0	0	0	0	0	0	1	-360	360;
	16	17	0.0465	0.0352	0	0	0	0	0	0	1	-360	360;
	17	18	0.0465	0.0352	0	0	0	0	0	0	1	-360	360;
	18	19	0.0465	0.0352	0	0	0	0	0	0	1	-360	360;
	19	20	0.0403	0.0305	0	0	0	0	0	0	1	-360	360;
	20	21	0.0403	0.0305	0	0	0	0	0	0	1	-360	360;
	21	22	0.0403	0.0305	0	0	0	0	0	0	1	-360	360;
	22	23	0.0403	0.0305	0	0	0	0	0	0	1	-360	360;
	23	24	0.0507	0.0384	0	0	0	0	0	0	1	-360	360;
	24	25	0.036	0.0273	0	0	0	0	0	0	1	-360	360;
	25	26	0.045	0.0341	0	0	0	0	0	0	1	-360	360;
	26	27	0.045	0.0341	0	0	0	0	0	0	1	-360	360;
	27	28	0.045	0.0341	0	0	0	0	0	0	1	-360	360;
	28	29	0.0289	0.0219	0	0	0	0	0	0	1	-360	360;
	29	30	0.0289	0.0219	0	0	0	0	0	0	1	-360	360;
	30	31	0.0431	0.0327	0	0	0	0	0	0	1	-360	360;
	31	32	0.0403	0.0305	0	0	0	0	0	0	1	-360	360;
	32	33	0.036	0.0273	0	0	0	0	0	0	1	-360	360;
	33	34	0.0431	0.0327	0	0	0	0	0	0	1	-360	360;
	34	35	0.0465	0.0352	0	0	0	0	0	0	1	-360	360;
	35	36	0.0465	0.0352	0	0	0	0	0	0	1	-360	360;
	36	37	0.0484	0.0366	0	0	0	0	0	0	1	-360	360;
	37	38	0.0484	0.0366	0	0	0	0	0	0	1	-360	360;
	38	39	0.0431	0.0327	0	0	0	0	0	0	1	-360	360;
	39	40	0.0431	0.0327	0	0	0	0	0	0	1	-360	360;
	40	41	0.036	0.0273	0	0	0	0	0	0	1	-360	360;
	41	42	0.036	0.0273	0	0	0	0	0	0	1	-360	360;
	42	43	0.0578	0.0438	0	0	0	0	0	0	1	-360	360;
	43	44	0.0465	0.0352	0	0	0	0	0	0	1	-360	360;
	44	45	0.0507	0.0384	0	0	0	0	0	0	1	-360	360;
	45	46	0.0536	0.0406	0	0	0	0	0	0	1	-360	360;
	46	47	0.0549	0.0417	0	0	0	0	0	0	1	-360	360;
	47	48	0.0564	0.0428	0	0	0	0	0	0	1	-360	360;
	48	49	0.0422	0.0319	0	0	0	0	0	0	1	-360	360;
	49	50	0.0431	0.0327	0	0	0	0	0	0	1	-360	360;
	50	51	0.0649	0.0492	0	0	0	0	0	0	1	-360	360;
	51	52	0.0621	0.0471	0	0	0	0	0	0	1	-360	360;
	52	53	0.0635	0.0481	0	0	0	0	0	0	1	-360	360;
	53	54	0.0635	0.0481	0	0	0	0	0	0	1	-360	360;
	54	55	0.0635	0.0481	0	0	0	0	0	0	1	-360	360;
	55	56	0.0318	0.0241	0	0	0	0	0	0	1	-360	360;
	56	57	0.0318	0.0241	0	0	0	0	0	0	1	-360	360;
	1	58	0.201	0.1523	0	0	0	0	0	0	1	-360	360;
	58	59	0.0721	0.0546	0	0	0	0	0	0	1	-360	360;
	59	60	0.0721	0.0546	0	0	0	0	0	0	1	-360	360;
	60	61	0.0721	0.0546	0	0	0	0	0	0	1	-360	360;
	61	62	0.0721	0.0546	0	0	0	0	0	0	1	-360	360;
	62	63	0.0721	0.0546	0	0	0	0	0	0	1	-360	360;
	63	64	0.0721	0.0546	0	0	0	0	0	0	1	-360	360;
	64	65	0.0721	0.0546	0	0	0	0	0	0	1	-360	360;
	65	66	0.0721	0.0546	0	0	0	0	0	0	1	-360	360;
	66	67	0.0721	0.0546	0	0	0	0	0	0	1	-360	360;
	67	68	0.0721	0.0546	0	0	0	0	0	0	1	-360	360;
	68	69	0.0721	0.0546	0	0	0	0	0	0	1	-360	360;
	69	70	0.0721	0.0546	0	0	0	0	0	0	1	-360	360;
	1	71	0.4741	0.3593	0	0	0	0	0	0	1	-360	360;
	71	72	0.2371	0.1797	0	0	0	0	0	0	1	-360	360;
	72	73	0.1659	0.1258	0	0	0	0	0	0	1	-360	360;
	73	74	0.1071	0.0812	0	0	0	0	0	0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	3	0	20	0;
];


% %% convert branch impedances from Ohms to p.u.
% [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
%     VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
% [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
%     TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
%     ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
% Vbase = mpc.bus(1, BASE_KV) * 1e3;      %% in Volts
% Sbase = mpc.baseMVA * 1e6;              %% in VA
% mpc.branch(:, [BR_R BR_X]) = mpc.branch(:, [BR_R BR_X]) / (Vbase^2 / Sbase);

% %% convert loads from kW to MW
% mpc.bus(:, [PD, QD]) = mpc.bus(:, [PD, QD]) / 1e3;
