function mpc = case51ga
%CASE51GA Power flow data for 51 bus distribution system from Gampa & Das
%   Please see CASEFORMAT for details on the case file format.
%
%   Data from ...
%       Gampa SR, Das D (2015) Optimum placement and sizing of DGs considering
%       average hourly variations of load. Int J Electr Power Energy Syst
%       66:25-40. doi: 10.1016/j.ijepes.2014.10.047
%       URL: https://doi.org/10.1016/j.ijepes.2014.10.047

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 1;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [ %% (Pd and Qd are specified in kW & kVAr here, converted to MW & MVAr below)
	1	3	0	0	0	0	1	1	0	11	1	1	1;
	2	1	40	30	0	0	1	1	0	11	1	1.1	0.9;
	3	1	60	40	0	0	1	1	0	11	1	1.1	0.9;
	4	1	20	10	0	0	1	1	0	11	1	1.1	0.9;
	5	1	80	60	0	0	1	1	0	11	1	1.1	0.9;
	6	1	38	18	0	0	1	1	0	11	1	1.1	0.9;
	7	1	20	15	0	0	1	1	0	11	1	1.1	0.9;
	8	1	60	40	0	0	1	1	0	11	1	1.1	0.9;
	9	1	70	45	0	0	1	1	0	11	1	1.1	0.9;
	10	1	60	35	0	0	1	1	0	11	1	1.1	0.9;
	11	1	80	50	0	0	1	1	0	11	1	1.1	0.9;
	12	1	10	5	0	0	1	1	0	11	1	1.1	0.9;
	13	1	25	15	0	0	1	1	0	11	1	1.1	0.9;
	14	1	55	45	0	0	1	1	0	11	1	1.1	0.9;
	15	1	120	80	0	0	1	1	0	11	1	1.1	0.9;
	16	1	40	25	0	0	1	1	0	11	1	1.1	0.9;
	17	1	35	25	0	0	1	1	0	11	1	1.1	0.9;
	18	1	60	30	0	0	1	1	0	11	1	1.1	0.9;
	19	1	80	50	0	0	1	1	0	11	1	1.1	0.9;
	20	1	60	35	0	0	1	1	0	11	1	1.1	0.9;
	21	1	50	30	0	0	1	1	0	11	1	1.1	0.9;
	22	1	50	30	0	0	1	1	0	11	1	1.1	0.9;
	23	1	80	60	0	0	1	1	0	11	1	1.1	0.9;
	24	1	45	25	0	0	1	1	0	11	1	1.1	0.9;
	25	1	38	18	0	0	1	1	0	11	1	1.1	0.9;
	26	1	78	48	0	0	1	1	0	11	1	1.1	0.9;
	27	1	16	8	0	0	1	1	0	11	1	1.1	0.9;
	28	1	18	10	0	0	1	1	0	11	1	1.1	0.9;
	29	1	40	30	0	0	1	1	0	11	1	1.1	0.9;
	30	1	40	30	0	0	1	1	0	11	1	1.1	0.9;
	31	1	20	15	0	0	1	1	0	11	1	1.1	0.9;
	32	1	30	20	0	0	1	1	0	11	1	1.1	0.9;
	33	1	36	26	0	0	1	1	0	11	1	1.1	0.9;
	34	1	50	40	0	0	1	1	0	11	1	1.1	0.9;
	35	1	27	18	0	0	1	1	0	11	1	1.1	0.9;
	36	1	33	16	0	0	1	1	0	11	1	1.1	0.9;
	37	1	42	22	0	0	1	1	0	11	1	1.1	0.9;
	38	1	55	30	0	0	1	1	0	11	1	1.1	0.9;
	39	1	44	26	0	0	1	1	0	11	1	1.1	0.9;
	40	1	80	70	0	0	1	1	0	11	1	1.1	0.9;
	41	1	60	30	0	0	1	1	0	11	1	1.1	0.9;
	42	1	45	30	0	0	1	1	0	11	1	1.1	0.9;
	43	1	48	28	0	0	1	1	0	11	1	1.1	0.9;
	44	1	68	38	0	0	1	1	0	11	1	1.1	0.9;
	45	1	77	23	0	0	1	1	0	11	1	1.1	0.9;
	46	1	60	30	0	0	1	1	0	11	1	1.1	0.9;
	47	1	40	20	0	0	1	1	0	11	1	1.1	0.9;
	48	1	45	45	0	0	1	1	0	11	1	1.1	0.9;
	49	1	70	50	0	0	1	1	0	11	1	1.1	0.9;
	50	1	30	20	0	0	1	1	0	11	1	1.1	0.9;
	51	1	35	30	0	0	1	1	0	11	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0	10	-10	1	100	1	10	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [  %% (r and x specified in ohms here, converted to p.u. below)
	1	2	0.2740	0.3560	0	0	0	0	0	0	1	-360	360;
	2	3	0.1370	0.1780	0	0	0	0	0	0	1	-360	360;
	3	4	0.3288	0.4272	0	0	0	0	0	0	1	-360	360;
	4	5	0.1096	0.1424	0	0	0	0	0	0	1	-360	360;
	5	6	0.5400	0.4356	0	0	0	0	0	0	1	-360	360;
	6	7	0.3600	0.2904	0	0	0	0	0	0	1	-360	360;
	7	8	0.3600	0.2904	0	0	0	0	0	0	1	-360	360;
	8	9	0.7200	0.5808	0	0	0	0	0	0	1	-360	360;
	9	10	2.7320	0.7792	0	0	0	0	0	0	1	-360	360;
	10	11	2.0490	0.5844	0	0	0	0	0	0	1	-360	360;
	11	12	2.0490	0.5844	0	0	0	0	0	0	1	-360	360;
	12	13	0.9562	0.2727	0	0	0	0	0	0	1	-360	360;
	13	14	1.0928	0.3117	0	0	0	0	0	0	1	-360	360;
	14	15	1.5026	0.4286	0	0	0	0	0	0	1	-360	360;
	15	16	3.0052	0.8571	0	0	0	0	0	0	1	-360	360;
	3	17	2.7320	0.7792	0	0	0	0	0	0	1	-360	360;
	17	18	0.8196	0.2338	0	0	0	0	0	0	1	-360	360;
	18	19	1.3660	0.3896	0	0	0	0	0	0	1	-360	360;
	19	20	1.3660	0.3896	0	0	0	0	0	0	1	-360	360;
	20	21	2.0490	0.5844	0	0	0	0	0	0	1	-360	360;
	21	22	1.5026	0.4286	0	0	0	0	0	0	1	-360	360;
	4	23	1.6392	0.4675	0	0	0	0	0	0	1	-360	360;
	23	24	1.7758	0.5065	0	0	0	0	0	0	1	-360	360;
	24	25	1.0928	0.3117	0	0	0	0	0	0	1	-360	360;
	25	26	0.8196	0.2338	0	0	0	0	0	0	1	-360	360;
	26	27	0.5464	0.1558	0	0	0	0	0	0	1	-360	360;
	27	28	1.0928	0.3117	0	0	0	0	0	0	1	-360	360;
	28	29	0.2732	0.0779	0	0	0	0	0	0	1	-360	360;
	6	30	0.7020	0.4774	0	0	0	0	0	0	1	-360	360;
	30	31	0.6480	0.4406	0	0	0	0	0	0	1	-360	360;
	31	32	0.6480	0.4406	0	0	0	0	0	0	1	-360	360;
	32	33	0.6480	0.4406	0	0	0	0	0	0	1	-360	360;
	33	34	0.5400	0.3672	0	0	0	0	0	0	1	-360	360;
	34	35	0.3240	0.2203	0	0	0	0	0	0	1	-360	360;
	35	36	0.3888	0.2644	0	0	0	0	0	0	1	-360	360;
	36	37	0.4320	0.2938	0	0	0	0	0	0	1	-360	360;
	37	38	0.5940	0.4039	0	0	0	0	0	0	1	-360	360;
	38	39	0.7020	0.4774	0	0	0	0	0	0	1	-360	360;
	7	40	1.9124	0.5454	0	0	0	0	0	0	1	-360	360;
	40	41	3.0052	0.8571	0	0	0	0	0	0	1	-360	360;
	41	42	2.4588	0.7013	0	0	0	0	0	0	1	-360	360;
	42	43	2.1856	0.6234	0	0	0	0	0	0	1	-360	360;
	43	44	2.1856	0.6234	0	0	0	0	0	0	1	-360	360;
	44	45	0.6830	0.1948	0	0	0	0	0	0	1	-360	360;
	9	46	0.9562	0.2727	0	0	0	0	0	0	1	-360	360;
	46	47	1.0245	0.2922	0	0	0	0	0	0	1	-360	360;
	47	48	1.2294	0.3506	0	0	0	0	0	0	1	-360	360;
	48	49	1.7758	0.5065	0	0	0	0	0	0	1	-360	360;
	49	50	1.6392	0.4675	0	0	0	0	0	0	1	-360	360;
	50	51	1.3660	0.3896	0	0	0	0	0	0	1	-360	360;
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
