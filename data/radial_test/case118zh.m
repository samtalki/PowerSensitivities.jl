function mpc = case118zh
%CASE118ZH

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 10;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	0	11	1	1	1;
	2	1	0.13384	0.10114	0	0	1	1	0	11	1	1.1	0.9;
	3	1	0.016214	0.011292	0	0	1	1	0	11	1	1.1	0.9;
	4	1	0.034315	0.021845	0	0	1	1	0	11	1	1.1	0.9;
	5	1	0.073016	0.063602	0	0	1	1	0	11	1	1.1	0.9;
	6	1	0.1442	0.068604	0	0	1	1	0	11	1	1.1	0.9;
	7	1	0.10447	0.061725	0	0	1	1	0	11	1	1.1	0.9;
	8	1	0.028547	0.011503	0	0	1	1	0	11	1	1.1	0.9;
	9	1	0.08756	0.051073	0	0	1	1	0	11	1	1.1	0.9;
	10	1	0.1982	0.10677	0	0	1	1	0	11	1	1.1	0.9;
	11	1	0.1468	0.075995	0	0	1	1	0	11	1	1.1	0.9;
	12	1	0.02604	0.018687	0	0	1	1	0	11	1	1.1	0.9;
	13	1	0.0521	0.02322	0	0	1	1	0	11	1	1.1	0.9;
	14	1	0.1419	0.1175	0	0	1	1	0	11	1	1.1	0.9;
	15	1	0.02187	0.02879	0	0	1	1	0	11	1	1.1	0.9;
	16	1	0.03337	0.02645	0	0	1	1	0	11	1	1.1	0.9;
	17	1	0.03243	0.02523	0	0	1	1	0	11	1	1.1	0.9;
	18	1	0.020234	0.011906	0	0	1	1	0	11	1	1.1	0.9;
	19	1	0.15694	0.078523	0	0	1	1	0	11	1	1.1	0.9;
	20	1	0.54629	0.3514	0	0	1	1	0	11	1	1.1	0.9;
	21	1	0.18031	0.1642	0	0	1	1	0	11	1	1.1	0.9;
	22	1	0.093167	0.054594	0	0	1	1	0	11	1	1.1	0.9;
	23	1	0.08518	0.03965	0	0	1	1	0	11	1	1.1	0.9;
	24	1	0.1681	0.095178	0	0	1	1	0	11	1	1.1	0.9;
	25	1	0.12511	0.15022	0	0	1	1	0	11	1	1.1	0.9;
	26	1	0.01603	0.02462	0	0	1	1	0	11	1	1.1	0.9;
	27	1	0.02603	0.02462	0	0	1	1	0	11	1	1.1	0.9;
	28	1	0.59456	0.52262	0	0	1	1	0	11	1	1.1	0.9;
	29	1	0.12062	0.059117	0	0	1	1	0	11	1	1.1	0.9;
	30	1	0.10238	0.099554	0	0	1	1	0	11	1	1.1	0.9;
	31	1	0.5134	0.3185	0	0	1	1	0	11	1	1.1	0.9;
	32	1	0.47525	0.45614	0	0	1	1	0	11	1	1.1	0.9;
	33	1	0.15143	0.13679	0	0	1	1	0	11	1	1.1	0.9;
	34	1	0.20538	0.083302	0	0	1	1	0	11	1	1.1	0.9;
	35	1	0.1316	0.093082	0	0	1	1	0	11	1	1.1	0.9;
	36	1	0.4484	0.36979	0	0	1	1	0	11	1	1.1	0.9;
	37	1	0.44052	0.32164	0	0	1	1	0	11	1	1.1	0.9;
	38	1	0.11254	0.055134	0	0	1	1	0	11	1	1.1	0.9;
	39	1	0.053963	0.038998	0	0	1	1	0	11	1	1.1	0.9;
	40	1	0.39305	0.3426	0	0	1	1	0	11	1	1.1	0.9;
	41	1	0.32674	0.27856	0	0	1	1	0	11	1	1.1	0.9;
	42	1	0.53626	0.24024	0	0	1	1	0	11	1	1.1	0.9;
	43	1	0.076247	0.066562	0	0	1	1	0	11	1	1.1	0.9;
	44	1	0.05352	0.03976	0	0	1	1	0	11	1	1.1	0.9;
	45	1	0.040328	0.031964	0	0	1	1	0	11	1	1.1	0.9;
	46	1	0.039653	0.020758	0	0	1	1	0	11	1	1.1	0.9;
	47	1	0.066195	0.042361	0	0	1	1	0	11	1	1.1	0.9;
	48	1	0.073904	0.051653	0	0	1	1	0	11	1	1.1	0.9;
	49	1	0.11477	0.057965	0	0	1	1	0	11	1	1.1	0.9;
	50	1	0.91837	1.2051	0	0	1	1	0	11	1	1.1	0.9;
	51	1	0.2103	0.14666	0	0	1	1	0	11	1	1.1	0.9;
	52	1	0.06668	0.056608	0	0	1	1	0	11	1	1.1	0.9;
	53	1	0.042207	0.040184	0	0	1	1	0	11	1	1.1	0.9;
	54	1	0.43374	0.28341	0	0	1	1	0	11	1	1.1	0.9;
	55	1	0.0621	0.02686	0	0	1	1	0	11	1	1.1	0.9;
	56	1	0.09246	0.08838	0	0	1	1	0	11	1	1.1	0.9;
	57	1	0.085188	0.055436	0	0	1	1	0	11	1	1.1	0.9;
	58	1	0.3453	0.3324	0	0	1	1	0	11	1	1.1	0.9;
	59	1	0.0225	0.01683	0	0	1	1	0	11	1	1.1	0.9;
	60	1	0.080551	0.049156	0	0	1	1	0	11	1	1.1	0.9;
	61	1	0.09586	0.090758	0	0	1	1	0	11	1	1.1	0.9;
	62	1	0.06292	0.0477	0	0	1	1	0	11	1	1.1	0.9;
	63	1	0.4788	0.46374	0	0	1	1	0	11	1	1.1	0.9;
	64	1	0.12094	0.052006	0	0	1	1	0	11	1	1.1	0.9;
	65	1	0.13911	0.10034	0	0	1	1	0	11	1	1.1	0.9;
	66	1	0.39178	0.1935	0	0	1	1	0	11	1	1.1	0.9;
	67	1	0.027741	0.026713	0	0	1	1	0	11	1	1.1	0.9;
	68	1	0.052814	0.025257	0	0	1	1	0	11	1	1.1	0.9;
	69	1	0.06689	0.038713	0	0	1	1	0	11	1	1.1	0.9;
	70	1	0.4675	0.39514	0	0	1	1	0	11	1	1.1	0.9;
	71	1	0.59485	0.23974	0	0	1	1	0	11	1	1.1	0.9;
	72	1	0.1325	0.084363	0	0	1	1	0	11	1	1.1	0.9;
	73	1	0.052699	0.022482	0	0	1	1	0	11	1	1.1	0.9;
	74	1	0.86979	0.614775	0	0	1	1	0	11	1	1.1	0.9;
	75	1	0.031349	0.029817	0	0	1	1	0	11	1	1.1	0.9;
	76	1	0.19239	0.12243	0	0	1	1	0	11	1	1.1	0.9;
	77	1	0.06575	0.04537	0	0	1	1	0	11	1	1.1	0.9;
	78	1	0.23815	0.22322	0	0	1	1	0	11	1	1.1	0.9;
	79	1	0.29455	0.16247	0	0	1	1	0	11	1	1.1	0.9;
	80	1	0.48557	0.43792	0	0	1	1	0	11	1	1.1	0.9;
	81	1	0.24353	0.18303	0	0	1	1	0	11	1	1.1	0.9;
	82	1	0.24353	0.18303	0	0	1	1	0	11	1	1.1	0.9;
	83	1	0.13425	0.11929	0	0	1	1	0	11	1	1.1	0.9;
	84	1	0.02271	0.02796	0	0	1	1	0	11	1	1.1	0.9;
	85	1	0.049513	0.026515	0	0	1	1	0	11	1	1.1	0.9;
	86	1	0.38378	0.25716	0	0	1	1	0	11	1	1.1	0.9;
	87	1	0.04964	0.0206	0	0	1	1	0	11	1	1.1	0.9;
	88	1	0.022473	0.011806	0	0	1	1	0	11	1	1.1	0.9;
	89	1	0.06293	0.04296	0	0	1	1	0	11	1	1.1	0.9;
	90	1	0.03067	0.03493	0	0	1	1	0	11	1	1.1	0.9;
	91	1	0.06253	0.06679	0	0	1	1	0	11	1	1.1	0.9;
	92	1	0.11457	0.081748	0	0	1	1	0	11	1	1.1	0.9;
	93	1	0.081292	0.066526	0	0	1	1	0	11	1	1.1	0.9;
	94	1	0.031733	0.01596	0	0	1	1	0	11	1	1.1	0.9;
	95	1	0.03332	0.06048	0	0	1	1	0	11	1	1.1	0.9;
	96	1	0.53128	0.22485	0	0	1	1	0	11	1	1.1	0.9;
	97	1	0.50703	0.36742	0	0	1	1	0	11	1	1.1	0.9;
	98	1	0.02639	0.0117	0	0	1	1	0	11	1	1.1	0.9;
	99	1	0.04599	0.030392	0	0	1	1	0	11	1	1.1	0.9;
	100	1	0.10066	0.047572	0	0	1	1	0	11	1	1.1	0.9;
	101	1	0.45648	0.3503	0	0	1	1	0	11	1	1.1	0.9;
	102	1	0.52256	0.44929	0	0	1	1	0	11	1	1.1	0.9;
	103	1	0.40843	0.16846	0	0	1	1	0	11	1	1.1	0.9;
	104	1	0.14148	0.13425	0	0	1	1	0	11	1	1.1	0.9;
	105	1	0.10443	0.066024	0	0	1	1	0	11	1	1.1	0.9;
	106	1	0.096793	0.083647	0	0	1	1	0	11	1	1.1	0.9;
	107	1	0.49392	0.41934	0	0	1	1	0	11	1	1.1	0.9;
	108	1	0.22538	0.13588	0	0	1	1	0	11	1	1.1	0.9;
	109	1	0.50921	0.38721	0	0	1	1	0	11	1	1.1	0.9;
	110	1	0.1885	0.17346	0	0	1	1	0	11	1	1.1	0.9;
	111	1	0.91803	0.89855	0	0	1	1	0	11	1	1.1	0.9;
	112	1	0.30508	0.21537	0	0	1	1	0	11	1	1.1	0.9;
	113	1	0.05438	0.04097	0	0	1	1	0	11	1	1.1	0.9;
	114	1	0.21114	0.1929	0	0	1	1	0	11	1	1.1	0.9;
	115	1	0.067009	0.053336	0	0	1	1	0	11	1	1.1	0.9;
	116	1	0.16207	0.090321	0	0	1	1	0	11	1	1.1	0.9;
	117	1	0.048785	0.029156	0	0	1	1	0	11	1	1.1	0.9;
	118	1	0.0339	0.01898	0	0	1	1	0	11	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0	10	-10	1	100	1	10	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.00297520661	0.00107107438	0	0	0	0	0	0	1	-360	360;
	2	3	0.00272727273	0.000981818182	0	0	0	0	0	0	1	-360	360;
	2	4	0.00371900826	0.00133884298	0	0	0	0	0	0	1	-360	360;
	4	5	0.00123966942	0.00446280992	0	0	0	0	0	0	1	-360	360;
	5	6	0.00123966942	0.00446280992	0	0	0	0	0	0	1	-360	360;
	6	7	0.00123966942	0.00103305785	0	0	0	0	0	0	1	-360	360;
	7	8	0.00148760331	0.00115702479	0	0	0	0	0	0	1	-360	360;
	8	9	0.00173553719	0.00520661157	0	0	0	0	0	0	1	-360	360;
	2	10	0.0137190083	0.011107438	0	0	0	0	0	0	1	-360	360;
	10	11	0.00925619835	0.00652066116	0	0	0	0	0	0	1	-360	360;
	11	12	0.0154545455	0.0258677686	0	0	0	0	0	0	1	-360	360;
	12	13	0.0117355372	0.0124958678	0	0	0	0	0	0	1	-360	360;
	13	14	0.0148760331	0.00975206612	0	0	0	0	0	0	1	-360	360;
	14	15	0.0123966942	0.00371900826	0	0	0	0	0	0	1	-360	360;
	15	16	0.0132231405	0.0148760331	0	0	0	0	0	0	1	-360	360;
	16	17	0.0129752066	0.0141322314	0	0	0	0	0	0	1	-360	360;
	11	18	0.0180165289	0.023553719	0	0	0	0	0	0	1	-360	360;
	18	19	0.00975206612	0.0152892562	0	0	0	0	0	0	1	-360	360;
	19	20	0.0132231405	0.0161983471	0	0	0	0	0	0	1	-360	360;
	20	21	0.00991735537	0.0156198347	0	0	0	0	0	0	1	-360	360;
	21	22	0.00991735537	0.00652066116	0	0	0	0	0	0	1	-360	360;
	22	23	0.116528926	0.0597520661	0	0	0	0	0	0	1	-360	360;
	23	24	0.024214876	0.0111404959	0	0	0	0	0	0	1	-360	360;
	24	25	0.0109917355	0.00859504132	0	0	0	0	0	0	1	-360	360;
	25	26	0.0147107438	0.0110743802	0	0	0	0	0	0	1	-360	360;
	26	27	0.0147107438	0.0110743802	0	0	0	0	0	0	1	-360	360;
	4	28	0.00123966942	0.00244628099	0	0	0	0	0	0	1	-360	360;
	28	29	0.000991735537	0.00228099174	0	0	0	0	0	0	1	-360	360;
	29	30	0.00991735537	0.0228595041	0	0	0	0	0	0	1	-360	360;
	30	31	0.0173553719	0.0200826446	0	0	0	0	0	0	1	-360	360;
	31	32	0.00991735537	0.00446280992	0	0	0	0	0	0	1	-360	360;
	32	33	0.0147107438	0.019338843	0	0	0	0	0	0	1	-360	360;
	33	34	0.0147107438	0.019338843	0	0	0	0	0	0	1	-360	360;
	34	35	0.0127272727	0.0133884298	0	0	0	0	0	0	1	-360	360;
	30	36	0.0154545455	0.0215702479	0	0	0	0	0	0	1	-360	360;
	36	37	0.0109917355	0.00818181818	0	0	0	0	0	0	1	-360	360;
	29	38	0.0272727273	0.0160330579	0	0	0	0	0	0	1	-360	360;
	38	39	0.0256198347	0.0160330579	0	0	0	0	0	0	1	-360	360;
	39	40	0.0107438017	0.0160330579	0	0	0	0	0	0	1	-360	360;
	40	41	0.0231404959	0.0123966942	0	0	0	0	0	0	1	-360	360;
	41	42	0.0975206612	0.0702479339	0	0	0	0	0	0	1	-360	360;
	42	43	0.0347107438	0.0201322314	0	0	0	0	0	0	1	-360	360;
	43	44	0.0223140496	0.00803305785	0	0	0	0	0	0	1	-360	360;
	44	45	0.0280165289	0.0100909091	0	0	0	0	0	0	1	-360	360;
	45	46	0.0223140496	0.0147024793	0	0	0	0	0	0	1	-360	360;
	35	47	0.0173553719	0.0114297521	0	0	0	0	0	0	1	-360	360;
	47	48	0.00991735537	0.00652066116	0	0	0	0	0	0	1	-360	360;
	48	49	0.0123966942	0.00815702479	0	0	0	0	0	0	1	-360	360;
	49	50	0.0123966942	0.00815702479	0	0	0	0	0	0	1	-360	360;
	50	51	0.0198347107	0.0130661157	0	0	0	0	0	0	1	-360	360;
	51	52	0.00991735537	0.00652066116	0	0	0	0	0	0	1	-360	360;
	52	53	0.0334710744	0.0120495868	0	0	0	0	0	0	1	-360	360;
	53	54	0.0334710744	0.0120495868	0	0	0	0	0	0	1	-360	360;
	29	55	0.0323140496	0.0116528926	0	0	0	0	0	0	1	-360	360;
	55	56	0.033553719	0.0120743802	0	0	0	0	0	0	1	-360	360;
	56	57	0.033553719	0.0120743802	0	0	0	0	0	0	1	-360	360;
	57	58	0.0583471074	0.0451322314	0	0	0	0	0	0	1	-360	360;
	58	59	0.0279338843	0.0100661157	0	0	0	0	0	0	1	-360	360;
	59	60	0.0279338843	0.0100661157	0	0	0	0	0	0	1	-360	360;
	60	61	0.017107438	0.00617355372	0	0	0	0	0	0	1	-360	360;
	61	62	0.0204132231	0.0737355372	0	0	0	0	0	0	1	-360	360;
	1	63	0.00231404959	0.00345454545	0	0	0	0	0	0	1	-360	360;
	63	64	0.00966942149	0.016661157	0	0	0	0	0	0	1	-360	360;
	64	65	0.0210743802	0.00758677686	0	0	0	0	0	0	1	-360	360;
	65	66	0.0173553719	0.00627272727	0	0	0	0	0	0	1	-360	360;
	66	67	0.0316528926	0.0114049587	0	0	0	0	0	0	1	-360	360;
	67	68	0.0416528926	0.0272975207	0	0	0	0	0	0	1	-360	360;
	68	69	0.033553719	0.0120743802	0	0	0	0	0	0	1	-360	360;
	69	70	0.0795041322	0.062892562	0	0	0	0	0	0	1	-360	360;
	70	71	0.0136363636	0.00495867769	0	0	0	0	0	0	1	-360	360;
	71	72	0.0250413223	0.00902479339	0	0	0	0	0	0	1	-360	360;
	72	73	0.0250413223	0.00902479339	0	0	0	0	0	0	1	-360	360;
	73	74	0.0170247934	0.0119008264	0	0	0	0	0	0	1	-360	360;
	74	75	0.0192561983	0.00694214876	0	0	0	0	0	0	1	-360	360;
	75	76	0.0488429752	0.0146528926	0	0	0	0	0	0	1	-360	360;
	76	77	0.0104132231	0.00374380165	0	0	0	0	0	0	1	-360	360;
	64	78	0.0461983471	0.0304710744	0	0	0	0	0	0	1	-360	360;
	78	79	0.0153719008	0.0101404959	0	0	0	0	0	0	1	-360	360;
	79	80	0.0153719008	0.0101404959	0	0	0	0	0	0	1	-360	360;
	80	81	0.0214876033	0.0114876033	0	0	0	0	0	0	1	-360	360;
	81	82	0.0127272727	0.012231405	0	0	0	0	0	0	1	-360	360;
	82	83	0.0190082645	0.0105785124	0	0	0	0	0	0	1	-360	360;
	83	84	0.0208264463	0.00876033058	0	0	0	0	0	0	1	-360	360;
	84	85	0.0148760331	0.012231405	0	0	0	0	0	0	1	-360	360;
	79	86	0.0132231405	0.0150413223	0	0	0	0	0	0	1	-360	360;
	86	87	0.0165289256	0.0190082645	0	0	0	0	0	0	1	-360	360;
	87	88	0.0132231405	0.0324793388	0	0	0	0	0	0	1	-360	360;
	65	89	0.0552892562	0.0199338843	0	0	0	0	0	0	1	-360	360;
	89	90	0.0219834711	0.0101404959	0	0	0	0	0	0	1	-360	360;
	90	91	0.0219834711	0.0101404959	0	0	0	0	0	0	1	-360	360;
	91	92	0.0219834711	0.0101404959	0	0	0	0	0	0	1	-360	360;
	92	93	0.0219834711	0.0101404959	0	0	0	0	0	0	1	-360	360;
	93	94	0.0192561983	0.00950413223	0	0	0	0	0	0	1	-360	360;
	94	95	0.0409917355	0.0114049587	0	0	0	0	0	0	1	-360	360;
	91	96	0.0161983471	0.0148760331	0	0	0	0	0	0	1	-360	360;
	96	97	0.0161983471	0.0148760331	0	0	0	0	0	0	1	-360	360;
	97	98	0.0154214876	0.0100826446	0	0	0	0	0	0	1	-360	360;
	98	99	0.00616528926	0.0262809917	0	0	0	0	0	0	1	-360	360;
	1	100	0.00516528926	0.00219008264	0	0	0	0	0	0	1	-360	360;
	100	101	0.0124049587	0.019338843	0	0	0	0	0	0	1	-360	360;
	101	102	0.0111322314	0.00733884298	0	0	0	0	0	0	1	-360	360;
	102	103	0.0190661157	0.00994214876	0	0	0	0	0	0	1	-360	360;
	103	104	0.0369421488	0.0132892562	0	0	0	0	0	0	1	-360	360;
	104	105	0.0134876033	0.00485950413	0	0	0	0	0	0	1	-360	360;
	105	106	0.0272727273	0.00818181818	0	0	0	0	0	0	1	-360	360;
	106	107	0.012892562	0.00463636364	0	0	0	0	0	0	1	-360	360;
	107	108	0.0315619835	0.0113553719	0	0	0	0	0	0	1	-360	360;
	108	109	0.0134380165	0.00483471074	0	0	0	0	0	0	1	-360	360;
	109	110	0.0315619835	0.0113553719	0	0	0	0	0	0	1	-360	360;
	110	111	0.0202066116	0.00726446281	0	0	0	0	0	0	1	-360	360;
	110	112	0.0172561983	0.0062231405	0	0	0	0	0	0	1	-360	360;
	112	113	0.0190165289	0.00684297521	0	0	0	0	0	0	1	-360	360;
	100	114	0.0504297521	0.0181487603	0	0	0	0	0	0	1	-360	360;
	114	115	0.0154214876	0.0104958678	0	0	0	0	0	0	1	-360	360;
	115	116	0.0308429752	0.0203305785	0	0	0	0	0	0	1	-360	360;
	116	117	0.0334710744	0.0303305785	0	0	0	0	0	0	1	-360	360;
	117	118	0.0404132231	0.0361983471	0	0	0	0	0	0	1	-360	360;
	46	27	0.0434545455	0.0241735537	0	0	0	0	0	0	0	-360	360;
	17	27	0.0434545455	0.0240991736	0	0	0	0	0	0	0	-360	360;
	8	24	0.0353057851	0.0127190083	0	0	0	0	0	0	0	-360	360;
	54	43	0.0396694215	0.0142809917	0	0	0	0	0	0	0	-360	360;
	62	49	0.0297520661	0.0107107438	0	0	0	0	0	0	0	-360	360;
	37	62	0.047107438	0.0472727273	0	0	0	0	0	0	0	-360	360;
	9	40	0.0438016529	0.0276694215	0	0	0	0	0	0	0	-360	360;
	58	96	0.0327024793	0.0117768595	0	0	0	0	0	0	0	-360	360;
	73	91	0.0561983471	0.053553719	0	0	0	0	0	0	0	-360	360;
	88	75	0.0335702479	0.0120991736	0	0	0	0	0	0	0	-360	360;
	99	77	0.038231405	0.0138347107	0	0	0	0	0	0	0	-360	360;
	108	83	0.0538016529	0.019338843	0	0	0	0	0	0	0	-360	360;
	105	86	0.0671487603	0.0241735537	0	0	0	0	0	0	0	-360	360;
	110	118	0.0585867769	0.0210991736	0	0	0	0	0	0	0	-360	360;
	25	35	0.041322314	0.041322314	0	0	0	0	0	0	0	-360	360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	3	0	20	0;
];
