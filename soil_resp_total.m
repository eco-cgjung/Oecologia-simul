%observed total soil respiration
%days	UC UW CC CW			
function y=soil_resp_total(No)				
y=[				
10	0.483666667	0.4731	0.34534	0.58448
92	0.853633333	1.42555	0.936616667	0.904366667
106	1.778875	1.54942	1.4213	2.286025
128	3.342133333	3.88725	2.91244	3.751733333
141	3.3458	4.38648	3.431366667	3.95392
162	5.123666667	6.044684	5.0081	5.950233333
175	5.72278	5.527416667	5.489533333	6.44682
200	4.191166667	4.651966667	3.980383333	4.59906
206	4.756055	4.258083333	4.55442	4.5481
208	3.134533333	3.61805	2.2064	2.518725
212	4.61062	4.345116667	3.4232	3.64402
228	5.04822	5.0154	4.0871	3.92166
235	3.774466667	4.49755	3.4696	4.193525
249	3.37645	3.509	3.24566	3.691825
268	2.39135	2.858366667	2.268333333	3.3113
289	1.283033333	1.348275	1.1604	1.764725
312	1.386075	2.006825	2.43426	2.848868
324	0.715604	0.810825	0.84844	0.935412
346	1.185333333	1.283	1.2255	1.772
354	1	1.1025	0.615	0.908
410	0.703666667	0.72225	0.795	0.96925
436	1.20825	1.265	0.83325	0.81
466	1.55	1.5375	1.20775	1.2984
508	3.574	4.2175	3.538333333	4.153333333
534	5.125	5.435	4.404	6.326666667
557	4.586	4.742	4.034	4.312
592	3.64	4.054	3.048333333	3.911666667
634	3.054	4.085	2.9475	3.553333333
662	1.625333333	1.8	1.4238	1.686666667
691	0.977	1.29	0.92325	1.025333333
724	0.5492	0.5335	0.55525	0.466166667
759	0.3984	0.3394	0.3598	0.359583333
771	0.598	0.5306	0.416	0.492875
802	0.6582	0.8628	0.7665	0.84825
823	1.53	1.832	1.359833333	1.3138
860	3.67	3.91	3.250833333	3.736666667
897	5.001	4.718	4.405	4.942
928	5.113333333	5.452	4.543333333	5.266
968	4.11	4.962	3.726666667	4.333333333
991	3.36	3.0325	2.096666667	2.7525
1024	1.4175	1.8	1.1612	1.7825
1068	0.442166667	0.52075	0.318026667	0.47108
1094	0.6265	0.632	0.6698	0.619
1121	0.436	0.5535	0.4915	0.5588
1147	0.91325	0.810166667	0.9484	0.8255
1200	1.58	1.885	1.646666667	1.9375
1213	2.2	2.588333333	2.12	2.313333333
1243	3.7225	4.14	3.69	4.1
1267	5.646666667	5.525	4.488333333	5.173333333
1284	5.408333333	5.57	4.49	5.21
1309	4.881666667	5.283333333	4.221666667	4.33
1332	4.864	5.888	3.996	4.99
1362	3.48	4.034	3.2	2.76
1394	1.194	1.458333333	1.3275	1.473333333
1420	0.6622	0.68025	0.774666667	0.949
1459	0.4302	0.359666667	0.474666667	0.5665
1490	0.431333333	0.7415	0.51075	0.7576
1515	0.487	0.73475	0.628333333	0.717
1547	0.7768	0.7788	0.80925	0.9145
1577	1.874	1.906	1.515	1.59375
1596	2.746	2.866666667	2.553333333	3.023333333
1610	4.21	5.138	3.682	4.0125
1637	4.294	4.934	4.431666667	5.1025
1658	3.29	3.478	2.8125	2.8225
1673	4.578333333	4.628333333	4.056666667	3.924
1711	3.301666667	3.314	2.491666667	2.9
1741	2.375	2.9325	2.55	2.708333333
1775	0.871	0.835	0.7808	1.351
1817	0.294166667	0.551166667	0.277666667	0.455166667
1843	0.3525	0.4195	0.405666667	0.540083333
1880	0.487	0.73475	0.634	0.717
1891	0.6488	0.7174	0.9334	0.7754
1927	0.985	1.0474	1.078333333	1.625
1944	2.041666667	2.635	2.446666667	3.115
1974	1.816666667	2.7176	1.363	2.0416
1997	5.646666667	5.625	4.488333333	5.438
2023	5.842	6.655	5.155833333	6.606666667
2062	4.531666667	5.44	4.491666667	5.301666667
2103	2.8	3.646666667	2.809166667	3.49
2129	0.989833333	1.386666667	1.141166667	1.843333333
2178	0.315	0.439	0.396833333	0.6445
2208	0.33825	0.4195	0.3525	0.540083333
2251	0.539	0.5172	0.5572	0.629833333
2281	1.275	1.52	1.184333333	1.7
2311	1.611666667	1.993333333	1.9775	2.273333333
2345	6.403333333	7.252	5.19	5.569166667
2365	5.816666667	6.536	5.27	6.293333333
2390	5.575	6.414166667	5.092	5.983333333
2424	2.664	4.57	2.402666667	3.906
2460	2.396666667	2.943333333	2.353333333	3.4
2494	0.791333333	0.829833333	0.797333333	1.0974
2525	0.3806	0.5302	0.503	0.679
2555	0.3806	0.5302	0.503	0.679
2638	0.7954	0.893833333	1.0182	1.305
2683	2.344	2.678333333	1.984	2.673333333
2703	2.796666667	3.258333333	2.521666667	4.284
2731	4.926666667	5.361666667	4.158	5.03
2750	3.838333333	5.008333333	3.47	4.73
2780	4.51	5.49	3.613333333	4.498333333
2796	4.123333333	4.571666667	3.41	3.911666667
2856	0.987166667	1.395	0.6956	1.06975
2968	0.29	0.4252	0.553333333	0.548
2997	0.66	0.76	0.542	0.9924
3037	1.156	2.6725	1.585	2.548
3057	2.528	3.53	2.575	4.324
3079	3.633333333	5.286	3.103333333	5.936
3101	5.633333333	4.8725	4.521666667	5.551666667
3124	4.402	5.485	5.91	7.41
3184	2.646	3.273333333	2.326	2.686666667
3242	0.513333333	0.956666667	0.622	0.706666667
3276	0.275	0.328	0.345	0.53
3307	0.245	0.505	0.3648	0.495
3338	0.346666667	0.526666667	0.5125	0.706
3367	0.705	0.746666667	0.956666667	1.034
3403	1.46	1.73	1.6828	1.696666667
3431	3.105	3.535	3.041666667	4.895833333
3448	3.485	3.83	4.49375	4.73
3468	2.86875	3.705	2.533333333	2.7525
3489	1.591666667	2.602	1.355	1.6175
3510	3.28	5.536	3.8925	5.7625
3531	2.38125	2.673333333	2.8825	2.955
3551	1.895833333	2.521666667	2.075	2.19375
3576	1.45	1.426	1.258333333	1.354166667
3601	0.475	0.74	0.4675	0.74
3629	0.3675	0.4625	0.4672	0.436666667
3665	0.3354	0.489666667	0.472333333	0.4108
3700	0.45825	0.598	0.3562	0.88075
3723	0.899166667	1.1414	0.619666667	1.612
3762	2.192666667	3.310666667	2.769	3.967166667
3790	5.722166667	6.255166667	5.1688	5.304
3819	5.975666667	7.186833333	5.663666667	6.773
3853	3.343166667	4.173	3.4242	3.527333333
3886	3.6205	4.333333333	2.808	3.642166667
3916	1.692166667	1.911	1.365	2.088666667
3952	0.936	1.133166667	0.728	1.144
3972	0.8086	0.875333333	0.8008	0.966333333
4003	0.6448	0.712833333	0.5824	0.6656
4035	0.39	0.325	0.348333333	0.376
4077	0.39	0.33	0.3125	0.446666667
4105	1.656666667	2.35	1.5325	2.524
4133	1.36	1.845	1.4	2.795
4164	3.094	2.676	2.3425	3.848333333
4168	1.975	3.266	2.382	3.99
4190	4.3	4.823333333	4.27	4.975
4228	5.025	4.9275	4.458	4.946666667
4265	2.873333333	3.304	2.946666667	3.644
4297	1.63	1.616	1.445	1.918
4392	0.302	0.488333333	0.342	0.451666667
4429	0.3825	0.426	0.345	0.498333333
4479	0.81	1.365	1.506666667	1.8625
4516	2.326666667	1.726666667	2.682	3.235
4532	5.1475	4.214	5.432	4.498
4566	7.53	4.89	5.386	5.426666667
4656	1.456666667	1.5125	0.723333333	3.338
4685	1.048	0.9775	0.6675	1.49
4761	0.4672	0.47	2.003333333	0.5575
4789	0.456	0.4725	0.4	0.51
4832	0.95	1.385	0.748	1.703333333
4860	2.343333333	5.243333333	2.964	4.235
4906	4.404	4.821666667	5.582222222	6.006
4955	6.89	7.856666667	8.01	8.591666667
4992	4.41	3.942	4.55	3.674
5037	1.658	2.1125	1.69	1.558
5062	1.745	2.156	1.584	2.068
5090	0.7075	1.145	0.82	0.96
5137	0.886666667	0.7425	0.52	0.87
5160	1.096	1.85	1.0875	1.468333333
5191	2.471666667	2.364	1.281666667	2.996666667
5223	1.565	4.84	3.19	4.366666667
5270	6.025	9.015	5.525	9.8175
5309	5.622	9.295	8.0225	8.35
5347	4.405	5.928	5.57	4.131666667
5375	4.612	5.2825	4.424	4.786666667
5411	1.973333333	2.9925	2.864	3.008333333
5431	1.718	2.27	1.97	2.468333333
5472	0.635	0.865	0.855	0.9475
5501	0.626666667	0.454	0.426	0.788
5537	0.905	1.576666667	0.993333333	2.073333333
5565	0.905	1.576666667	0.993333333	2.073333333
5589	1.68	2.6475	2.752	3.841666667
5620	3.55	4.2	3.75	5.644
5655	4.728	6.248	5.44	5.2875
5676	7.0575	7.521666667	9.9	8.5225
5702	7.4575	7.11	6.092	10.8475
5738	5.755	3.2	4.725	5.866666667
5750	2.46	2.58	2.068	2.446
];