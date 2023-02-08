/* *****************************************************************************
 * CERN HTS Modeling
 * Volumetric heat capacity Cv of silver.
 *
 * Author: Erik Schnaubelt
 *
 ******************************************************************************/

Function{
	temp() = {
10,//artificial entry to extend to below 15K
15,
20,
25,
30,
35,
40,
45,
50,
55,
60,
65,
70,
75,
77.5,
80,
85,
86,
300
};

	Ic_temp() = { // Cv (J/m^3/K)
368.34,//artificial entry to extend to below 15k
368.34,
280.63,
209.29,
141.46,
102.31,
66.02,
46.78,
33.77,
24.43,
17.94,
12.53,
8.05,
4.48,
3.07,
1.82,
0.07,
0,
0
	};

	Ic()   = ListAlt[temp(), Ic_temp()];

	Ic[]  = $1 < 10? Ic_temp(0): InterpolationLinear[$1]{ Ic() };
	dIc[] = dInterpolationLinear[$1]{ Ic() } * SquDyadicProduct[$1];
}
