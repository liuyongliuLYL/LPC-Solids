#include "../data/use_data.h"

DATATYPE reconstruction_1order_upwind(DATATYPE a1)
{
	return a1;
}

DATATYPE reconstruction_3order_MUSCL(bool flag, DATATYPE U1, DATATYPE U2, DATATYPE U3)
{
	DATATYPE epsilon = 1e-6;
	/*
	*		U1	U2 U3
	        I-1 I I+1
	*/
	if (flag) {
		DATATYPE s1 = (2.0 * (U2 - U1) * (U3 - U2) + epsilon) / ((U2 - U1) * (U2 - U1) + (U3 - U2) * (U3 - U2) + epsilon);
		return U2 + s1 / 4.0 * ((1.0 - s1 / 3.0) * (U2 - U1) + (1.0 + s1 / 3.0) * (U3 - U2));
	}
	/*
	*		U1	U2 U3
	        I  I+1 I+2
	*/
	else {
		DATATYPE s2 = (2.0 * (U2 - U1) * (U3 - U2) + epsilon) / ((U2 - U1) * (U2 - U1) + (U3 - U2) * (U3 - U2) + epsilon);
		return U2 - s2 / 4.0 * ((1.0 - s2 / 3.0) * (U3 - U2) + (1.0 + s2 / 3.0) * (U2 - U1));
	}
}

DATATYPE min_mod(DATATYPE a, DATATYPE b) {
	if ((a * b) < 0) return 0;
	else {
		if (abs(a) > abs(b)) return b; else return a;
	}
}
DATATYPE reconstruction_2order_NND(bool flag, DATATYPE U1, DATATYPE U2, DATATYPE U3)
{
	/*
	*		U1	U2 U3
	        I-1 I I+1
	*/
	if (flag) {
		return U2 + 0.5 * min_mod(U3 - U2, U2 - U1);
	}
	/*
	*		U1	U2 U3
	         I  I+1 I+2
	*/
	else {
		return U2 - 0.5 * min_mod(U2 - U1, U3 - U2);
	}
}
