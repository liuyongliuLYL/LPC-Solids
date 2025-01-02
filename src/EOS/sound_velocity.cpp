#include "../data/use_data.h"

DATATYPE sound_velocity(DATATYPE d, DATATYPE p, DATATYPE dst, DATATYPE e, int i)
{
	if (para.EOS_set == 2) {
		return sqrt(M.big_gamma0[i] * p / d);
	}
	//Cheng, J.2015 A high-order cell-centered Lagrangian scheme for one-dimensional    3.1.1 （18）
	DATATYPE _a0 = M.c0[i]; DATATYPE _s = M.s[i];
	DATATYPE _mu = M.mu[i]; DATATYPE _big_gamma0 = M.big_gamma0[i]; DATATYPE _density0 = M.density_ini[i];
	DATATYPE ita = d / _density0;
	//asq
	DATATYPE asq = _a0 * _a0 * (ita + (_s - _big_gamma0) * (ita - 1.0)) / (pow((ita - _s * (ita - 1.0)), 3)) + _big_gamma0 * p / (_density0 * (ita * ita));
	if (asq <= 0) {
		asq = _a0 * _a0;
		//system("pause");
	}

	//return sqrt(abs(asq - _density0 / (d * d) * _big_gamma0 * dst + 4.0 / 3.0 * _mu / d));
	//return sqrt(asq + 4.0 / 3.0 * mu / d);
	return sqrt(asq);
}

DATATYPE sound_c(int i)
{
	return sqrt(M.big_gamma0[i] * C.pressure[i] / C.density[i]);
}



