#include "../data/use_data.h"
#include <iostream>
inline DATATYPE get_Hugoniot_PH(DATATYPE c0, DATATYPE s, DATATYPE v0, DATATYPE v) {
    return c0 * c0 * (v0 - v) / pow(v0 - s * (v0 - v), 2.0);
}
inline DATATYPE get_Hugoniot_eH(DATATYPE PH, DATATYPE v0, DATATYPE v) {
    return PH * (v0 - v) / 2.0;
}
inline DATATYPE big_gamma_empirical_formula(DATATYPE _big_gamma0, DATATYPE rho, DATATYPE rho0)
{
    DATATYPE ita = rho / rho0;
    return _big_gamma0 / ita;
}

DATATYPE Gruneisen_EOS(int flag, int i, Cell_Value& C, Material_Constant& M)
{
    DATATYPE Pr = get_Hugoniot_PH(M.c0[i], M.s[i], 1.0 / M.density_ini[i], 1.0 / C.density[i]);
    DATATYPE er = get_Hugoniot_eH(Pr, 1.0 / M.density_ini[i], 1.0 / C.density[i]);

    DATATYPE big_gamma = big_gamma_empirical_formula(M.big_gamma0[i], C.density[i], M.density_ini[i]);

    if (flag == 1) {
        return (C.pressure[i] - Pr) / (big_gamma * C.density[i]) + er;
    }
    else if (flag == 2) {
        return Pr + big_gamma * C.density[i] * (C.i_energy[i] - er);
    }
    else {
        abort();
        return -1;
    }
}

DATATYPE Gruneisen_EOS1(int flag, int i, Cell_Value& C, Material_Constant& M)
{
    DATATYPE rho = C.density[i];
    DATATYPE _rho0 = M.density_ini[i];
    DATATYPE _big_gamma0 = M.big_gamma0[i];
    DATATYPE _s = M.s[i];
    DATATYPE _a0 = M.c0[i];

    DATATYPE ita = rho / _rho0;
    DATATYPE f = (ita - 1.0) * (ita - _big_gamma0 * (ita - 1.0) / 2.0) / ((ita - _s * (ita - 1.0)) * (ita - _s * (ita - 1.0)));
    if (flag == 1) {
        DATATYPE p = C.pressure[i];
        return (p - _rho0 * _a0 * _a0 * f) / (_rho0 * _big_gamma0);
    }
    else if (flag == 2) {
        DATATYPE e = C.i_energy[i];
        return e * (_rho0 * _big_gamma0) + _rho0 * _a0 * _a0 * f;
    }
    else {
        abort();
        return -1;
    }
}

DATATYPE Perfect_Gas_EOS(int flag, int i, Cell_Value& C, Material_Constant& M)
{
    if (flag == 1) {
        return C.pressure[i] / ((M.big_gamma0[i] - 1.0) * C.density[i]);
    }
    else if (flag == 2) {
        return (M.big_gamma0[i] - 1.0) * C.density[i] * C.i_energy[i];
    }
    else if (flag == 3) {
        return 0;
    }
    else {
        abort();
        return -1;
    }
}


DATATYPE cal_eos(int flag, int index, Cell_Value& C, Material_Constant& M)
{
    //if (flag == 3) {
    //	system("pause");
    //}
    if (para.EOS_set == 1) {//Gruneisen_EOS
        return Gruneisen_EOS(flag, index, C, M);
    }
    else if (para.EOS_set == 2) {
        return Perfect_Gas_EOS(flag, index, C, M);
    }
    else {
        abort();
        return -1;
    }
}

