#include "../data/use_data.h"
#include <algorithm>
using std::max;
using std::min;
DATATYPE sound_velocity(DATATYPE d, DATATYPE p, DATATYPE dst, DATATYPE e, int i);

void Wave_Speed_Estimates(DATATYPE& sel, DATATYPE& ser, Interface& reL_I, Interface& reR_I, int index)
{
	DATATYPE dl = reL_I.density[index]; DATATYPE ul = reL_I.velocity[index]; DATATYPE pl = reL_I.pressure[index]; DATATYPE dstl = reL_I.s_11[index]; DATATYPE sigmaL = reL_I.sigma_11[index]; DATATYPE el = reL_I.i_energy[index]; DATATYPE El = reL_I.t_energy[index];
	DATATYPE dr = reR_I.density[index]; DATATYPE ur = reR_I.velocity[index]; DATATYPE pr = reR_I.pressure[index]; DATATYPE dstr = reR_I.s_11[index]; DATATYPE sigmaR = reR_I.sigma_11[index]; DATATYPE er = reR_I.i_energy[index]; DATATYPE Er = reR_I.t_energy[index];

	DATATYPE cel = sound_velocity(dl, pl, dstl, el, index);
	DATATYPE cer = sound_velocity(dr, pr, dstr, er, index + 1);

	//DATATYPE sel = min(ul - cel, ur - cer);
	//DATATYPE ser = max(ul + cel, ur + cer);
	//sel = min(ul - cel, ur - cer);
	//ser = max(ul + cel, ur + cer);

	sel = ul - cel;
	ser = ur + cer;

	//sel *= 2;	ser *= 2;
}

void riemann_ee_HLLC1(Interface& II, Interface& reL_I, Interface& reR_I, int index)
{
	DATATYPE dl = reL_I.density[index]; DATATYPE ul = reL_I.velocity[index]; DATATYPE pl = reL_I.pressure[index]; DATATYPE dstl = reL_I.s_11[index]; DATATYPE sigmaL = reL_I.sigma_11[index]; DATATYPE el = reL_I.i_energy[index]; DATATYPE El = reL_I.t_energy[index];
	DATATYPE dr = reR_I.density[index]; DATATYPE ur = reR_I.velocity[index]; DATATYPE pr = reR_I.pressure[index]; DATATYPE dstr = reR_I.s_11[index]; DATATYPE sigmaR = reR_I.sigma_11[index]; DATATYPE er = reR_I.i_energy[index]; DATATYPE Er = reR_I.t_energy[index];

	DATATYPE sel; DATATYPE ser;
	Wave_Speed_Estimates(sel, ser, reL_I, reR_I, index);

	II.sigma_11[index] = (1.0 / (sel * dl - ul * dl - ser * dr + ur * dr)) *
		(-sel * ser * ul * dl * dr + ser * ul * ul * dl * dr +
			sel * ser * ur * dl * dr + sel * ul * ur * dl * dr -
			ser * ul * ur * dl * dr - ul * ul * ur * dl * dr -
			sel * ur * ur * dl * dr + ul * ur * ur * dl * dr -
			ser * dr * sigmaL + ur * dr * sigmaL +
			sel * dl * sigmaR - ul * dl * sigmaR);
	II.velocity[index] = (sel * ul * dl - ul * ul * dl - ser * ur * dr +
		ur * ur * dr + sigmaL - sigmaR) / (
			sel * dl - ul * dl - ser * dr + ur * dr);
	II.sigma_velocity[index] = II.velocity[index] * II.sigma_11[index];
}

void riemann_ee_HLLC(Interface& II, Interface& reL_I, Interface& reR_I, int index, int fix)
{
	DATATYPE dl = reL_I.density[index]; DATATYPE ul = reL_I.velocity[index]; DATATYPE pl = reL_I.pressure[index]; DATATYPE dstl = reL_I.s_11[index]; DATATYPE sigmaL = reL_I.sigma_11[index]; DATATYPE el = reL_I.i_energy[index]; DATATYPE El = reL_I.t_energy[index];
	DATATYPE dr = reR_I.density[index]; DATATYPE ur = reR_I.velocity[index]; DATATYPE pr = reR_I.pressure[index]; DATATYPE dstr = reR_I.s_11[index]; DATATYPE sigmaR = reR_I.sigma_11[index]; DATATYPE er = reR_I.i_energy[index]; DATATYPE Er = reR_I.t_energy[index];

	DATATYPE sel; DATATYPE ser;
	Wave_Speed_Estimates(sel, ser, reL_I, reR_I, index);

	II.velocity[index] = (sigmaR - sigmaL + dl * ul * (ul - sel) + dr * ur * (ser - ur))
		/ (dl * (ul - sel) + dr * (ser - ur));
	II.sigma_11[index] = (sigmaL * dr * (ser - ur) + sigmaR * dl * (ul - sel) + dl * (ul - sel) * dr * (ser - ur) * (ur - ul))
		/ (dr * (ser - ur) + dl * (ul - sel));


	II.sigma_22[index] = (reL_I.sigma_22[index] * dr * (ser - ur) + reR_I.sigma_22[index] * dl * (ul - sel) + dl * (ul - sel) * dr * (ser - ur) * (ur - ul))
		/ (dr * (ser - ur) + dl * (ul - sel));

	II.sigma_velocity[index] = II.velocity[index] * II.sigma_11[index];

	II.density_L[index] = dl * (sel - ul) / (sel - II.velocity[index]);
	II.density_R[index] = dr * (ser - ur) / (ser - II.velocity[index]);

}



void riemann_ee_HLL(Interface& II, Interface& reL_I, Interface& reR_I, int index)
{
	DATATYPE dl = reL_I.density[index]; DATATYPE ul = reL_I.velocity[index]; DATATYPE pl = reL_I.pressure[index]; DATATYPE dstl = reL_I.s_11[index]; DATATYPE sigmaL = reL_I.sigma_11[index]; DATATYPE el = reL_I.i_energy[index]; DATATYPE El = reL_I.t_energy[index];
	DATATYPE dr = reR_I.density[index]; DATATYPE ur = reR_I.velocity[index]; DATATYPE pr = reR_I.pressure[index]; DATATYPE dstr = reR_I.s_11[index]; DATATYPE sigmaR = reR_I.sigma_11[index]; DATATYPE er = reR_I.i_energy[index]; DATATYPE Er = reR_I.t_energy[index];

	DATATYPE sel; DATATYPE ser;
	Wave_Speed_Estimates(sel, ser, reL_I, reR_I, index);

	DATATYPE alphaL = -dl * (sel - ul);
	DATATYPE alphaR = dr * (ser - ur);

	//II.velocity[index] = (alphaR * sel + alphaL * ser) / (alphaR + alphaL);//HLL
	II.velocity[index] = (sigmaR - sigmaL + dl * ul * (ul - sel) + dr * ur * (ser - ur))
		/ (dl * (ul - sel) + dr * (ser - ur));//HLLC

	II.sigma_11[index] = (alphaR * sigmaL + alphaL * sigmaR - alphaL * alphaR * (ul - ur)) / (alphaL + alphaR);

	II.sigma_velocity[index] = (alphaR * sigmaL * ul + alphaL * sigmaR * ur - alphaL * alphaR * (El - Er)) / (alphaL + alphaR);

}



