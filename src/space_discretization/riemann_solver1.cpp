#include "../data/use_data.h"
#include <algorithm>
using std::max;
using std::min;
DATATYPE sound_velocity(DATATYPE d, DATATYPE p, DATATYPE dst, DATATYPE e, int i);

//DATATYPE riemann_acoustic(DATATYPE pl,DATATYPE ul)

void Wave_Speed_Estimates_1(DATATYPE& sel, DATATYPE& ser, Interface& reL_I, Interface& reR_I, int index)
{
	DATATYPE dl = reL_I.density[index]; DATATYPE ul = reL_I.velocity[index]; DATATYPE pl = reL_I.pressure[index]; DATATYPE dstl = reL_I.s_11[index]; DATATYPE sigmaL = reL_I.sigma_11[index]; DATATYPE el = reL_I.i_energy[index]; DATATYPE El = reL_I.t_energy[index];
	DATATYPE dr = reR_I.density[index]; DATATYPE ur = reR_I.velocity[index]; DATATYPE pr = reR_I.pressure[index]; DATATYPE dstr = reR_I.s_11[index]; DATATYPE sigmaR = reR_I.sigma_11[index]; DATATYPE er = reR_I.i_energy[index]; DATATYPE Er = reR_I.t_energy[index];

	dl = C.density[index];	dr = C.density[index + 1];
	pl = C.pressure[index];	pr = C.pressure[index + 1];
	ul = C.velocity[index];	ur = C.velocity[index + 1];

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
void riemann_ee_HLLC_1(Interface& II, Interface& reL_I, Interface& reR_I, int index)
{
	DATATYPE dl = reL_I.density[index]; DATATYPE ul = reL_I.velocity[index]; DATATYPE pl = reL_I.pressure[index]; DATATYPE dstl = reL_I.s_11[index]; DATATYPE sigmaL = reL_I.sigma_11[index]; DATATYPE el = reL_I.i_energy[index]; DATATYPE El = reL_I.t_energy[index];
	DATATYPE dr = reR_I.density[index]; DATATYPE ur = reR_I.velocity[index]; DATATYPE pr = reR_I.pressure[index]; DATATYPE dstr = reR_I.s_11[index]; DATATYPE sigmaR = reR_I.sigma_11[index]; DATATYPE er = reR_I.i_energy[index]; DATATYPE Er = reR_I.t_energy[index];

	DATATYPE sel; DATATYPE ser;
	Wave_Speed_Estimates_1(sel, ser, reL_I, reR_I, index);

	dl = C.density[index];	dr = C.density[index + 1];

	II.velocity[index] = (sigmaR - sigmaL + dl * ul * (ul - sel) + dr * ur * (ser - ur))
		/ (dl * (ul - sel) + dr * (ser - ur));
	II.sigma_11[index] = (sigmaL * dr * (ser - ur) + sigmaR * dl * (ul - sel) + dl * (ul - sel) * dr * (ser - ur) * (ur - ul))
		/ (dr * (ser - ur) + dl * (ul - sel));
	II.sigma_velocity[index] = II.velocity[index] * II.sigma_11[index];

	II.sigma_velocity[index] = II.velocity[index] * II.sigma_11[index];

}

void riemann_ee_HLL_1(Interface& IIx, Interface& reL_I, Interface& reR_I, int index)
{
	//DATATYPE CL = sqrt(M.big_gamma0[index] * (-C.sigma_11[index]) * C.density[index]);
	//DATATYPE CR = sqrt(M.big_gamma0[index + 1] * (-C.sigma_11[index + 1]) * C.density[index + 1]);

	DATATYPE CL = C.density[index] * sound_velocity(C.density[index], C.pressure[index], 0, 0, index);
	DATATYPE CR = C.density[index + 1] * sound_velocity(C.density[index + 1], C.pressure[index + 1], 0, 0, index + 1);
	DATATYPE Cstar = sqrt(max(1.0 / C.density[index] * CL * CL, 1.0 / C.density[index + 1] * CR * CR) / max(1.0 / C.density[index], 1.0 / C.density[index + 1]));

	IIx.sigma_11[index] = -(((-reR_I.sigma_11[index]) + (-reL_I.sigma_11[index])) / (2)
		+ Cstar / (2) * (reL_I.velocity[index] - reR_I.velocity[index]));
	IIx.sigma_22[index] = -(((-reR_I.sigma_22[index]) + (-reL_I.sigma_22[index])) / (2)
		+ Cstar / 2 * (reL_I.velocity[index] - reR_I.velocity[index]));

	IIx.velocity[index] = (reL_I.velocity[index] + reR_I.velocity[index]) / (2)
		+ 1.0 / (2 * Cstar) * (reR_I.sigma_11[index] - reL_I.sigma_11[index]);

	IIx.sigma_velocity[index] = IIx.sigma_11[index] * IIx.velocity[index];

	IIx.big_omega_1[index] = (reL_I.big_omega_1[index] + reR_I.big_omega_1[index]) / (2)
		+ 1.0 / (2 * Cstar) * (reR_I.big_phi_1[index] - reL_I.big_phi_1[index]);
	IIx.big_phi_1[index] = -(((-reR_I.big_phi_1[index]) + (-reL_I.big_phi_1[index])) / (2)
		+ Cstar / (2) * (reL_I.big_omega_1[index] - reR_I.big_omega_1[index]));
	//bug fix
	if (M.Y0[index] == 0) IIx.big_phi_1[index] = 0;
	//IIx.big_omega_1[index] = (CL * reL_I.big_omega_1[index] + CR * reR_I.big_omega_1[index]) / (CL + CR);
	//IIx.big_phi_1[index] = -(CL * (-reR_I.big_phi_1[index]) + CR * (-reL_I.big_phi_1[index])) / (CL + CR);

	//if (index == 0 && para.alpha != 1) {
	//	IIx.sigma_11[index] = 0;
	//	IIx.sigma_22[index] = 0;
	//	IIx.velocity[index] = 0;
	//	IIx.sigma_velocity[index] = 0;
	//	IIx.big_omega_1[index] = 0;
	//	IIx.big_phi_1[index] = 0;
	//}

}

void riemann_acoustic_1(Interface& IIx, Interface& reL_I, Interface& reR_I, int index)
{
	//DATATYPE CL = sqrt(M.big_gamma0[index] * (-C.sigma_11[index]) * C.density[index]);
	//DATATYPE CR = sqrt(M.big_gamma0[index + 1] * (-C.sigma_11[index + 1]) * C.density[index + 1]);

	DATATYPE CL = C.density[index] * sound_velocity(C.density[index], C.pressure[index], 0, 0, index);
	DATATYPE CR = C.density[index + 1] * sound_velocity(C.density[index + 1], C.pressure[index + 1], 0, 0, index + 1);

	IIx.sigma_11[index] = -((CL * (-reR_I.sigma_11[index]) + CR * (-reL_I.sigma_11[index])) / (CL + CR)
		+ CL * CR / (CL + CR) * (reL_I.velocity[index] - reR_I.velocity[index]));
	IIx.sigma_22[index] = -((CL * (-reR_I.sigma_22[index]) + CR * (-reL_I.sigma_22[index])) / (CL + CR)
		+ CL * CR / (CL + CR) * (reL_I.velocity[index] - reR_I.velocity[index]));
	IIx.velocity[index] = (CL * reL_I.velocity[index] + CR * reR_I.velocity[index]) / (CL + CR)
		+ 1.0 / (CL + CR) * (reR_I.sigma_11[index] - reL_I.sigma_11[index]);
	IIx.sigma_velocity[index] = IIx.sigma_11[index] * IIx.velocity[index];

	//CL = 1; CR = 1;
	IIx.big_omega_1[index] = (CL * reL_I.big_omega_1[index] + CR * reR_I.big_omega_1[index]) / (CL + CR)
		+ 1.0 / (CL + CR) * (reR_I.big_phi_1[index] - reL_I.big_phi_1[index]);
	IIx.big_phi_1[index] = -((CL * (-reR_I.big_phi_1[index]) + CR * (-reL_I.big_phi_1[index])) / (CL + CR)
		+ CL * CR / (CL + CR) * (reL_I.big_omega_1[index] - reR_I.big_omega_1[index]));
	//bug fix
	if (M.Y0[index] == 0) IIx.big_phi_1[index] = 0;

	//
	if (para.rbc == 2.0 && (index == para.cells || index == 0)) {
		IIx.sigma_11[index] = 0;
		IIx.sigma_22[index] = 0;
		//IIx.velocity[index] = 0;
		IIx.sigma_velocity[index] = 0;
		//IIx.big_omega_1[index] = 0;
		IIx.big_phi_1[index] = 0;
	}
	//
	if (para.rbc == 1.0 && (index == para.cells || index == 0)) {
		//IIx.sigma_11[index] = 0;
		//IIx.sigma_22[index] = 0;
		IIx.velocity[index] = 0;
		IIx.sigma_velocity[index] = 0;
		IIx.big_omega_1[index] = 0;
		//IIx.big_phi_1[index] = 0;
	}

	//if ((index == 0) && para.alpha != 1 && para.offset == 0.0) {
	//	IIx.sigma_11[index] = 0;
	//	IIx.sigma_22[index] = 0;
	//	IIx.velocity[index] = 0;
	//	IIx.sigma_velocity[index] = 0;
	//	IIx.big_omega_1[index] = 0;
	//	IIx.big_phi_1[index] = 0;
	//}
}


void riemann_acoustic_1_2order(Interface& II, Interface& reL_I, Interface& reR_I, int index)
{
	//DATATYPE CL = C.density[index] * sound_velocity(C.density[index], C.pressure[index], 0, 0, index);
	//DATATYPE CR = C.density[index + 1] * sound_velocity(C.density[index + 1], C.pressure[index + 1], 0, 0, index + 1);

	//if (para.now_step == 1) {
	//	riemann_acoustic_1(II, reL_I, reR_I, index);
	//}
	//else {
	//	DATATYPE fai; //Van Leer limiter
	//	DATATYPE D1, D2, D3;//Dj+1/2  Dj+1  Dj
	//	DATATYPE MM;//mj+1 - mj
	//	if (para.alpha == 1)
	//		MM = M.density_ini[index] * (mesh.r_half[index] - mesh.r[index]) + M.density_ini[index + 1] * (mesh.r[index + 1] - mesh.r_half[index]);
	//	else if (para.alpha == 2)
	//		MM = M.density_ini[index] * (pow(mesh.r_half[index], 2.0) - pow(mesh.r[index], 2.0))
	//		+ M.density_ini[index + 1] * (pow(mesh.r[index + 1], 2.0) - pow(mesh.r_half[index], 2.0));
	//	else if (para.alpha == 3)
	//		MM = M.density_ini[index] * (pow(mesh.r_half[index], 3.0) - pow(mesh.r[index], 3.0)) / 3.0
	//		+ M.density_ini[index + 1] * (pow(mesh.r[index + 1], 3.0) - pow(mesh.r_half[index], 3.0)) / 3.0;
	//	D1 = (C.velocity[index + 1] - C.velocity[index]) / MM;
	//	D2 = (C.velocity[index + 1] - C.velocity[index]) / mesh.dM[index+1];
	//}
}
