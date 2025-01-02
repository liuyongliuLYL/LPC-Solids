
#include "../data/use_data.h"
DATATYPE reconstruction_1order_upwind(DATATYPE a1);
DATATYPE reconstruction_3order_MUSCL(bool flag, DATATYPE U1, DATATYPE U2, DATATYPE U3);
DATATYPE reconstruction_2order_NND(bool flag, DATATYPE U1, DATATYPE U2, DATATYPE U3);
DATATYPE cal_eos(int flag, int index, Cell_Value& C, Material_Constant& M);

void riemann_ee_HLLC(Interface& II, Interface& reL_I, Interface& reR_I, int index, int fix);
void riemann_ee_HLL(Interface& II, Interface& reL_I, Interface& reR_I, int index);
void riemann_ee_HLL_1(Interface& II, Interface& reL_I, Interface& reR_I, int index);
void riemann_acoustic_1(Interface& II, Interface& reL_I, Interface& reR_I, int index);

inline void reconstruction1_1upwind(Interface& reL_I, Interface& reR_I, int i)
{
	reL_I.velocity[i] = reconstruction_1order_upwind(C.velocity[i]);
	reR_I.velocity[i] = reconstruction_1order_upwind(C.velocity[i + 1]);
	reL_I.density[i] = reconstruction_1order_upwind(C.density[i]);
	reR_I.density[i] = reconstruction_1order_upwind(C.density[i + 1]);
	reL_I.pressure[i] = reconstruction_1order_upwind(C.pressure[i]);
	reR_I.pressure[i] = reconstruction_1order_upwind(C.pressure[i + 1]);
	reL_I.s_11[i] = reconstruction_1order_upwind(C.s_11[i]);
	reR_I.s_11[i] = reconstruction_1order_upwind(C.s_11[i + 1]);
	if (abs(reL_I.s_11[i]) > 2.0 * M.Y0[i] / 3.0) reL_I.s_11[i] = reL_I.s_11[i] / abs(reL_I.s_11[i]) * 2.0 * M.Y0[i] / 3.0;
	if (abs(reR_I.s_11[i]) > 2.0 * M.Y0[i] / 3.0) reR_I.s_11[i] = reR_I.s_11[i] / abs(reR_I.s_11[i]) * 2.0 * M.Y0[i] / 3.0;
	reL_I.sigma_11[i] = reL_I.s_11[i] - reL_I.pressure[i];
	reR_I.sigma_11[i] = reR_I.s_11[i] - reR_I.pressure[i];
	reL_I.sigma_11[i] = reconstruction_1order_upwind(C.sigma_11[i]);
	reR_I.sigma_11[i] = reconstruction_1order_upwind(C.sigma_11[i + 1]);
	reL_I.sigma_22[i] = reconstruction_1order_upwind(C.sigma_22[i]);
	reR_I.sigma_22[i] = reconstruction_1order_upwind(C.sigma_22[i + 1]);

	reL_I.i_energy[i] = reconstruction_1order_upwind(C.i_energy[i]);
	reR_I.i_energy[i] = reconstruction_1order_upwind(C.i_energy[i + 1]);
	reL_I.t_energy[i] = reL_I.i_energy[i] + 0.5 * pow(reL_I.velocity[i], 2.0);
	reR_I.t_energy[i] = reR_I.i_energy[i] + 0.5 * pow(reR_I.velocity[i], 2.0);
}
inline void reconstruction2_2NND(Interface& reL_I, Interface& reR_I, int i)
{
	reL_I.velocity[i] = reconstruction_2order_NND(1, C.velocity[i - 1], C.velocity[i], C.velocity[i + 1]);
	reR_I.velocity[i] = reconstruction_2order_NND(0, C.velocity[i], C.velocity[i + 1], C.velocity[i + 2]);
	reL_I.density[i] = reconstruction_2order_NND(1, C.density[i - 1], C.density[i], C.density[i + 1]);
	reR_I.density[i] = reconstruction_2order_NND(0, C.density[i], C.density[i + 1], C.density[i + 2]);
	reL_I.pressure[i] = reconstruction_2order_NND(1, C.pressure[i - 1], C.pressure[i], C.pressure[i + 1]);
	reR_I.pressure[i] = reconstruction_2order_NND(0, C.pressure[i], C.pressure[i + 1], C.pressure[i + 2]);
	reL_I.s_11[i] = reconstruction_2order_NND(1, C.s_11[i - 1], C.s_11[i], C.s_11[i + 1]);
	reR_I.s_11[i] = reconstruction_2order_NND(0, C.s_11[i], C.s_11[i + 1], C.s_11[i + 2]);
	if (abs(reL_I.s_11[i]) > 2.0 * M.Y0[i] / 3.0) reL_I.s_11[i] = reL_I.s_11[i] / abs(reL_I.s_11[i]) * 2.0 * M.Y0[i] / 3.0;
	if (abs(reR_I.s_11[i]) > 2.0 * M.Y0[i] / 3.0) reR_I.s_11[i] = reR_I.s_11[i] / abs(reR_I.s_11[i]) * 2.0 * M.Y0[i] / 3.0;
	reL_I.sigma_11[i] = reL_I.s_11[i] - reL_I.pressure[i];
	reR_I.sigma_11[i] = reR_I.s_11[i] - reR_I.pressure[i];
	reL_I.sigma_11[i] = reconstruction_2order_NND(1, C.sigma_11[i - 1], C.sigma_11[i], C.sigma_11[i + 1]);
	reR_I.sigma_11[i] = reconstruction_2order_NND(0, C.sigma_11[i], C.sigma_11[i + 1], C.sigma_11[i + 2]);
	//reL_I.s_11[i] = reL_I.sigma_11[i] + reL_I.pressure[i];
	//reR_I.s_11[i] = reR_I.sigma_11[i] + reR_I.pressure[i];
	//if (abs(reL_I.s_11[i]) > 2.0 * M.Y0[i] / 3.0) reL_I.s_11[i] = reL_I.s_11[i] / abs(reL_I.s_11[i]) * 2.0 * M.Y0[i] / 3.0;
	//if (abs(reR_I.s_11[i]) > 2.0 * M.Y0[i] / 3.0) reR_I.s_11[i] = reR_I.s_11[i] / abs(reR_I.s_11[i]) * 2.0 * M.Y0[i] / 3.0;

	reL_I.sigma_22[i] = reconstruction_2order_NND(1, C.sigma_22[i - 1], C.sigma_22[i], C.sigma_22[i + 1]);
	reR_I.sigma_22[i] = reconstruction_2order_NND(0, C.sigma_22[i], C.sigma_22[i + 1], C.sigma_22[i + 2]);

	reL_I.i_energy[i] = reconstruction_2order_NND(1, C.i_energy[i - 1], C.i_energy[i], C.i_energy[i + 1]);
	reR_I.i_energy[i] = reconstruction_2order_NND(0, C.i_energy[i], C.i_energy[i + 1], C.i_energy[i + 2]);
	reL_I.t_energy[i] = reL_I.i_energy[i] + 0.5 * pow(reL_I.velocity[i], 2.0);
	reR_I.t_energy[i] = reR_I.i_energy[i] + 0.5 * pow(reR_I.velocity[i], 2.0);
}

void cell_interface(Cell_Value& C, Interface& II)
{
	Interface reL_I;	reL_I._malloc(para.cells + 1);
	Interface reR_I;	reR_I._malloc(para.cells + 1);
	if (para.mode_space == 1) {
#pragma omp parallel for num_threads(para.n_threads)
		for (int i = 1; i < para.cells; ++i) {
			reconstruction1_1upwind(reL_I, reR_I, i);
		}
	}
	else if (para.mode_space == 2) {
#pragma omp parallel for num_threads(para.n_threads)
		for (int i = 2; i < para.cells - 1; ++i) {
			reconstruction2_2NND(reL_I, reR_I, i);
		}
		reconstruction1_1upwind(reL_I, reR_I, 1);
		reconstruction1_1upwind(reL_I, reR_I, para.cells - 1);
	}


	M.copy_i_to_j(1, 0);
	C.copy_i_to_j(1, 0);
	reconstruction1_1upwind(reL_I, reR_I, 0);
	if (para.lbc == 0.0);
	else if (para.lbc == 1.0) reL_I.velocity[0] *= (-1.0);
	else if (para.lbc == 2.0) {
		reL_I.sigma_11[0] *= (-1.0);
		reL_I.sigma_22[0] *= (-1.0);
	}
	//else if (para.lbc == 3.0) reL_I.sigma_11[0] = 2.0 * (-loading(para.now_time)) - reR_I.sigma_11[0];
	//else if (para.lbc == 4.0) reL_I.velocity[0] = 2.0 * (loading(para.now_time)) - reR_I.velocity[0];

	M.copy_i_to_j(para.cells, para.cells + 1);
	C.copy_i_to_j(para.cells, para.cells + 1);
	reconstruction1_1upwind(reL_I, reR_I, para.cells);
	if (para.rbc == 0.0);
	else if (para.rbc == 1.0) reR_I.velocity[para.cells] *= (-1.0);
	else if (para.rbc == 2.0) {
		reR_I.sigma_11[para.cells] *= (-1.0);
		reR_I.sigma_22[para.cells] *= (-1.0);
	}
	//else if (para.rbc == 3.0) reR_I.sigma_11[para.cells] = 2.0 * (-loading(para.now_time)) - reL_I.sigma_11[para.cells];
	//else if (para.rbc == 4.0) reR_I.velocity[para.cells] = 2.0 * (loading(para.now_time)) - reL_I.velocity[para.cells];

	if (para.mode_riemann == 0) {
#pragma omp parallel for num_threads(para.n_threads)
		for (int i = 0; i <= para.cells; ++i) {
			II.sigma_11[i] = (reR_I.sigma_11[i] + reL_I.sigma_11[i]) / 2.0;
			II.velocity[i] = (reR_I.velocity[i] + reL_I.velocity[i]) / 2.0;
			II.sigma_velocity[i] = (reR_I.velocity[i] * reR_I.sigma_11[i] + reL_I.velocity[i] * reL_I.sigma_11[i]) / 2.0;
		}

	}
	if (para.mode_riemann == 1) {//HLLC
#pragma omp parallel for num_threads(para.n_threads)
		for (int i = 0; i <= para.cells; ++i) {
			riemann_ee_HLLC(II, reL_I, reR_I, i, para.fix);
		}
	}
	else if (para.mode_riemann == 2) {//HLL
#pragma omp parallel for num_threads(para.n_threads)
		for (int i = 0; i <= para.cells; ++i) {
			riemann_ee_HLL(II, reL_I, reR_I, i);
		}
	}
//	else if (para.mode_riemann == 3) {//
//#pragma omp parallel for num_threads(para.n_threads)
//		for (int i = 0; i <= para.cells; ++i) {
//			riemann_ee_HLL_HLLC(II, reL_I, reR_I, i);
//		}
//	}
	else if (para.mode_riemann == 4) {//
#pragma omp parallel for num_threads(para.n_threads)
		for (int i = 0; i <= para.cells; ++i) {
			riemann_acoustic_1(II, reL_I, reR_I, i);
			//riemann_ee_HLL_1(I, reL_I, reR_I, i);
		}
	}

}



