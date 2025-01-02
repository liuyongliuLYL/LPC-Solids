/*
	I_i表示C_i+1/2界面
	I0  I1  I2  I3  I4	     In-1 In	单元边界I
	|---|---|---|---|.....|---|---|
	 C1  C2  C3  C4        Cn-1 Cn		单元中心C
*/
#include "../data/use_data.h"
DATATYPE reconstruction_1order_upwind(DATATYPE a1);
DATATYPE reconstruction_2order_NND(bool flag, DATATYPE U1, DATATYPE U2, DATATYPE U3);
DATATYPE reconstruction_3order_MUSCL(bool flag, DATATYPE U1, DATATYPE U2, DATATYPE U3);

void riemann_ee_HLLC_1(Interface& II, Interface& reL_I, Interface& reR_I, int index);
void riemann_ee_HLL_1(Interface& II, Interface& reL_I, Interface& reR_I, int index);
void riemann_ee_HLL_HLLC(Interface& II, Interface& reL_I, Interface& reR_I, int index);
void riemann_acoustic_1(Interface& II, Interface& reL_I, Interface& reR_I, int index);

inline void reconstruction1_1upwind_1(Interface& reL_I, Interface& reR_I, int i)
{
	//最简单的一阶迎风重构
	reL_I.velocity[i] = reconstruction_1order_upwind(C_1.velocity[i]);
	reR_I.velocity[i] = reconstruction_1order_upwind(C_1.velocity[i + 1]);
	reL_I.density[i] = reconstruction_1order_upwind(C_1.density[i]);
	reR_I.density[i] = reconstruction_1order_upwind(C_1.density[i + 1]);
	reL_I.pressure[i] = reconstruction_1order_upwind(C_1.pressure[i]);
	reR_I.pressure[i] = reconstruction_1order_upwind(C_1.pressure[i + 1]);
	reL_I.s_11[i] = reconstruction_1order_upwind(C_1.s_11[i]);
	reR_I.s_11[i] = reconstruction_1order_upwind(C_1.s_11[i + 1]);
	//if (abs(reL_I.s_11[i]) > 2.0 * M.Y0[i] / 3.0) reL_I.s_11[i] = reL_I.s_11[i] / abs(reL_I.s_11[i]) * 2.0 * M.Y0[i] / 3.0;
	//if (abs(reR_I.s_11[i]) > 2.0 * M.Y0[i] / 3.0) reR_I.s_11[i] = reR_I.s_11[i] / abs(reR_I.s_11[i]) * 2.0 * M.Y0[i] / 3.0;
	//reL_I.sigma_11[i] = reL_I.s_11[i] - reL_I.pressure[i]; //柯西应力 = 偏应力 - 流体静水压
	//reR_I.sigma_11[i] = reR_I.s_11[i] - reR_I.pressure[i];

	reL_I.sigma_11[i] = reconstruction_1order_upwind(C_1.sigma_11[i]);
	reR_I.sigma_11[i] = reconstruction_1order_upwind(C_1.sigma_11[i + 1]);

	reL_I.sigma_22[i] = reconstruction_1order_upwind(C_1.sigma_22[i]);
	reR_I.sigma_22[i] = reconstruction_1order_upwind(C_1.sigma_22[i + 1]);
	reL_I.big_omega_1[i] = reconstruction_1order_upwind(C_1.big_omega_1[i]);
	reR_I.big_omega_1[i] = reconstruction_1order_upwind(C_1.big_omega_1[i + 1]);
	reL_I.big_phi_1[i] = reconstruction_1order_upwind(C_1.big_phi_1[i]);
	reR_I.big_phi_1[i] = reconstruction_1order_upwind(C_1.big_phi_1[i + 1]);

	reL_I.i_energy[i] = reconstruction_1order_upwind(C_1.i_energy[i]);
	reR_I.i_energy[i] = reconstruction_1order_upwind(C_1.i_energy[i + 1]);
	reL_I.t_energy[i] = reL_I.i_energy[i] + 0.5 * pow(reL_I.velocity[i], 2.0);
	reR_I.t_energy[i] = reR_I.i_energy[i] + 0.5 * pow(reR_I.velocity[i], 2.0);
}
inline void reconstruction2_2NND_1(Interface& reL_I, Interface& reR_I, int i)
{
	reL_I.velocity[i] = reconstruction_2order_NND(1, C_1.velocity[i - 1], C_1.velocity[i], C_1.velocity[i + 1]);
	reR_I.velocity[i] = reconstruction_2order_NND(0, C_1.velocity[i], C_1.velocity[i + 1], C_1.velocity[i + 2]);
	reL_I.density[i] = reconstruction_2order_NND(1, C_1.density[i - 1], C_1.density[i], C_1.density[i + 1]);
	reR_I.density[i] = reconstruction_2order_NND(0, C_1.density[i], C_1.density[i + 1], C_1.density[i + 2]);
	reL_I.pressure[i] = reconstruction_2order_NND(1, C_1.pressure[i - 1], C_1.pressure[i], C_1.pressure[i + 1]);
	reR_I.pressure[i] = reconstruction_2order_NND(0, C_1.pressure[i], C_1.pressure[i + 1], C_1.pressure[i + 2]);
	reL_I.s_11[i] = reconstruction_2order_NND(1, C_1.s_11[i - 1], C_1.s_11[i], C_1.s_11[i + 1]);
	reR_I.s_11[i] = reconstruction_2order_NND(0, C_1.s_11[i], C_1.s_11[i + 1], C_1.s_11[i + 2]);
	//if (abs(reL_I.s_11[i]) > 2.0 * M.Y0[i] / 3.0) reL_I.s_11[i] = reL_I.s_11[i] / abs(reL_I.s_11[i]) * 2.0 * M.Y0[i] / 3.0;
	//if (abs(reR_I.s_11[i]) > 2.0 * M.Y0[i] / 3.0) reR_I.s_11[i] = reR_I.s_11[i] / abs(reR_I.s_11[i]) * 2.0 * M.Y0[i] / 3.0;
	//reL_I.sigma_11[i] = reL_I.s_11[i] - reL_I.pressure[i]; //柯西应力 = 偏应力 - 流体静水压
	//reR_I.sigma_11[i] = reR_I.s_11[i] - reR_I.pressure[i];

	reL_I.sigma_11[i] = reconstruction_2order_NND(1, C_1.sigma_11[i - 1], C_1.sigma_11[i], C_1.sigma_11[i + 1]);
	reR_I.sigma_11[i] = reconstruction_2order_NND(0, C_1.sigma_11[i], C_1.sigma_11[i + 1], C_1.sigma_11[i + 2]);


	reL_I.sigma_22[i] = reconstruction_2order_NND(1, C_1.sigma_22[i - 1], C_1.sigma_22[i], C_1.sigma_22[i + 1]);
	reR_I.sigma_22[i] = reconstruction_2order_NND(0, C_1.sigma_22[i], C_1.sigma_22[i + 1], C_1.sigma_22[i + 2]);
	reL_I.big_omega_1[i] = reconstruction_2order_NND(1, C_1.big_omega_1[i - 1], C_1.big_omega_1[i], C_1.big_omega_1[i + 1]);
	reR_I.big_omega_1[i] = reconstruction_2order_NND(0, C_1.big_omega_1[i], C_1.big_omega_1[i + 1], C_1.big_omega_1[i + 2]);
	reL_I.big_phi_1[i] = reconstruction_2order_NND(1, C_1.big_phi_1[i - 1], C_1.big_phi_1[i], C_1.big_phi_1[i + 1]);
	reR_I.big_phi_1[i] = reconstruction_2order_NND(0, C_1.big_phi_1[i], C_1.big_phi_1[i + 1], C_1.big_phi_1[i + 2]);

	reL_I.i_energy[i] = reconstruction_2order_NND(1, C_1.i_energy[i - 1], C_1.i_energy[i], C_1.i_energy[i + 1]);
	reR_I.i_energy[i] = reconstruction_2order_NND(0, C_1.i_energy[i], C_1.i_energy[i + 1], C_1.i_energy[i + 2]);
	reL_I.t_energy[i] = reL_I.i_energy[i] + 0.5 * pow(reL_I.velocity[i], 2.0);
	reR_I.t_energy[i] = reR_I.i_energy[i] + 0.5 * pow(reR_I.velocity[i], 2.0);
}

void cell_interface_1(Cell_Value& C_1, Interface& I_1)
{
	Interface reL_I_1;	reL_I_1._malloc(para.cells + 1);
	Interface reR_I_1;	reR_I_1._malloc(para.cells + 1);

	if (para.mode_space == 1) {
#pragma omp parallel for num_threads(para.n_threads)
		for (int i = 1; i < para.cells; ++i) {
			reconstruction1_1upwind_1(reL_I_1, reR_I_1, i);
		}
	}
	else if (para.mode_space == 2)
	{
#pragma omp parallel for num_threads(para.n_threads)
		for (int i = 2; i < para.cells - 1; ++i) {
			reconstruction2_2NND_1(reL_I_1, reR_I_1, i);
		}
		reconstruction1_1upwind_1(reL_I_1, reR_I_1, 1);
		reconstruction1_1upwind_1(reL_I_1, reR_I_1, para.cells - 1);
	}


	M_1.copy_i_to_j(1, 0);
	C_1.copy_i_to_j(1, 0);
	reconstruction1_1upwind_1(reL_I_1, reR_I_1, 0);


	M_1.copy_i_to_j(para.cells, para.cells + 1);
	C_1.copy_i_to_j(para.cells, para.cells + 1);
	reconstruction1_1upwind_1(reL_I_1, reR_I_1, para.cells);

	if (para.mode_riemann == 0) {//LF
#pragma omp parallel for num_threads(para.n_threads)
		for (int i = 0; i <= para.cells; ++i) {
			I_1.sigma_11[i] = (reR_I_1.sigma_11[i] + reL_I_1.sigma_11[i]) / 2.0;
			I_1.velocity[i] = (reR_I_1.velocity[i] + reL_I_1.velocity[i]) / 2.0;
			I_1.sigma_velocity[i] = (reR_I_1.velocity[i] * reR_I_1.sigma_11[i] + reL_I_1.velocity[i] * reL_I_1.sigma_11[i]) / 2.0;
		}
	}
	else {
#pragma omp parallel for num_threads(para.n_threads)
		for (int i = 0; i <= para.cells; ++i) {
			riemann_acoustic_1(I_1, reL_I_1, reR_I_1, i);//传入I_1  
			//riemann_ee_HLL_1(I_1, reL_I_1, reR_I_1, i);
		}
	}
}
