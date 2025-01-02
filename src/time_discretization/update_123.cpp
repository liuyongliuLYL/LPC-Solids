
#include "../data/use_data.h"
#include <algorithm>
#include <iostream>
using std::max;
using std::min;
#define Pi    3.14159265358979323846

void cell_interface(Cell_Value& C, Interface& II);
void cell_interface_1(Cell_Value& C_1, Interface& I_1);
DATATYPE sound_velocity(DATATYPE d, DATATYPE p, DATATYPE dst, DATATYPE e, int i);
DATATYPE cal_eos(int flag, int index, Cell_Value& C, Material_Constant& M);

void update_dt();
void update_r();
void update_firstorder();


int flag111 = 0;
int flag_g = 0;

Cell_Value C_next;
Mesh mesh_next;
Cell_Value C_1_next;
Mesh mesh_1_next;


void update_1euler_3()
{
	C_next = C; //C_next._malloc(para.cells + 2);
	mesh_next = mesh;
	cell_interface(C, II);

	if (para.fo) {
		C_1_next = C_1;
		mesh_1_next = mesh_1;
		cell_interface_1(C_1, I_1);
	}

	update_dt();
	update_r();

	if (!para.RTI) {
#pragma omp parallel for num_threads(para.n_threads)
		for (int i = 1; i <= para.cells; ++i) {
			////////////////////////////////base flow/////////////////////////////////////////////////////////////////////
			DATATYPE delta = mesh.dM[i];

			C_next.density[i] = 1.0 /
				(1.0 / C.density[i] + ((mesh.A[i] * II.velocity[i] - mesh.A[i - 1] * II.velocity[i - 1]) / delta) * para.dt);
			C_next.volume[i] = 1.0 / C_next.density[i];
			C_next.velocity[i] = (mesh.A[i] * II.sigma_11[i] - mesh.A[i - 1] * II.sigma_11[i - 1]) * para.dt / delta + C.velocity[i]
				//+para.dt * (-II.sigma_22[i] - II.sigma_22[i - 1]) * (mesh.A[i] - mesh.A[i - 1]) / (2.0 * delta);//S
				+ para.dt * (-C.sigma_22[i]) * (mesh.A[i] - mesh.A[i - 1]) / (delta);
			C_next.velocity[i] = C_next.velocity[i] + para.dt * para.g * flag_g;

			C_next.t_energy[i] = ((mesh.A[i] * II.sigma_velocity[i] - mesh.A[i - 1] * II.sigma_velocity[i - 1]) / delta) * para.dt + C.t_energy[i];
			C_next.t_energy[i] = C_next.t_energy[i] + para.dt * para.g * C.velocity[i] * flag_g;

			C_next.s_11[i] = C.s_11[i] + para.dt * (4.0 / 3) * M.mu[i] * (II.velocity[i] - II.velocity[i - 1]) / mesh.dR[i];
			C_next.s_22[i] = C.s_22[i] + para.dt * (-2.0 / 3) * M.mu[i] * (II.velocity[i] - II.velocity[i - 1]) / mesh.dR[i];
			if (para.alpha == 2.0) {
				C_next.s_11[i] = C_next.s_11[i] + para.dt * (-2.0 / 3) * M.mu[i] * C_next.velocity[i] / mesh_next.r[i];
				C_next.s_22[i] = C_next.s_22[i] + para.dt * (4.0 / 3) * M.mu[i] * C_next.velocity[i] / mesh_next.r[i];
			}
			else if (para.alpha == 3.0) {
				C_next.s_11[i] = C_next.s_11[i] + para.dt * (-4.0 / 3) * M.mu[i] * C_next.velocity[i] / mesh_next.r[i];
				C_next.s_22[i] = C_next.s_22[i] + para.dt * (2.0 / 3) * M.mu[i] * C_next.velocity[i] / mesh_next.r[i];
			}

			if (pow(C_next.s_11[i], 2.0) >= pow(M.Y0[i], 2.0) / 3.0) {
				if (C_next.s_11[i] > M.Y0[i] / sqrt(3.0)) C_next.s_11[i] = M.Y0[i] / sqrt(3.0);
				else C_next.s_11[i] = -M.Y0[i] / sqrt(3.0);

				if (C_next.s_22[i] > 0.0) C_next.s_22[i] = M.Y0[i] / sqrt(3.0);
				else C_next.s_22[i] = -M.Y0[i] / sqrt(3.0);
			}


			C_next.i_energy[i] = C_next.t_energy[i] - 0.5 * C_next.velocity[i] * C_next.velocity[i];
			C_next.pressure[i] = cal_eos(2, i, C_next, M);

			if (M.Y0[i] == 0) {
				C_next.sigma_11[i] = -C_next.pressure[i];
				C_next.sigma_22[i] = -C_next.pressure[i];
			}
			else {
				C_next.sigma_11[i] = C_next.s_11[i] - C_next.pressure[i];
				C_next.sigma_22[i] = C_next.s_22[i] - C_next.pressure[i];
			}
		}
	}
	else {
		C_next = C; //C_next._malloc(para.cells + 2);
		mesh_next = mesh;
	}

	if (para.fo) update_firstorder();

	C = C_next;			mesh = mesh_next;
	if (para.fo) { C_1 = C_1_next;		mesh_1 = mesh_1_next; }
}

void update_r()
{
	mesh_next.R_half[0] = mesh.R_half[0] + para.dt * II.velocity[0];
	if (para.fo) mesh_1_next.R_half[0] = mesh_1.R_half[0] + para.dt * I_1.velocity[0];

	for (int i = 1; i <= para.cells; ++i) {
		mesh_next.R_half[i] = mesh.R_half[i] + para.dt * II.velocity[i];		
		if (para.fo)mesh_1_next.R_half[i] = mesh_1.R_half[i] + para.dt * I_1.velocity[i];
		mesh_next.dR[i] = mesh_next.R_half[i] - mesh_next.R_half[i - 1];
		mesh_next.R[i] = mesh_next.R_half[i - 1] + mesh_next.dR[i] / 2.0;	
		if (para.alpha != 1.0) mesh_next.dm[i] = C.density[i] * mesh.dR[i];

		if (para.alpha == 3.0) {
			mesh.A[i] = (pow(mesh_next.R_half[i], 2) + mesh_next.R_half[i] * mesh.R_half[i] + pow(mesh.R_half[i], 2)) / 3.0;
			mesh.A[i - 1] = (pow(mesh_next.R_half[i - 1], 2) + mesh_next.R_half[i - 1] * mesh.R_half[i - 1] + pow(mesh.R_half[i - 1], 2)) / 3.0;

			if (para.fo)mesh_1.A[i] = (2.0 * mesh_next.R_half[i] * mesh_1_next.R_half[i] + mesh_next.R_half[i] * mesh_1.R_half[i] + mesh.R_half[i] * mesh_1_next.R_half[i] + 2.0 * mesh.R_half[i] * mesh_1.R_half[i]) / 3.0;
			if (para.fo)mesh_1.A[i - 1] = (2.0 * mesh_next.R_half[i - 1] * mesh_1_next.R_half[i - 1] + mesh_next.R_half[i - 1] * mesh_1.R_half[i - 1] + mesh.R_half[i - 1] * mesh_1_next.R_half[i - 1] + 2.0 * mesh.R_half[i - 1] * mesh_1.R_half[i - 1]) / 3.0;
		}
		else if (para.alpha == 2.0) {
			mesh.A[i] = (mesh_next.R_half[i] + mesh.R_half[i]) / 2.0;
			mesh.A[i - 1] = (mesh_next.R_half[i - 1] + mesh.R_half[i - 1]) / 2.0;
			if (para.fo)mesh_1.A[i] = (mesh_1_next.R_half[i] + mesh_1.R_half[i]) / 2.0;		mesh_1.A[i - 1] = (mesh_1_next.R_half[i - 1] + mesh_1.R_half[i - 1]) / 2.0;
		}
		else if (para.alpha == 1.0) {
			mesh.A[i] = 1.0;	mesh.A[i - 1] = 1.0;
			if (para.fo)mesh_1.A[i] = 0.0;	if (para.fo)mesh_1.A[i - 1] = 0.0;
		}
	}
	if (para.lbc == 0.0) mesh_next.dR[1] = mesh_next.dR[2];
	if (para.rbc == 0.0) mesh_next.dR[para.cells] = mesh_next.dR[para.cells - 1];
}

void update_firstorder()
{
#pragma omp parallel for num_threads(para.n_threads)
	for (int i = 1; i <= para.cells; ++i) {
		DATATYPE delta = mesh.dM[i];
		///////////////////////////////////first order/////////////////////////////////////////////////////////////////
		/*34*/
		int flag_V = 3;
		if (flag_V == 1) {
			if (para.alpha == 3.0) {
				C_1_next.big_omega_1[i] = (
					(mesh.R[i]) * C_1.big_omega_1[i] +
					para.omega * para.dt * (
						1.0 / C.density[i] * (-C_1.sigma_22[i]) / mesh.R[i] +
						(mesh.R_half[i] * mesh_1.R_half[i] + mesh.R_half[i - 1] * mesh_1.R_half[i - 1]) / 2.0 *
						(II.sigma_22[i] - II.sigma_22[i - 1]) / delta
						)
					+ para.dt * (
						(mesh.A[i] + mesh.A[i - 1]) / 2.0 *
						(I_1.big_phi_1[i] * mesh_next.R_half[i] - I_1.big_phi_1[i - 1] * mesh_next.R_half[i - 1]) / (delta)
						+3.0 * C_next.volume[i] * mesh_next.R[i] * C_1.big_phi_1[i]
						- C_next.velocity[i] * C_1.big_omega_1[i] 
						)
					) / mesh_next.R[i];
			}
			else if (para.alpha == 2.0) {
				C_1_next.big_omega_1[i] = (
					(mesh.R_half[i] + mesh.R_half[i - 1]) * 0.5 * C_1.big_omega_1[i] +
					para.omega * para.dt * (
						1.0 / C.density[i] * (-C_1.sigma_22[i]) / mesh_next.r[i] +
						(mesh_1.R_half[i] + mesh_1.R_half[i - 1]) / 2.0 *
						(II.sigma_22[i] - II.sigma_22[i - 1]) / delta
						)
					+ para.dt * (
						mesh.A[i] *
						(I_1.big_phi_1[i] * mesh.A[i] - I_1.big_phi_1[i - 1] * mesh.A[i - 1]) / (delta)
						+2.0 * C_next.volume[i] * C_1.big_phi_1[i]
						- C_next.velocity[i] * C_1.big_omega_1[i]
						)
					) / ((mesh_next.R_half[i] + mesh_next.R_half[i - 1]) * 0.5);
			}

			else if (para.alpha == 1.0)
				C_1_next.big_omega_1[i] = C_1.big_omega_1[i]
				- para.k * para.k * para.dt * (
					1.0 / C.density[i] * C_1.sigma_22[i] -
					(mesh_1.R_half[i] + mesh_1.R_half[i - 1]) / 2.0 * (II.sigma_22[i] - II.sigma_22[i - 1]) / (delta))
				+ para.dt * (I_1.big_phi_1[i] - I_1.big_phi_1[i - 1]) / (delta)
				;
		}
		else if (flag_V == 2) {
			if (para.alpha == 3.0) {
				C_1_next.big_omega_1[i] = (
					((pow(mesh.R_half[i], 2) + mesh.R_half[i] * mesh.R_half[i - 1] + pow(mesh.R_half[i], 2)) / 3.0) * C_1.big_omega_1[i] +
					para.omega * para.dt * (
						1.0 / C.density[i] * (-C_1.sigma_22[i]) +
						(mesh.A[i] * mesh_1.R_half[i] + mesh.A[i - 1] * mesh_1.R_half[i - 1]) / 2.0 *
						(II.sigma_22[i] - II.sigma_22[i - 1]) / delta
						)
					+ para.dt * (
						(mesh.A[i] * mesh_next.R_half[i] + mesh.A[i - 1] * mesh_next.R_half[i - 1]) / 2.0 *
						(I_1.big_phi_1[i] * mesh_next.R_half[i] - I_1.big_phi_1[i - 1] * mesh_next.R_half[i - 1]) / (delta)
						+3.0 * C_next.volume[i] * mesh_next.R[i] * C_1.big_phi_1[i]
						)
					) / ((pow(mesh_next.R_half[i], 2) + mesh_next.R_half[i] * mesh_next.R_half[i - 1] + pow(mesh_next.R_half[i], 2)) / 3.0);
			}
			else if (para.alpha == 2.0) {
				C_1_next.big_omega_1[i] = (
					((pow(mesh.R_half[i], 2) + mesh.R_half[i] * mesh.R_half[i - 1] + pow(mesh.R_half[i], 2)) / 3.0) * C_1.big_omega_1[i] +
					para.omega * para.dt * (
						1.0 / C.density[i] * (-C_1.sigma_22[i]) +
						(mesh_1.R_half[i] * mesh.R_half[i] + mesh_1.R_half[i - 1] * mesh.R_half[i - 1]) / 2.0 *
						(II.sigma_22[i] - II.sigma_22[i - 1]) / delta
						)
					+ para.dt * (
						pow((mesh.A[i] + mesh.A[i - 1]) / 2.0, 2.0) *
						(I_1.big_phi_1[i] * mesh.A[i] - I_1.big_phi_1[i - 1] * mesh.A[i - 1]) / (delta)
						+2.0 * C_next.volume[i] * C_1.big_phi_1[i] * (mesh.A[i] + mesh.A[i - 1]) / 2.0
						)
					) / ((pow(mesh_next.R_half[i], 2) + mesh_next.R_half[i] * mesh_next.R_half[i - 1] + pow(mesh_next.R_half[i], 2)) / 3.0);
			}
			else if (para.alpha == 1.0)
				C_1_next.big_omega_1[i] = C_1.big_omega_1[i]
				- para.k * para.k * para.dt * (
					1.0 / C.density[i] * C_1.sigma_22[i] -
					(mesh_1.R_half[i] + mesh_1.R_half[i - 1]) / 2.0 * (II.sigma_22[i] - II.sigma_22[i - 1]) / (delta))
				+ para.dt * (I_1.big_phi_1[i] - I_1.big_phi_1[i - 1]) / (delta)
				;
		}
		else if (flag_V == 3) {
			para.beta = 2;
			if (para.alpha == 1) para.beta = 0;
			REAL temp1 = 1, temp2 = 1; //r^β
			if (para.alpha != 1) {
				temp1 = ((pow(mesh.R_half[i], 2) + mesh.R_half[i] * mesh.R_half[i - 1] + pow(mesh.R_half[i], 2)) / 3.0);
				temp2 = ((pow(mesh_next.R_half[i], 2) + mesh_next.R_half[i] * mesh_next.R_half[i - 1] + pow(mesh_next.R_half[i], 2)) / 3.0);
			}
			C_1_next.big_omega_1[i] = (
				temp1
				* C_1.big_omega_1[i]
				+ para.omega * para.dt * (
					1.0 / C.density[i] * (-C_1.sigma_22[i]) +
					//mesh_1.R[i] * pow(mesh.R[i], para.alpha - 1) *//
					(mesh_1.R_half[i] * mesh.A[i] + mesh_1.R_half[i - 1] * mesh.A[i - 1]) / 2.0 *
					(II.sigma_22[i] - II.sigma_22[i - 1]) / delta
					)
				+ para.dt * (I_1.big_phi_1[i] * pow(mesh.R_half[i], para.alpha - 1 + para.beta) - I_1.big_phi_1[i - 1] * pow(mesh.R_half[i - 1], para.alpha - 1 + para.beta)) / (delta)
				) /
				temp2
				;
		}

		/*125*/
		{
			C_1_next.big_lambda_1[i] = C_1.big_lambda_1[i] + para.dt * C_1_next.big_omega_1[i];

			C_1_next.volume[i] =
				C_1.volume[i] - M_1.density_ini[i] / M.density_ini[i] * (1.0 / C_next.density[i] - 1.0 / C.density[i])
				//C_1.volume[i] + M_1.volume_ini[i] * M.density_ini[i] * (1.0 / C_next.density[i] - 1.0 / C.density[i])
				+ para.dt / delta * (mesh_1.A[i] * II.velocity[i] - mesh_1.A[i - 1] * II.velocity[i - 1])
				+ para.dt / delta * (mesh.A[i] * I_1.velocity[i] - mesh.A[i - 1] * I_1.velocity[i - 1])
				+ para.dt / delta * (mesh.A[i] * II.velocity[i] - mesh.A[i - 1] * II.velocity[i - 1]) * C_1_next.big_lambda_1[i]
				+ para.dt * C_1_next.big_omega_1[i] / C.density[i]
				;
			C_1_next.velocity[i] = (
				C_1.velocity[i] - M_1.density_ini[i] / M.density_ini[i] * (C_next.velocity[i] - C.velocity[i])
				//C_1.velocity[i] + M_1.volume_ini[i] * M.density_ini[i] * (C_next.velocity[i] - C.velocity[i])
				+ para.dt / delta * (mesh_1.A[i] * II.sigma_11[i] - mesh_1.A[i - 1] * II.sigma_11[i - 1])
				+ para.dt / delta * (mesh.A[i] * I_1.sigma_11[i] - mesh.A[i - 1] * I_1.sigma_11[i - 1])
				+ para.dt / delta * (mesh.A[i] * II.sigma_11[i] - mesh.A[i - 1] * II.sigma_11[i - 1]) * C_1_next.big_lambda_1[i]
				+ para.dt * (//Σ
					C_1.big_phi_1[i] / C.density[i]
					+ (-C_1.sigma_22[i]) * (mesh.A[i] - mesh.A[i - 1]) / (delta)
					//+ (-I_1.sigma_22[i] - I_1.sigma_22[i - 1]) * (mesh.A[i] - mesh.A[i - 1]) / (2.0 * delta)

					+(-C.sigma_22[i]) * (mesh_1.A[i] - mesh_1.A[i - 1]) / (delta)
					//+ (-II.sigma_22[i] - II.sigma_22[i - 1]) * (mesh_1.A[i] - mesh_1.A[i - 1]) / (2.0 * delta)
					+C_1_next.big_lambda_1[i] * (-II.sigma_22[i] - II.sigma_22[i - 1]) * (mesh.A[i] - mesh.A[i - 1]) / (2.0 * delta)
					)
				);

			C_1_next.t_energy[i] = (
				C_1.t_energy[i] - M_1.density_ini[i] / M.density_ini[i] * (C_next.t_energy[i] - C.t_energy[i])
				//C_1.t_energy[i] + M_1.volume_ini[i] * M.density_ini[i] * (C_next.t_energy[i] - C.t_energy[i])
				+ para.dt / delta * (mesh_1.A[i] * II.sigma_11[i] * II.velocity[i] - mesh_1.A[i - 1] * II.sigma_11[i - 1] * II.velocity[i - 1])
				+ para.dt / delta * (mesh.A[i] * (I_1.sigma_11[i] * II.velocity[i] + II.sigma_11[i] * I_1.velocity[i]) - mesh.A[i - 1] * (I_1.sigma_11[i - 1] * II.velocity[i - 1] + II.sigma_11[i - 1] * I_1.velocity[i - 1]))
				+ para.dt / delta * (mesh.A[i] * II.sigma_11[i] * II.velocity[i] - mesh.A[i - 1] * II.sigma_11[i - 1] * II.velocity[i - 1]) * C_1_next.big_lambda_1[i]
				+ para.dt * (C_1.big_phi_1[i] * C.velocity[i] + C.sigma_22[i] * C_1_next.big_omega_1[i]) / C.density[i]
				+ para.dt * (para.g * C_1.velocity[i]) * flag_g
				);
		}

		if (M.Y0[i] != 0) {
			REAL W1, W2, W3, W4, W4_;
			/*****s11*****/
			W1 = C_1.s_11[i];	W2 = 4.0 / 3 * M.mu[i];	W3 = (I_1.velocity[i] - I_1.velocity[i - 1]);
			if (para.alpha == 1)
				W4 = -2.0 / 3 * M.mu[i] * C_1.big_omega_1[i];
			else if (para.alpha == 2)
				W4 = -2.0 / 3 * M.mu[i] * C_1.big_omega_1[i] * mesh.R[i] / mesh.r[i] - 2.0 / 3 * M.mu[i] * C_1.velocity[i] / mesh.r[i];
			else
				W4 = -2.0 / 3 * M.mu[i] * C_1.big_omega_1[i] * mesh.R[i] / mesh.r[i] - 4.0 / 3 * M.mu[i] * C_1.velocity[i] / mesh.r[i];
			if (para.alpha == 1)
				W4_ = -2.0 / 3 * M.mu[i] * C_1_next.big_omega_1[i];
			else if (para.alpha == 2)
				W4_ = -2.0 / 3 * M.mu[i] * C_1_next.big_omega_1[i] * mesh_next.R[i] / mesh_next.r[i] - 2.0 / 3 * M.mu[i] * C_1_next.velocity[i] / mesh_next.r[i];
			else
				W4_ = -2.0 / 3 * M.mu[i] * C_1_next.big_omega_1[i] * mesh_next.R[i] / mesh_next.r[i] - 4.0 / 3 * M.mu[i] * C_1_next.velocity[i] / mesh_next.r[i];

			C_1_next.s_11[i] = W1 + para.dt / mesh.dr[i] * W2 * W3 + para.dt * W4;

			/*****s22*****/
			REAL tempX = 4.0 / 3; //4/3  1/3
			W1 = C_1.s_22[i];	W2 = -2.0 / 3 * M.mu[i];	W3 = (I_1.velocity[i] - I_1.velocity[i - 1]);
			if (para.alpha == 1)
				W4 = tempX * M.mu[i] * C_1.big_omega_1[i];
			else if (para.alpha == 2)
				W4 = tempX * M.mu[i] * C_1.big_omega_1[i] * mesh.R[i] / mesh.r[i] + 4.0 / 3 * M.mu[i] * C_1.velocity[i] / mesh.r[i];
			else
				W4 = tempX * M.mu[i] * C_1.big_omega_1[i] * mesh.R[i] / mesh.r[i] + 2.0 / 3 * M.mu[i] * C_1.velocity[i] / mesh.r[i];
			if (para.alpha == 1)
				W4_ = tempX * M.mu[i] * C_1_next.big_omega_1[i];
			else if (para.alpha == 2)
				W4_ = tempX * M.mu[i] * C_1_next.big_omega_1[i] * mesh_next.R[i] / mesh_next.r[i] + 4.0 / 3 * M.mu[i] * C_1_next.velocity[i] / mesh_next.r[i];
			else
				W4_ = tempX * M.mu[i] * C_1_next.big_omega_1[i] * mesh_next.R[i] / mesh_next.r[i] + 2.0 / 3 * M.mu[i] * C_1_next.velocity[i] / mesh_next.r[i];
			C_1_next.s_22[i] = W1 + para.dt / mesh.dr[i] * W2 * W3 + para.dt * W4;

			/*****s12*****/
			W1 = pow(mesh.R[i], para.beta / 2) * C_1.big_phi_1[i];	W2 = (C.s_11[i] - C.s_22[i]) / 2.0 + M.mu[i];
			W3 = (I_1.big_omega_1[i] * pow(mesh.R_half[i], para.beta / 2) - I_1.big_omega_1[i - 1] * pow(mesh.R_half[i - 1], para.beta / 2));
			if (para.alpha == 1)
				W4 = ((C.s_11[i] - C.s_22[i]) / 2.0 - M.mu[i]) * para.omega * C_1.velocity[i];
			else
				W4 = ((C.s_11[i] - C.s_22[i]) / 2.0 - M.mu[i]) * (para.omega * C_1.velocity[i] / mesh.r[i] + C_1.big_omega_1[i] * mesh.R[i] / mesh.r[i]);
			if (para.alpha == 1)
				W4_ = ((C_next.s_11[i] - C_next.s_22[i]) / 2.0 - M.mu[i]) * para.omega * C_1_next.velocity[i];
			else
				W4_ = ((C_next.s_11[i] - C_next.s_22[i]) / 2.0 - M.mu[i]) * (para.omega * C_1_next.velocity[i] / mesh_next.r[i] + C_1_next.big_omega_1[i] * mesh_next.R[i] / mesh_next.r[i]);

			C_1_next.big_phi_1[i] = (W1 + para.dt / mesh.dr[i] * W2 * W3 + para.dt * W4) / pow(mesh_next.R[i], para.beta / 2);
		}



		if (para.EOS_set == 2) {
			//C_1_next.pressure[i] = ((M.big_gamma0[i] - 1.0) * (C_1_next.t_energy[i] - C_next.velocity[i] * C_1_next.velocity[i]) - C_next.pressure[i] * C_1_next.volume[i]) * C_next.density[i];
			C_1_next.pressure[i] = -(M.big_gamma0[i] - 1.0) * C_next.density[i] * C_next.density[i] * C_next.i_energy[i] * C_1_next.volume[i]
				+ (M.big_gamma0[i] - 1.0) * C_next.density[i] * (C_1_next.t_energy[i] - C_next.velocity[i] * C_1_next.velocity[i]);
		}
		else if (para.EOS_set == 1) {
			C_1_next.pressure[i] = M.density_ini[i] * M.c0[i] * M.c0[i] * C_1_next.volume[i] * (
				(M.density_ini[i] * (-1.0 - M.s[i] + M.big_gamma0[i] + (M.s[i] - M.big_gamma0[i]) * M.density_ini[i] * C_next.volume[i]))
				/ pow(1.0 + M.s[i] * (-1.0 + M.density_ini[i] * C_next.volume[i]), 3)
				)
				+ M.density_ini[i] * M.big_gamma0[i] * (C_1_next.t_energy[i] - C_next.velocity[i] * C_1_next.velocity[i])
				;
		}
		if (M.Y0[i] == 0) {
			C_1_next.sigma_11[i] = -C_1_next.pressure[i];
			C_1_next.sigma_22[i] = -C_1_next.pressure[i];
		}
		else {
			C_1_next.sigma_11[i] = C_1_next.s_11[i] - C_1_next.pressure[i];
			C_1_next.sigma_22[i] = C_1_next.s_22[i] - C_1_next.pressure[i];
		}
	}
}


