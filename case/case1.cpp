#include "../src/data/use_data.h"
DATATYPE cal_eos(int flag, int index, Cell_Value& C, Material_Constant& M);
#define Pi    3.14159265358979323846
void _malloc(int cells)
{
	mesh._malloc(cells + 2);		if (para.fo) mesh_1._malloc(cells + 2);
	C._malloc(cells + 2);			if (para.fo) C_1._malloc(cells + 2);
	M._malloc(cells + 2);			if (para.fo) M_1._malloc(cells + 2);
	II._malloc(cells + 1);			if (para.fo) I_1._malloc(cells + 1);
}
void uniform_grid()
{
	DATATYPE dx = para.domain_len.back() / (DATATYPE)para.cells;
	for (int i = 1; i <= para.cells; ++i) {
		mesh.R_half[i] = (DATATYPE)i * dx + para.offset;
		mesh.R[i] = (DATATYPE)i * dx - 0.5 * dx + para.offset;
		mesh.dR[i] = dx;

		mesh.r_half[i] = (DATATYPE)i * dx + para.offset;
		mesh.r[i] = (DATATYPE)i * dx - 0.5 * dx + para.offset;
		mesh.dr[i] = dx;

		M.density_ini[i] = C.density[i];
	}
}
void uniform_grid_m()
{
	for (int i = 1; i <= para.cells; ++i) {
		M.density_ini[i] = C.density[i];

		mesh.dm[i] = C.density[i] * mesh.dR[i];

		//DATATYPE A = 1;
		//if (para.alpha == 1.0) A = 1;
		//else if (para.alpha == 2.0) A = 2 * Pi * mesh.R[i];
		//else if (para.alpha == 3.0) A = 4 * Pi * mesh.R[i] * mesh.R[i];
		//mesh.dM[i] = A * C.density[i] * mesh.dR[i];

		if (para.alpha == 3.0)
			mesh.dM[i] = C.density[i] * mesh.dR[i] * (pow(mesh.R_half[i], 2) + mesh.R_half[i] * mesh.R_half[i - 1] + pow(mesh.R_half[i - 1], 2)) / 3.0;
		else if (para.alpha == 2.0) {
			mesh.dM[i] = C.density[i] * mesh.dR[i] * (mesh.R_half[i] + mesh.R_half[i - 1]) / 2.0;
		}
		else if (para.alpha == 1.0) {
			mesh.dM[i] = C.density[i] * mesh.dR[i];
		}
	}
}


//Jaouen S. A purely Lagrangian method for computing linearly-perturbed flows in spherical geometry[J]. Journal of Computational Physics, 2007, 225(1): 464-490.
void case_RMI_1D()
{
	para.domain_len = { 3,
						4,
						10 };

	para.cells = 1000;
	para.fo = 1;
	para.k = 10;

	_malloc(para.cells);
	para.cfl = 0.45;
	para.timeout = 4;
	para.lbc = 1;
	para.rbc = 0;
	para.out_file_dir = "./Data/Out/";
	para.IOmode = "bin";	para.store_steps = 1000;	para.printing = 100;

	para.alpha = 3;	
	para.EOS_set = 2;
	para.mode_space = 2;
	para.mode_time = 1;	
	para.mode_riemann = 1;
	para.n_threads = 8;

	if (para.alpha == 3) para.omega = para.k * (para.k + 1.0);
	else para.omega = para.k * para.k;


	uniform_grid();
	for (int i = 1; i <= para.cells; ++i) {
		if (mesh.R[i] <= para.domain_len[0]) {
			C.density[i] = 4;
			C.velocity[i] = 0;
			C.pressure[i] = 1;
			M.big_gamma0[i] = 3;
		}
		else if (mesh.R[i] <= para.domain_len[1]) {
			C.density[i] = 1;
			C.velocity[i] = 0;
			C.pressure[i] = 1;
			M.big_gamma0[i] = 1.5;
		}
		else {
			DATATYPE s = 0.5;
			C.pressure[i] = 1.0 / (1.0 - s);
			C.density[i] = (2.5 * C.pressure[i] + 0.5) / (0.5 * C.pressure[i] + 2.5);
			C.velocity[i] = -sqrt((C.pressure[i] - 1.0) * (1.0 - 1.0 / C.density[i]));
			M.big_gamma0[i] = 1.5;
			//C.density[i] = 5.5/3.5;
			//C.velocity[i] = -sqrt(2.0 / 2.5);
			//C.pressure[i] = 2;
			//M.big_gamma0[i] = 1.5;
		}
		C.volume[i] = 1.0 / C.density[i];
		mesh.dm[i] = C.density[i] * mesh.dr[i];
		C.s_11[i] = 0;	C.sigma_11[i] = -C.pressure[i]; C.sigma_22[i] = -C.pressure[i];
		C.i_energy[i] = cal_eos(1, i, C, M);		C.t_energy[i] = C.i_energy[i] + 0.5 * pow(C.velocity[i], 2.0);
	}
	uniform_grid_m();

	mesh_1.R_half[300] = 1;
	for (int i = 1; i <= para.cells; ++i) {
		//mesh_1.R_half[i] = exp(-10.0 * pow(mesh.R[i] - 3.0, 2));
		//mesh_1.R_half[300] = 1;
		if (para.alpha == 3.0)
			C_1.big_lambda_1[i] = -M.density_ini[i] * (pow(mesh.R_half[i], 2) * mesh_1.R_half[i] - pow(mesh.R_half[i - 1], 2) * mesh_1.R_half[i - 1]) / mesh.dM[i];
		else if (para.alpha == 2.0)
			C_1.big_lambda_1[i] = -M.density_ini[i] * (pow(mesh.R_half[i], 1) * mesh_1.R_half[i] - pow(mesh.R_half[i - 1], 1) * mesh_1.R_half[i - 1]) / mesh.dM[i];
		else
			C_1.big_lambda_1[i] = -M.density_ini[i] * (mesh_1.R_half[i] - mesh_1.R_half[i - 1]) / mesh.dM[i];
	}
}


void case_RTI_EE_planar_geometry()
{
	REAL	rho1 = 24.3;// 24.3 10.8   6.3
	REAL	rho2 = 2.7;//
	REAL	At = (rho1 - rho2) / (rho1 + rho2);
	REAL	gg = 40;//
	REAL	G1 = rho1;//
	REAL	G2 = rho1 * 0.1;
	REAL	T = G2 / G1;
	REAL	n = 5;//
	REAL	R = 10;//
	REAL	Br = R * rho1 * gg / G1;

	REAL lam = 2 * 3.14 / n;

	REAL tempx = 2;
	para.domain_len = { R,
						tempx * R };		

	para.cells = 1000;		
	para.fo = 1;	
	para.k = n;	//mm-1   1~10


	_malloc(para.cells);
	para.cfl = 0.4;	
	para.timeout = 2;	
	para.lbc = 0;
	para.rbc = 0;

	para.out_file_dir = "./Data/out/";
	para.IOmode = "bin";	para.store_steps = 1000;	para.printing = 1000;

	para.alpha = 1;				
	para.EOS_set = 1;			
	para.mode_space = 2;		
	para.mode_time = 1;			
	para.mode_riemann = 1;
	para.n_threads = 4;
	para.RTI = 1;

	if (para.alpha == 3) para.k--;

	if (para.alpha == 3) para.omega = para.k * (para.k + 1.0);
	else para.omega = para.k * para.k;

	uniform_grid();
	para.g = -gg;
	for (int i = 1; i <= para.cells; ++i) {
		if (mesh.R[i] <= para.domain_len[0]) {
			M.set_constant(i, rho2, 1.97, 5.333, 1.356, G2, 999999999e9);
			C.velocity[i] = 0;
		}
		else {
			M.set_constant(i, rho1, 1.97, 5.333, 1.356, G1, 999999999e9);//Al
			C.velocity[i] = 0;
		}
	}
	for (int i = para.cells; i >= 1; --i) {
		C.pressure[i - 1] = C.pressure[i] - M.density_ini[i] * mesh.dR[i] * para.g;
		//M.Y0[i] = 0;
		C.density[i] = M.density_ini[i];	C.volume[i] = 1.0 / C.density[i];
		mesh.dm[i] = C.density[i] * mesh.dr[i];
		C.s_11[i] = 0;	C.sigma_11[i] = -C.pressure[i]; C.s_22[i] = 0; C.sigma_22[i] = -C.pressure[i];
		C.i_energy[i] = cal_eos(1, i, C, M);		C.t_energy[i] = C.i_energy[i] + 0.5 * pow(C.velocity[i], 2.0);
	}

	uniform_grid_m();

	mesh_1.R_half[para.cells / 2] = 1;
	for (int i = 1; i <= para.cells; ++i) {
		//mesh_1.R_half[i] = exp(-tempx * 10 * pow(mesh.R[i] - para.domain_len[0], 2));
		if (para.alpha == 3.0)
			C_1.big_lambda_1[i] = -M.density_ini[i] * (pow(mesh.R_half[i], 2) * mesh_1.R_half[i] - pow(mesh.R_half[i - 1], 2) * mesh_1.R_half[i - 1]) / mesh.dM[i];
		else if (para.alpha == 2.0)
			C_1.big_lambda_1[i] = -M.density_ini[i] * (pow(mesh.R_half[i], 1) * mesh_1.R_half[i] - pow(mesh.R_half[i - 1], 1) * mesh_1.R_half[i - 1]) / mesh.dM[i];
		else if (para.alpha == 1.0)
			C_1.big_lambda_1[i] = -M.density_ini[i] * (mesh_1.R_half[i] - mesh_1.R_half[i - 1]) / mesh.dM[i];
	}
}

void case_RTI_EE_cylindrical_and_spherical_geometry()
{
	REAL	rho1 = 24.3;// 24.3 10.8   6.3
	REAL	rho2 = 2.7;//
	REAL	At = (rho1 - rho2) / (rho1 + rho2);
	REAL	gg = 1;//
	REAL	G1 = rho1;//
	REAL	G2 = rho1 * 0.1;//
	REAL	T = G2 / G1;
	REAL	n = 5;//
	REAL	R = 25.0;//
	REAL	Br = R * rho1 * gg / G1;

	REAL tempx = 1.5;
	REAL lam = 2 * 3.1415926 / n; if (para.alpha != 1) lam *= R;


	para.domain_len = { R,
						tempx * R };

	para.cells = 1000;		
	para.fo = 1;
	//para.omega = 10 * 11;
	para.k = n;	//mm-1   1~10


	_malloc(para.cells);	
	para.cfl = 0.45;	
	para.timeout = 40;		/*tm=10/γ */
	para.lbc = 0;			
	para.rbc = 2;
	para.out_file_dir = "./Data/out/";
	para.IOmode = "bin";	para.store_steps = 1000;	para.printing = 5000;

	para.alpha = 2;				
	para.EOS_set = 1;			
	para.mode_space = 2;		
	para.mode_time = 1;			
	para.mode_riemann = 1;		
	para.n_threads = 4;
	para.RTI = 1;


	if (para.alpha == 3) para.omega = para.k * (para.k + 1.0);
	else para.omega = para.k * para.k;

	uniform_grid();

	para.g = -gg;
	for (int i = 1; i <= para.cells; ++i) {
		if (mesh.R[i] <= para.domain_len[0]) {
			M.set_constant(i, rho2, 1.97, 5.333, 1.356, G2, 999999999);
			C.velocity[i] = 0;
		}
		else {
			//std::cout << i << std::endl;
			//system("pause");
			M.set_constant(i, rho1, 1.97, 5.333, 1.356, G1, 999999999);//Al
			C.velocity[i] = 0;
		}
	}
	for (int i = para.cells; i >= 1; --i) {
		C.pressure[i - 1] = C.pressure[i] - M.density_ini[i] * mesh.dR[i] * para.g;
		//M.Y0[i] = 0;
		C.density[i] = M.density_ini[i];	C.volume[i] = 1.0 / C.density[i];
		mesh.dm[i] = C.density[i] * mesh.dr[i];
		C.s_11[i] = 0;	C.sigma_11[i] = -C.pressure[i]; C.s_22[i] = 0; C.sigma_22[i] = -C.pressure[i];
		C.i_energy[i] = cal_eos(1, i, C, M);		C.t_energy[i] = C.i_energy[i] + 0.5 * pow(C.velocity[i], 2.0);
	}

	uniform_grid_m();

	mesh_1.R_half[(int)(para.cells / tempx)] = 1;  // 1000 667
	//mesh_1.R_half[200] = 1;
	for (int i = 1; i <= para.cells; ++i) {
		//mesh_1.R_half[i] = exp(-1 * pow(mesh.R[i] - para.domain_len[0], 2));
		if (para.alpha == 3.0)
			C_1.big_lambda_1[i] = -M.density_ini[i] * (pow(mesh.R_half[i], 2) * mesh_1.R_half[i] - pow(mesh.R_half[i - 1], 2) * mesh_1.R_half[i - 1]) / mesh.dM[i];
		else if (para.alpha == 2.0)
			C_1.big_lambda_1[i] = -M.density_ini[i] * (pow(mesh.R_half[i], 1) * mesh_1.R_half[i] - pow(mesh.R_half[i - 1], 1) * mesh_1.R_half[i - 1]) / mesh.dM[i];
		else if (para.alpha == 1.0)
			C_1.big_lambda_1[i] = -M.density_ini[i] * (mesh_1.R_half[i] - mesh_1.R_half[i - 1]) / mesh.dM[i];
	}
}



