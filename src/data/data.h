#pragma once
#include <string>
#include <vector>
#include <complex>
using std::string;
using std::vector;
typedef double REAL;
typedef double DATATYPE;

struct Para {
	int cells;				
	int now_step = 0;
	vector<REAL> domain_len;
	REAL cfl;
	REAL timeout;
	REAL now_time = 0;
	REAL lbc;
	REAL rbc;
	REAL alpha = 1;	
	REAL beta;
	int EOS_set = 1;

	/*IO*/
	string IOmode = "bin";	
	int store_steps;
	string out_file_dir;
	string loading_file_dir;
	REAL change_t = 1;
	REAL change_p = 1;
	int printing = 1;
	/*计算方案*/
	int mode_space = 1;	
	int mode_riemann = 1;
	REAL q;		
	int fix = 0;			
	int mode_time = 1;		

	int pt = 0;				
	int cl = 0;				
	int tem = 0;			

	int fo = 0;			
	REAL omega;			
	REAL k;			
	REAL g = 0;		
	REAL dt = 1;			
	REAL max_ce;			
	int RTI = 0;

	REAL offset = 0.0;		

	int n_threads;			
};

/*
	I0  I1  I2  I3  I4	     In-1 In	
	|---|---|---|---|.....|---|---|
	 C1  C2  C3  C4        Cn-1 Cn		
*/

struct Mesh {
	vector<DATATYPE> R;		
	vector<DATATYPE> R_half;	
	vector<DATATYPE> A;		
	vector<DATATYPE> dR;		

	vector<DATATYPE> r;		
	vector<DATATYPE> r_half;	
	vector<DATATYPE> dr;	

	vector<DATATYPE> m;
	vector<DATATYPE> m_half;
	vector<DATATYPE> dm;		

	vector<DATATYPE> M;
	vector<DATATYPE> M_half;
	vector<DATATYPE> dM;	

	void _malloc(int n) {
		R.assign(n, 0);		R_half.assign(n - 1, 0);	dR.assign(n, 0);	A.assign(n, 0);
		r.assign(n, 0);		r_half.assign(n - 1, 0);	dr.assign(n, 0);
		m.assign(n, 0);		m_half.assign(n - 1, 0);	dm.assign(n, 0);
		M.assign(n, 0);		M_half.assign(n - 1, 0);	dM.assign(n, 0);
	}

};

struct Cell_Value {
	/*density = dmess / delta_x*/
	vector<DATATYPE> density;		
	vector<DATATYPE> volume;	
	/*velocity*/
	vector<DATATYPE> velocity;		

	vector<DATATYPE> sigma_11;		
	vector<DATATYPE> s_11;			
	vector<DATATYPE> sigma_22;		
	vector<DATATYPE> s_22;			
	vector<DATATYPE> sigma_12;		
	vector<DATATYPE> s_12;			
	vector<DATATYPE> pressure;		

	vector<DATATYPE> i_energy;		
	vector<DATATYPE> t_energy;		

	vector<DATATYPE> Tem;			
	vector<DATATYPE> Tem2;		
	vector<DATATYPE> state;		

	vector<DATATYPE> is_Y;	


	vector<DATATYPE> big_lambda_1;	//Λ  
	//vector<DATATYPE> lambda_1;		//Λ=y_1*ik  
	vector<DATATYPE> big_omega_1;	//Ω 
	vector<DATATYPE> big_phi_1;		//Φ=σ12_1*ik

	vector<DATATYPE> ce;
	void _malloc(int n) {
		density.assign(n, 0);	volume.assign(n, 0);
		velocity.assign(n, 0);
		sigma_11.assign(n, 0);	s_11.assign(n, 0);
		sigma_22.assign(n, 0);	s_22.assign(n, 0);
		sigma_12.assign(n, 0);	s_12.assign(n, 0);
		pressure.assign(n, 0);
		i_energy.assign(n, 0);
		t_energy.assign(n, 0);
		Tem.assign(n, 298);	Tem2.assign(n, 298); state.assign(n, -1);
		big_lambda_1.assign(n, 0);	big_omega_1.assign(n, 0);
		//lambda_1.assign(n, 0);	omega_1.assign(n, 0);
		big_phi_1.assign(n, 0);
		is_Y.assign(n, 0);
		ce.assign(n, 0);
	}
	void copy_i_to_j(int i, int j) {
		density[j] = density[i];	volume[j] = volume[i];
		velocity[j] = velocity[i];
		sigma_11[j] = sigma_11[i];	s_11[j] = s_11[i];
		sigma_22[j] = sigma_22[i];	s_22[j] = s_22[i];
		sigma_12[j] = sigma_12[i];	s_12[j] = s_12[i];
		pressure[j] = pressure[i];
		i_energy[j] = i_energy[i];
		t_energy[j] = t_energy[i];
		Tem[j] = Tem[i]; Tem2[j] = Tem2[i]; state[j] = state[i];
		big_lambda_1[j] = big_lambda_1[i];    big_omega_1[j] = big_omega_1[i];
		//lambda_1[j] = lambda_1[i];    omega_1[j] = omega_1[i];
		big_phi_1[j] = big_phi_1[i];
		is_Y[j] = is_Y[i];
		ce[j] = ce[i];
	}
};

/*材料常数*/
struct Material_Constant {
	vector<DATATYPE> density_ini;		
	vector<DATATYPE> volume_ini;
	vector<DATATYPE> c0;			
	vector<DATATYPE> big_gamma0;		
	vector<DATATYPE> s;					
	vector<DATATYPE> mu;				
	vector<DATATYPE> Y0;				

	vector<DATATYPE> Cv;				
	vector<DATATYPE> beta;	

	void _malloc(int n) {
		density_ini.assign(n, 0); volume_ini.assign(n, 0);
		c0.assign(n, 0);
		big_gamma0.assign(n, 0);
		s.assign(n, 0);
		mu.assign(n, 0);
		Y0.assign(n, 0);
		Cv.assign(n, 0);
		beta.assign(n, 0);
	}
	void copy_i_to_j(int i, int j) {
		density_ini[j] = density_ini[i]; volume_ini[j] = volume_ini[i];
		c0[j] = c0[i];
		big_gamma0[j] = big_gamma0[i];
		s[j] = s[i];
		mu[j] = mu[i];
		Y0[j] = Y0[i];
		Cv[j] = Cv[i];
		beta[j] = beta[i];
	}
	void set_constant(int i, DATATYPE _rho0, DATATYPE _g0, DATATYPE _c0, DATATYPE _s, DATATYPE _mu, DATATYPE _Y0) {
		density_ini[i] = _rho0; volume_ini[i] = 1.0 / _rho0;
		big_gamma0[i] = _g0;
		c0[i] = _c0;
		s[i] = _s;
		mu[i] = _mu;
		Y0[i] = _Y0;
	}
	void set_constant(int i, string m) {
		if (m == "Al") set_constant(i, 2698, 2.18, 5333, 1.356, 27.6e9, 3.0e8);//Al
		else if (m == "Ta") set_constant(i, 16754, 1.82, 3430, 1.19, 65.25e9, 0.77e9);//Ta
		else if (m == "Cu") set_constant(i, 8932, 1.97, 3899, 1.534, 44e9, 33.3e6);//Cu
		else if (m == "LiF") set_constant(i, 2638, 1.63, 5130, 1.31, 55.14e9, 1.8e8);//LiF
	}
};
/*
*				  rho0  g0    c0     s      mu      Y0
M.set_constant(i, 2698, 2.18, 5333, 1.356, 27.6e9, 3.0e8);//Al
M.set_constant(i, 16754, 1.82, 3430, 1.19, 65.25e9, 0.77e9);//Ta
M.set_constant(i, 8932, 1.97, 3899, 1.534, 44e9, 33.3e6);//Cu
M.set_constant(i, 2638, 1.63, 5130, 1.31, 55.14e9, 1.8e8);//LiF
*/
//F(U)
struct Interface {
	vector<DATATYPE> velocity;
	vector<DATATYPE> density;
	vector<DATATYPE> density_L;
	vector<DATATYPE> density_R;
	vector<DATATYPE> pressure;
	vector<DATATYPE> s_11;
	vector<DATATYPE> sigma_11;
	vector<DATATYPE> sigma_22;
	vector<DATATYPE> big_omega_1;
	vector<DATATYPE> big_phi_1;		//Φ=σ12_1*ik
	vector<DATATYPE> sigma_velocity;
	vector<DATATYPE> i_energy;
	vector<DATATYPE> t_energy;
	void _malloc(int n) {
		velocity.assign(n, 0);
		density.assign(n, 0); density_L.assign(n, 0); density_R.assign(n, 0);
		pressure.assign(n, 0);
		s_11.assign(n, 0);
		sigma_11.assign(n, 0);
		sigma_velocity.assign(n, 0);
		i_energy.assign(n, 0);
		t_energy.assign(n, 0);
		sigma_22.assign(n, 0); big_omega_1.assign(n, 0); big_phi_1.assign(n, 0);
	}
	void copy_i_to_j(int i, int j) {
		velocity[j] = velocity[i];
		density[j] = density[i];
		density_L[j] = density_L[i];
		density_R[j] = density_R[i];
		pressure[j] = pressure[i];
		s_11[j] = s_11[i];
		sigma_11[j] = sigma_11[i];
		sigma_velocity[j] = sigma_velocity[i];
		i_energy[j] = i_energy[i];
		t_energy[j] = t_energy[i];
		sigma_22[j] = sigma_22[i]; big_omega_1[j] = big_omega_1[i]; big_phi_1[j] = big_phi_1[i];
	}
};
