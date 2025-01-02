/*********************************************************************************************
LPC-Elastic   liuyongliu
*********************************************************************************************/
#include "data/data.h"
#include <time.h>
#include <iostream>
#include <omp.h>
using std::cout;

Para para;
Mesh mesh;			Mesh mesh_1;
Cell_Value C;		Cell_Value C_1;
Material_Constant M;	Material_Constant M_1;
Interface II;		Interface I_1;


void open_files(string IOmode);
void close_files();
void Write_current_layer_data_to_file(int ss, string mode);
void update_1euler_3();
void case_RMI_1D();
void case_RTI_EE_planar_geometry();
void case_RTI_EE_cylindrical_and_spherical_geometry();

int main()
{
	clock_t start, now, end;
	start = clock();

	case_RMI_1D();
	//case_RTI_EE_planar_geometry();
	//case_RTI_EE_cylindrical_and_spherical_geometry();


	open_files(para.IOmode);
	int ii = 0;
	while (para.now_time < para.timeout && para.now_time >= 0.0 && para.dt > 0.0)
	{
		if ((para.now_time / para.timeout) >= ((REAL)ii / para.store_steps)) {
			ii++;
			Write_current_layer_data_to_file(1, para.IOmode);
		}



		if (para.mode_time == 1)
			update_1euler_3();
		//else if (para.mode_time == 3)
			//update_3RK();


		if (para.now_step % para.printing == 0)
		{
			now = clock();
			printf("now_step=%d\t now_time=%e\t %.1f%%\t dt=%e\t max_c=%e\t%.fs\n", para.now_step, para.now_time, 100.0 * para.now_time / para.timeout, para.dt, para.max_ce, REAL(now - start) / CLOCKS_PER_SEC);
		}
	}
	Write_current_layer_data_to_file(1, para.IOmode);
	printf("now_step=%d\t now_time=%e\t %.1f%%\t dt=%e\t max_ce=%e\n", para.now_step, para.now_time, 100.0 * para.now_time / para.timeout, para.dt, para.max_ce);

	end = clock();
	cout << "computational time = " << REAL(end - start) / CLOCKS_PER_SEC << "s" << '\t' << REAL(end - start) / CLOCKS_PER_SEC / 60 << "min" << '\t' << REAL(end - start) / CLOCKS_PER_SEC / 60 / 60 << "h" << std::endl;
	cout << "cells=" << para.cells << " CFL=" << para.cfl << " timeout=" << para.timeout << std::endl;

	close_files();

	system("pause");
	return 0;
}

