#include "../data/use_data.h"
#include <algorithm>
#include <iostream>
using std::max;
using std::min;
#define Pi    3.14159265358979323846
DATATYPE sound_velocity(DATATYPE d, DATATYPE p, DATATYPE dst, DATATYPE e, int i);

void update_dt()
{
	para.dt = 1e10; para.max_ce = 0;

	for (int i = 1; i <= para.cells; ++i) {
		DATATYPE ce = sound_velocity(C.density[i], C.pressure[i], C.s_11[i], C.i_energy[i], i);	para.max_ce = std::max(para.max_ce, ce);
		C.ce[i] = ce;
		para.dt = min(para.dt, para.cfl * (mesh.dR[i]) / (abs(C.velocity[i]) + abs(ce)));
		//para.dt = min(para.dt, para.cfl * (mesh.dR[i]) / (abs(ce)));

		if (para.fo) {
			para.dt = min(para.dt, para.cfl * (2 * 3.14159 / para.k) / (400 * abs(ce)));
		}
	}

	REAL max_uu = -1.0;
	for (int i = 1; i <= para.cells; ++i) {
		if (II.velocity[i] < II.velocity[i - 1]) {
			max_uu = std::max(max_uu, abs(II.velocity[i] - II.velocity[i - 1]) / mesh.dR[i]);
			//if (para.cfl / max_uu < para.dt) {
			//	para.dt = para.cfl / max_uu;
			//	//std::cout << "max_uu" << std::endl;
			//	//system("pause");
			//}
			para.dt = min(para.dt, para.cfl / max_uu);
		}
	}

	//para.dt = 1e-11;
	para.now_time += para.dt;	para.now_step++;
}