#include "../data/use_data.h"
#include <iostream>
#include <fstream>
std::ofstream outfile1, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7, outfile8, outfile9, outfile10, outfile11;
std::ofstream outfile12, outfile13, outfile14, outfile15, outfile16, outfile17;
std::ofstream outfile18, outfile19, outfile20, outfile21, outfile22, outfile23;
std::ofstream outfile24, outfile25, outfile26, outfile27, outfile28, outfile29;
using namespace std;

void open_files(string IOmode)
{
	if (IOmode == "bin") {
		outfile1.open(para.out_file_dir + "density_all.dat", ios::out | ios::binary);
		outfile2.open(para.out_file_dir + "velocity_all.dat", ios::out | ios::binary);
		outfile3.open(para.out_file_dir + "pressure_all.dat", ios::out | ios::binary);
		outfile4.open(para.out_file_dir + "s_x_all.dat", ios::out | ios::binary);
		outfile5.open(para.out_file_dir + "sigma_x_all.dat", ios::out | ios::binary);
		outfile6.open(para.out_file_dir + "xpoint_all.dat", ios::out | ios::binary);
		outfile7.open(para.out_file_dir + "t_all.dat", ios::out | ios::binary);
		outfile8.open(para.out_file_dir + "i_energy_all.dat", ios::out | ios::binary);
		outfile9.open(para.out_file_dir + "t_energy_all.dat", ios::out | ios::binary);
		outfile12.open(para.out_file_dir + "Tem_all.dat", ios::out | ios::binary);
		outfile13.open(para.out_file_dir + "mesh.dat", ios::out | ios::binary);
		outfile14.open(para.out_file_dir + "Tem2_all.dat", ios::out | ios::binary);
		outfile15.open(para.out_file_dir + "Cv_all.dat", ios::out | ios::binary);

		outfile16.open(para.out_file_dir + "r_1.dat", ios::out | ios::binary);
		outfile17.open(para.out_file_dir + "p_1.dat", ios::out | ios::binary);
		outfile18.open(para.out_file_dir + "velocity_1.dat", ios::out | ios::binary);
		outfile19.open(para.out_file_dir + "volume_1.dat", ios::out | ios::binary);
		outfile20.open(para.out_file_dir + "big_omega_1.dat", ios::out | ios::binary);
		outfile21.open(para.out_file_dir + "big_lambda_1.dat", ios::out | ios::binary);

		outfile22.open(para.out_file_dir + "s_y_all.dat", ios::out | ios::binary);
		outfile23.open(para.out_file_dir + "is_Y.dat", ios::out | ios::binary);

		outfile24.open(para.out_file_dir + "big_phi_1.dat", ios::out | ios::binary);
		outfile25.open(para.out_file_dir + "ce.dat", ios::out | ios::binary);

		outfile26.open(para.out_file_dir + "sigma11_1.dat", ios::out | ios::binary);
		outfile27.open(para.out_file_dir + "sigma22_1.dat", ios::out | ios::binary);


		if (outfile1.is_open() == 0) {
			cout << "error\n";
			system("pause");
			return;
		}
	}

}
void close_files()
{
	outfile1.close(); outfile2.close(); outfile3.close(); outfile4.close(); outfile5.close(); outfile6.close(); outfile7.close(); outfile8.close(); outfile9.close(); outfile10.close(), outfile11.close();
	outfile12.close(); outfile13.close(); outfile14.close(); outfile15.close(); outfile16.close(); outfile17.close();
	outfile18.close(); outfile19.close(); outfile20.close(); outfile21.close(); outfile22.close(); outfile23.close();
	outfile24.close(); outfile25.close(); outfile26.close(); outfile27.close(); outfile28.close();
}

void Write_current_layer_data_to_file(int ss, string mode)
{
	if (para.now_step % ss != 0) return;

	if (mode == "bin") {
		outfile1.write(((char*)&C.density[1]), para.cells * sizeof(DATATYPE));
		outfile2.write(((char*)&C.velocity[1]), para.cells * sizeof(DATATYPE));
		outfile3.write(((char*)&C.pressure[1]), para.cells * sizeof(DATATYPE));
		outfile4.write(((char*)&C.s_11[1]), para.cells * sizeof(DATATYPE));
		outfile5.write(((char*)&C.sigma_11[1]), para.cells * sizeof(DATATYPE));
		outfile6.write(((char*)&mesh.R[1]), para.cells * sizeof(DATATYPE));
		vector<DATATYPE> nowtime;	nowtime.assign(para.cells, para.now_time);
		outfile7.write(((char*)&nowtime[0]), para.cells * sizeof(DATATYPE));
		outfile8.write(((char*)&C.i_energy[1]), para.cells * sizeof(DATATYPE));
		outfile9.write(((char*)&C.t_energy[1]), para.cells * sizeof(DATATYPE));
		outfile12.write(((char*)&C.Tem[1]), para.cells * sizeof(DATATYPE));
		outfile13.write(((char*)&mesh.R_half[0]), (para.cells + 1) * sizeof(DATATYPE));
		outfile14.write(((char*)&C.Tem2[1]), para.cells * sizeof(DATATYPE));
		outfile15.write(((char*)&M.Cv[1]), para.cells * sizeof(DATATYPE));

		if (para.fo) outfile16.write(((char*)&mesh_1.R_half[1]), para.cells * sizeof(DATATYPE));
		if (para.fo) outfile17.write(((char*)&C_1.pressure[1]), para.cells * sizeof(DATATYPE));
		if (para.fo) outfile18.write(((char*)&C_1.velocity[1]), para.cells * sizeof(DATATYPE));
		if (para.fo) outfile19.write(((char*)&C_1.volume[1]), para.cells * sizeof(DATATYPE));
		if (para.fo) outfile20.write(((char*)&C_1.big_omega_1[1]), para.cells * sizeof(DATATYPE));
		if (para.fo) outfile21.write(((char*)&C_1.big_lambda_1[1]), para.cells * sizeof(DATATYPE));

		outfile22.write(((char*)&C.s_22[1]), para.cells * sizeof(DATATYPE));
		outfile23.write(((char*)&C.is_Y[1]), para.cells * sizeof(DATATYPE));

		if (para.fo) outfile24.write(((char*)&C_1.big_phi_1[1]), para.cells * sizeof(DATATYPE));
		outfile25.write(((char*)&C.ce[1]), para.cells * sizeof(DATATYPE));

		if (para.fo) outfile26.write(((char*)&C_1.sigma_11[1]), para.cells * sizeof(DATATYPE));
		if (para.fo) outfile27.write(((char*)&C_1.sigma_22[1]), para.cells * sizeof(DATATYPE));


	}
	else if (mode == "ascii") {

	}

}

