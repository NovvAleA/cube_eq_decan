#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <fstream>
#include "float.h"
#include "fluid.h"
#include "functions_for_modeling.h"
#include "phase_diagram.h"


int main() {
    std::string gas = "decan";
    std::vector<double> P_V1_V2_V3;
    fluid_data decan {0.4890,617.7,2119999.96};
    std::string gas1 = "gecsan";
    fluid_data gecsan {0.3003,507.6,3.03e+6};
    double T = 520;
    double P = 1000000;
    std::vector<double> rt = solver_P_of_V(decan, T, P);
    //saturation_pressure(decan,T);
    double T_array[] {-100.0,-50.0,0.0, 50.0, 100.0,150.0,200.0,250.0}; //C

    std::vector<fluid_data> fluids(2);
    fluids[0] = gecsan;
    fluids[1] = decan;
    compute_phase_diagram(fluids,520);
    std::ofstream fout;
    /*
    fout.open("/home/ome123/Документы/Наука/Курсовая/saturation_pressure.txt");
    fout << "Decan C10, Saturation Pressure (Pa), V1, V1, V3" <<std::endl;


    for (int i = 0; i < std::size(T_array); i ++) {
        P_V1_V2_V3 = saturation_pressure(decan,T_array[i] + 273.15);;
        fout << T_array[i] << "  " <<P_V1_V2_V3[0]<< "  " << P_V1_V2_V3[1]<< "  " << P_V1_V2_V3[2]<<"  " <<  P_V1_V2_V3[3] << std::endl;

    }
    fout<<std::endl;*/

    /*gas = "C16";
    fluid_data C16 {0.742,722,1409977.86};
    fout << "C16: T, Saturation Pressure (Pa), V1, V1, V3" <<std::endl;


    for (int i = 0; i < std::size(T_array); i ++) {
        P_V1_V2_V3 = saturation_pressure(C16,T_array[i] + 273.15);;
        fout << T_array[i] << "  " <<P_V1_V2_V3[0]<< "  " << P_V1_V2_V3[1]<< "  " << P_V1_V2_V3[2]<<"  " <<  P_V1_V2_V3[3] << std::endl;

    }
    fout<<std::endl;
    fout.close();*/

    fout.open("/home/ome123/Документы/Наука/Курсовая/cube_res_points.txt");
    for (int i = 0; i < 100; i ++) {
        rt = solver_P_of_V(decan, T, 3e+5 + double(30000*i));
        for (int j = 0; j < rt.size(); j++) {
            fout << rt[j]<< " " << 3e+5 + double(30000*i)<< std::endl;
        }
    }
    fout.close();

    /*
    std::cout << "Temperature_(K) Pressure_(bar) Molar_volume_(m3/mol)" <<std::endl;
    for (int i = 0; i < rt.v.size(); i++) {
       std::cout << T << "  " << P/100000 << "  " << rt.v[i]<< std::endl;
    }
    std::ofstream resout;          // поток для записи
    resout.open("/home/ome123/Документы/Наука/Курсовая/PR_solve.txt");
    resout << "Temperature_(K) Pressure_(bar) Molar_volume_(m3/mol)" <<std::endl;
    for (int i = 0; i < rt.v.size(); i++) {
        resout << T << "  " << P/100000 << "  " << rt.v[i]<< std::endl;
    }
    resout.close();
    std::ofstream fout;          // поток для записи
    fout.open("/home/ome123/Документы/Наука/Курсовая/cube_res_points.txt");

    for (int i = 0; i < 40; i ++) {
        rt = solver_P_of_V(decan, T, 1700000.0 + double(15000*i));
        for (int j = 0; j < rt.v.size(); j++) {
            fout << rt.v[j]<< " " << 1700000.0 + double(15000*i)<< std::endl;
        }
    }
    fout.close();

    std::ofstream resout1;          // поток для записи
    resout1.open("/home/ome123/Документы/Наука/Курсовая/PR_solve_1.txt");
    double P_arr[] {1000000,1100000,1200000,1300000,1318300,1400000,1500000, 1600000 };
    resout1 << "Temperature_(K) Pressure_(bar) Molar_volume_(m3/mol)" <<std::endl;
    for (int k = 0; k < 8; k ++) {
        rt = solver_P_of_V(decan, T, P_arr[k]);
        for (int i = 0; i < rt.v.size(); i++) {
            resout1 << T << "  " << P_arr[k]/100000 << "  " << rt.v[i]<< std::endl;
        }
    }
    resout1.close();
   */

    return 0;
}

// TIP See CLion help at <a
// href="https://www.jetbrains.com/help/clion/">jetbrains.com/help/clion/</a>.
//  Also, you can try interactive lessons for CLion by selecting
//  'Help | Learn IDE Features' from the main menu.