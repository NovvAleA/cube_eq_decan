//
// Created by ome123 on 22.05.25.
//

#ifndef FLUID_H
#define FLUID_H
const double R_gas_const = 8.31446261815324;
struct fluid_data {
public:
    //const std::string name;
    double R = 8.31446261815324;
    double w;
    double T_c;
    double P_c;
    double b;
    fluid_data() {
        w = 1;
        T_c = 0;
        P_c = 0;
        b = 1;
    }
    fluid_data( double w1, double T_c1, double P_c1) {
        w = w1;
        T_c = T_c1;
        P_c = P_c1;
        b = 0.07780 * R * T_c1 / P_c1;
    }
};



#endif //FLUID_H
