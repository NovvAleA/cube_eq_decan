//
// Created by ome123 on 22.05.25.
//

#ifndef PHASE_DIAGRAM_H
#define PHASE_DIAGRAM_H
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <fstream>
#include "float.h"
#include "fluid.h"
#include "functions_for_modeling.h"


double Ki_from_Ps(double P, double Ps) {
    return Ps / P;
}

// Шаг 0: Начальные приближения K_i по формуле Вильсона
void wilsonKi(const std::vector<fluid_data> &fluids, std::vector<double> &K,  double T, double P) {
    for (int i = 0; i < fluids.size(); i++) {
        K[i] =  (fluids[i].P_c / P) * exp(5.31 * (1 + fluids[i].w) * (1 - fluids[i].T_c / T)) ;
    }
}

//Шаг 1
//Уравнение Рачфорда–Райса для случая двухкомпонентной смеси
//Здесь V - молярная доля
double analytic_V(double K1,double K2,double z1,double z2) {
    double V = 0;
    double denominator = (K1 - 1) * (K2 - 1) * (z1 + z2);
    if (abs(denominator) < 1e-12) {
        std::cout<<"Деление на почти ноль в уранвении Рачфорда–Райса";
    }
    V = ((1 - K1) * z1 + (1 - K2) * z2)/ denominator;
    return V;
}

double computeSi(const std::vector<std::vector<double>>& Aij, const std::vector<double>& z, size_t i) {
    double Si = 0.0;
    size_t N = z.size();
    for (size_t j = 0; j < N; ++j) {
        Si += Aij[i][j] * z[j];
    }
    return Si;
}

void computePhaseConcentrations(
    const std::vector<double>& z,   // исходные молярные доли компонентов
    const std::vector<double>& K,   // коэффициенты распределения K_i
    double V,                       // мольная доля газовой фазы//
    std::vector<double>& x,
    std::vector<double>& y
) {
    std::vector<double> D = std::vector<double>(z.size(), 0);
    for (size_t i = 0; i < x.size(); ++i) {
        D[i] = 1.0 + (K[i] - 1.0) * V;
        x[i] = z[i] / D[i];
        y[i] = K[i] * x[i];
    }
}
double a_mixture(const std::vector<fluid_data> &fluids, const std::vector<double> &z,const double T){
    double a = 0;
    double k_ij = 0;
    std::vector<double> a_i(fluids.size());

    for (size_t i = 0; i < fluids.size(); ++i) {
        a_i[i] = a_PR(fluids[i], T);
    }
    std::vector<std::vector<double>> a_ij( fluids.size(), std::vector<double>(fluids.size(), 0));
    for (size_t i = 0; i < fluids.size(); ++i) {
        for (size_t j = 0; j < fluids.size(); ++j) {
            a += (1 - k_ij) * sqrt(a_i[i]*a_i[j])*z[i]*z[j];
        }
    }
    return a;
}

double b_mixture(const std::vector<fluid_data> &fluids, const std::vector<double> &z){
    double b = 0;
    for (size_t i = 0; i < fluids.size(); ++i) {
        b += fluids[i].b*z[i];
    }
    return b;
}

double S_i(const std::vector<fluid_data> &fluids, const std::vector<double> &z,const size_t i, const double T){
    double s_i = 0;
    double k_ij = 0;
    std::vector<double> a_i(fluids.size());

    for (size_t i = 0; i < fluids.size(); ++i) {
        a_i[i] = a_PR(fluids[i], T);
    }
    std::vector<std::vector<double>> a_ij( fluids.size(), std::vector<double>(fluids.size(), 0));
    for (size_t j = 0; j < fluids.size(); ++j) {
        s_i += (1 - k_ij) * sqrt(a_i[i]*a_i[j])*z[j];
    }

    return s_i;
}
// Вычисление логарифма коэффициента летучести ln(phi_i)
double chemical_potential(
    double Z,
    double a,
    double b,
    double bi,
    double Si
) {
    if (Z - b <= 0.0 || (Z + (1 - sqrt(2)) * b) <= 0.0 || (Z + (1 + sqrt(2)) * b) <= 0.0) {
        throw std::domain_error("Invalid arguments for logarithm in ln(phi_i) formula");
    }

    double ln_phi = 0.0;

    ln_phi += -std::log(Z - b);
    ln_phi += (a / (2*sqrt(2) * b))*((2.0 * Si / a - bi / b) *
              std::log((Z + (1 - sqrt(2)) * b) / (Z + (1 + sqrt(2)) * b)));
    ln_phi += bi * (Z - 1.0) / b;

    return ln_phi;
}



//z - молярные доли
void compute_phase_diagram(const std::vector<fluid_data> &fluids,const double T) {
    std::vector<double> z(fluids.size(), 0);

    double quant_iter = 0;
    std::vector<double> phi_L = std::vector<double>(fluids.size(), 0);
    std::vector<double> phi_V = std::vector<double>(fluids.size(), 0);
    std::vector<double> Z(fluids.size());
    std::vector<double> x = std::vector<double>(fluids.size(), 0);
    std::vector<double> y = std::vector<double>(fluids.size(), 0);
    std::vector<double> K = std::vector<double>(fluids.size(), 0);

    double V = 0;
    double A = 0;
    double B = 0;
    double a_ec = 0; //ecuasion coef
    double b_ec = 0;
    double c_ec = 0;
    double Si = 0;
    double diff = 0;

    double a = -1;
    double b = -1;

    //double a_mixture = a_mixture(fluids, z, T);
    //double b_mixture = b_mixture(fluids, z);

    std::ofstream phaseConcentrations;
    phaseConcentrations.open("phase_diagram.txt");
    //step 0

    int N_nodes = 1000;
    double P = 1000000;
    wilsonKi(fluids, K, T, P);
    for (int сoncentration_iter = 0; сoncentration_iter < N_nodes; сoncentration_iter++) {
        z[0] = double (сoncentration_iter)/double(N_nodes);
        z[1] = 1 - z[0];
        a = a_mixture(fluids, z, T);
        b = b_mixture(fluids, z);
        while (quant_iter < 10000) {
            //step 1
            V = analytic_V(K[0],K[1],z[0],z[1]);

            //step 2
            computePhaseConcentrations(z, K, V, x,y);

            //step 3-4
            for (size_t i = 0; i < fluids.size(); ++i) {
                A = (P*a)/(pow(R_gas_const*T,2));
                B = P*b/(T*R_gas_const);
                a_ec =(B - 1)/3;
                b_ec =(A - 2*B - 3*B*B)/3;
                c_ec = -A*B + B*B + B*B*B;
                Z = solve_cubic_equation(a_ec,b_ec,c_ec);

                Si = S_i(fluids, x, i, T);
                phi_L[i] = chemical_potential(Z[Z.size()-1], a,b, fluids[i].b,Si);//самый левый корень
                Si = S_i(fluids, y, i, T);
                phi_V[i] = chemical_potential(Z[0],  a,b, fluids[i].b,Si);//
            }
            for (size_t i = 0; i < fluids.size(); ++i) {
                diff += (x[i]*phi_L[i])/(y[i]*phi_V[i]);
            }
            if (abs(diff - 1.0) < 1e6) {
                break;
            }
            diff = 0;

            //step 6
            for (size_t i = 0; i < fluids.size(); ++i) {
                K[i]= phi_L[i]/phi_V[i];
            }
            quant_iter++;
        }
        phaseConcentrations << z[0] << " " << Z[0] << " " << Z[Z.size() -1] << std::endl;
    }



}

#endif //PHASE_DIAGRAM_H
