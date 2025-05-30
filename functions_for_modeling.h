
#ifndef FUNCTIONS_FOR_MODELING_H
#define FUNCTIONS_FOR_MODELING_H
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <fstream>
#include "float.h"
#include "fluid.h"

double alpha(double T, double Tc,double w) {
    //w > 0.49!!!
    double x = 1 + (0.37464 + 1.54226 * w - 0.2699 * w * w)*(1 - sqrt((T/Tc)));
    return pow(x,2);
}

double a_PR( fluid_data fluid, double T ) {
    return 0.45724*(fluid.R*fluid.R*fluid.T_c*fluid.T_c/(fluid.P_c))*alpha(T,fluid.T_c,fluid.w);
}



std::vector<double> solve_cubic_equation(double a, double b, double c) {
    //y^3 + 3ay^2 + 3b^y+c = 0
    std::vector<double> roots(3);
    double p = pow(a,2)-b;//3*b в дипломе было
    double q = pow(a,3) + (c - 3*a*b)/2;



    double Q = q*q - p*p*p;

    if (Q > 0) {
        double sqQ = sqrt(Q);
        roots[0] = cbrt(-q + sqQ) - cbrt(q + sqQ);
        roots.pop_back();
        roots.pop_back();

    }
    if (Q < 0) {
        double sq_p = sqrt(p);
        double phi = acos(sqrt((q*q)/(p*p*p)));

        roots[0] = 2*sq_p*cos(phi/3);
        roots[1] = 2*sq_p*cos((phi + 2*M_PI)/3);
        roots[2] = 2*sq_p*cos((phi + 2*M_PI*2)/3);

    }
    if (Q == 0) {
        if (p == 0 && q==0) {
            roots[0] = cbrt(-q) - cbrt(q);
            roots[1] = roots[0];
            roots[2] = roots[0];
        }
        else{
            double phi = acos(-q/pow(p,1.5));
            double sq_p = sqrt(p);
            roots[0] = 2*sq_p*cos(phi/3);
            roots[1] = 2*sq_p*cos((phi + 2*M_PI)/3);
            roots[2] = 2*sq_p*cos((phi + 2*M_PI*2)/3);
        }
    }
    //если мы уже вызваем решатель второй раз, когда у нас корни очень отличаются, то нам не надо идти дальше и решать уравнение с заменой 1/y

    for (size_t i = 0; i < roots.size(); i++) {
        roots[i] -= a;
    }

    std::sort(roots.rbegin(), roots.rend());

    return roots;
}
double P_PR (fluid_data& fluid, double V, double T){
    double a = 0.45724*(fluid.R*fluid.R*fluid.T_c*fluid.T_c/(fluid.P_c))*alpha(T,fluid.T_c,fluid.w);
    double P = fluid.R*T/(V - fluid.b) - a/(V*V + 2*V*fluid.b - fluid.b*fluid.b);
    return P;
}
std::vector<double> solver_P_of_V(fluid_data& fluid, double T, double P) {
    std::vector<double> rt, rt1;
    int quantity = 3;
    double y1 = 0;
    double y2 = 0;
    double y3 = 0;

    double v1 = 0;
    double v2 = 0;
    double v3 = 0;

    double A = 0.45724*alpha(T,fluid.T_c,fluid.w)*P*fluid.T_c*fluid.T_c/(T*T*fluid.P_c);
    double B = P*fluid.b/(T*fluid.R);
    double b_gas_param = 0.07780*fluid.R*fluid.T_c/(fluid.P_c);

    //coefficients of the cubic polynomial
    double a =(B - 1)/3;
    double b =(A - 2*B - 3*B*B)/3;
    double c = -A*B + B*B + B*B*B;

    rt = solve_cubic_equation(a,b,c);
    /*rt[0] = (fluid.R*T)*(rt[0] -a)/P;
    if (rt.size() != 1) {
        rt[1] = (fluid.R*T)*(rt[1] -a)/P;
        rt[2] = (fluid.R*T)*(rt[2] -a)/P;
    }*/
    rt[0] = (fluid.R*T)*rt[0]/P;
    if (rt.size() != 1) {
        rt[1] = (fluid.R*T)*rt[1]/P;
        rt[2] = (fluid.R*T)*rt[2]/P;
    }

    if (rt.size() > 1) {
        std::sort(rt.begin(), rt.end());
        if (P < 10) {
            if (c != 0) {
                std::vector<double> roots1(3);
                roots1 = solve_cubic_equation(b/c, a/c, 1/c);
                /*for (int i = 0 ; i  < roots1.size(); i++ ) {
                    roots1[i] = 1/(roots1[i] - b/c);
                }*/
                for (int i = 0 ; i  < roots1.size(); i++ ) {
                    roots1[i] = 1/(roots1[i]);
                }
                roots1[0] = (fluid.R*T)*roots1[0]/P;
                if (rt.size() != 1) {
                    roots1[1] = (fluid.R*T)*roots1[1]/P;
                    roots1[2] = (fluid.R*T)*roots1[2]/P;
                }
                std::sort(roots1.begin(), roots1.end());
                rt[0] = roots1[0];
                rt[1] = roots1[1];
            }
        }
    }

    //вот зщесь корень меньше б получается
    std::sort(rt.rbegin(), rt.rend());

    for (int i = rt.size() - 1; i  >= 0; i-- ) {
        if (abs(P - P_PR(fluid,rt[i],T))/P > 0.05 ) {
            std::cout << "residual > 0.001" << std::endl;
        }
    }

    for (int i = rt.size() - 1; i  >= 0; i-- ) {
        if (rt[i] <= b_gas_param) {
            rt.pop_back();
        }
    }



    return rt;
}



std::pair<double,double> max_min_pressure(fluid_data& fluid, double T) {
    std::pair<double,double> P = {0,0};

    double V_left_for_Pmax=0, V_right_for_Pmax=0,Vmiddle=0,V_left_for_Pmin=0,V_right_for_Pmin = 0;
    double P_left, P_right, Pmin, Pmax;
    double left,right;
    double P_i_m1=0, P_i=0, P_i_p1=0;

    double max_Vright = 0.042;//value bigger then Vmax
    double h = 0.0000001;
    int N = int((max_Vright - fluid.b)/h);


    P_i_m1 = P_PR(fluid,fluid.b + h, T);
    P_i = P_PR(fluid,fluid.b + 2*h, T);

    for (int i = 2; i < N-1; i++) {
        P_i_p1 = P_PR(fluid,fluid.b + h*(i+1), T);
        if (P_i_m1 <= P_i && P_i >= P_i_p1) {
            V_left_for_Pmax = fluid.b + h*(i - 1);
            V_right_for_Pmax = fluid.b + h*(i + 1);
        }
        if (P_i_m1 >= P_i && P_i <= P_i_p1) {
            V_left_for_Pmin = fluid.b + h*(i - 1);
            V_right_for_Pmin = fluid.b + h*(i + 1);
        }
        P_i_m1 = P_i;
        P_i = P_i_p1;
    }



    int i = 0;
    P_left = P_PR(fluid,V_left_for_Pmin ,T );
    P_right = P_PR(fluid,V_right_for_Pmin,T);


    while ((P_left - P_right) != 0) {
        if (i > 10000000) {
            std::cout<<"i = 10000 in  max_min_pressure"<< std::endl;
            break;
        }
        Vmiddle = (V_left_for_Pmin + V_right_for_Pmin)/2;
        //Pmin = ;
        P_left = P_PR(fluid,Vmiddle - DBL_EPSILON,T);
        P_right = P_PR(fluid,Vmiddle + DBL_EPSILON,T);
        if (P_left >P_right) {
            V_left_for_Pmin = Vmiddle;
        }
        else if (P_left < P_right) {
            V_right_for_Pmin = Vmiddle;
        }
        i++;
    }

    P.first = P_PR(fluid,Vmiddle,T);


    //find P max
    //Pmax = 3.49153*101325.0;
    i = 0;
    P_left = P_PR(fluid,V_left_for_Pmax,T );
    P_right = P_PR(fluid,V_right_for_Pmax,T);

    while ((P_left - P_right) != 0) {
        if (i > 10000000) {
            std::cout<<"i = 10000 in  max_min_pressure"<< std::endl;
            break;
        }
        Vmiddle = (V_left_for_Pmax + V_right_for_Pmax)/2;
        //Pmin = ;
        P_left = P_PR(fluid,Vmiddle - DBL_EPSILON,T);
        P_right = P_PR(fluid,Vmiddle + DBL_EPSILON,T);
        if (P_left <P_right) {
            V_left_for_Pmax = Vmiddle;
        }
        else if (P_left > P_right) {
            V_right_for_Pmax = Vmiddle;
        }
        i++;
    }


    P.second = P_PR(fluid,Vmiddle,T);
    return std::make_pair(P.first,P.second);
}

double integral_from_pressure(fluid_data& fluid, double T, double Va,double Vb, double a) {

    return fluid.R*T*log((Vb - fluid.b)/(Va - fluid.b)) - (a/(2*sqrt(2)*fluid.b))*log((((Vb - (-1 + sqrt(2))*fluid.b)*(Va - (-1 - sqrt(2))*fluid.b))/((Va - (-1 + sqrt(2))*fluid.b)*(Vb - (-1 - sqrt(2))*fluid.b))));
}

std::vector<double> saturation_pressure(fluid_data& fluid, double T) {
    std::vector<double> P_V1_V2_V3(4);
    double Ps, Pmin,Pmax;
    std::pair<double,double> Pminmax;
    std::vector<double>  rt;
    double a = 0.45724*(fluid.R*fluid.R*fluid.T_c*fluid.T_c/(fluid.P_c))*alpha(T,fluid.T_c,fluid.w);
    double S_BDE = 0;
    double S_EFG = 1;


    Pminmax = max_min_pressure(fluid,T);
    Pmin = Pminmax.first;
    Pmax = Pminmax.second;
    if (Pmin<0) {
        Pmin = 0;
    }

    int i = 0;
    while (S_BDE != S_EFG) {
        if (i > 10000) {
            std::cout<<"i = 10000 in  saturation_pressure"<< std::endl;
            break;
        }
        i++;

        Ps = (Pmin + Pmax)/2;
        rt = solver_P_of_V(fluid,T,Ps);
        std::sort(rt.begin(), rt.end());


        if (rt.size() < 3) {
            std::cout<<"rt.v.size() < 3 in saturation_pressure"<< std::endl;
        }

        //S_BDE = (Ps - Pmin_)*(rt.v[1] - rt.v[0]) - integral_from_pressure(fluid,T,rt.v[0],rt.v[1],a);
        S_BDE = Ps*(rt[1] - rt[0]) - integral_from_pressure(fluid,T,rt[0],rt[1],a);
        S_EFG = integral_from_pressure(fluid,T,rt[1],rt[2],a) - Ps*(rt[2] - rt[1]);

        if (S_BDE > S_EFG) {
            Pmax = Ps;
        }
        else {
            Pmin = Ps;
        }
    }

    P_V1_V2_V3[0] = Ps;
    P_V1_V2_V3[1] = rt[0];
    P_V1_V2_V3[2] = rt[1];
    P_V1_V2_V3[3] = rt[2];

    return P_V1_V2_V3;
}



#endif //FUNCTIONS_FOR_MODELING_H
