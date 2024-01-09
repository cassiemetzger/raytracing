#include <iostream>
#include <cmath>
#include "rk4_adaptive.h"
#include "EoM.h"
#include <fstream>
#include "dormand_prince_scary.h"
int main() {
    //first test case
    
    double a=0;
    double M=1.0;
    double eta = 0.48857874;
    double R = 10.0 * M;
    boyer_lindquist_metric metric1(a, M);
    auto dydx = [&metric1, &eta] (double x, const std::vector<double> &y) {
        metric1.compute_metric(y[0], y[1]); // y[0] is r, y[1] is theta
        // compute the right hand sides of equations (17)-(22)
        std::vector<double> equations(6);
        double ut = std::sqrt((metric1.gamma11*y[3]*y[3]) + (metric1.gamma22 * y[4] * y[4]) + (metric1.gamma33*y[5]*y[5]))/metric1.alpha;
        equations[0] = metric1.gamma11 * y[3]/ut;
        equations[1] = metric1.gamma22 * y[4]/ut;
        equations[2] = metric1.gamma33 * y[5]/ut - metric1.beta3;
        equations[3] = -metric1.alpha * ut * metric1.d_alpha_dr + y[5]*metric1.d_beta3_dr - (1.0/(2.0*ut))*(y[3]*y[3]*metric1.d_gamma11_dr + y[4]*y[4]*metric1.d_gamma22_dr + y[5]*y[5]*metric1.d_gamma33_dr);
        equations[4] = -metric1.alpha * ut * metric1.d_alpha_dth + y[5]*metric1.d_beta3_dth - (1.0/(2.0*ut))*(y[3]*y[3]*metric1.d_gamma11_dth + y[4]*y[4]*metric1.d_gamma22_dth + y[5]*y[5]*metric1.d_gamma33_dth);
        equations[5] = 0;
        return equations;
    };
    auto stop_test = [] (double x, const std::vector<double>& y) {return (y[2]>2*M_PI || y[0] < 2.01);};
    double h = 0.001;
    double x0 = 0.0;
    double th0=M_PI/2.0;
    double phi0=0.0;
    metric1.compute_metric(R, th0);
    double u_r = -std::cos(eta)*std::sqrt(metric1.g_11);
    double u_th = 0.0;
    double u_phi = std::sin(eta) * (std::sqrt(metric1.g_33));
    std::vector<double> y0 = {R,th0 ,phi0,u_r,u_th,u_phi};

    rk4_adaptive rk4(6, 1e-11, 1e-11);
    rk4.integrate(dydx, stop_test, h, x0, y0);
    std::ofstream output_file("output_test_case_1.csv");
    for (int i = 0; i < rk4.xs.size(); i++) {
        output_file << rk4.xs[i] << "," << rk4.result[i][0] <<","<< rk4.result[i][1]<<","<<rk4.result[i][2] <<std::endl;
    };
    output_file.close();

    std::ofstream output_file2("output_test_case_1xy.csv");
    for (int i = 0; i < rk4.xs.size(); i++) {
        output_file2 << rk4.result[i][0] * std::cos(rk4.result[i][2]) << ", " << rk4.result[i][0] * std::sin(rk4.result[i][2]) << std::endl;
    };
    output_file2.close();
    //second orbit
     eta = 0.24904964;
    R = 20.0 * M;
    auto dydx7 = [&metric1, &eta] (double x, const std::vector<double> &y) {
        metric1.compute_metric(y[0], y[1]); // y[0] is r, y[1] is theta
        // compute the right hand sides of equations (17)-(22)
        std::vector<double> equations(6);
        double ut = std::sqrt((metric1.gamma11*y[3]*y[3]) + (metric1.gamma22 * y[4] * y[4]) + (metric1.gamma33*y[5]*y[5]))/metric1.alpha;
        equations[0] = metric1.gamma11 * y[3]/ut;
        equations[1] = metric1.gamma22 * y[4]/ut;
        equations[2] = metric1.gamma33 * y[5]/ut - metric1.beta3;
        equations[3] = -metric1.alpha * ut * metric1.d_alpha_dr + y[5]*metric1.d_beta3_dr - (1.0/(2.0*ut))*(y[3]*y[3]*metric1.d_gamma11_dr + y[4]*y[4]*metric1.d_gamma22_dr + y[5]*y[5]*metric1.d_gamma33_dr);
        equations[4] = -metric1.alpha * ut * metric1.d_alpha_dth + y[5]*metric1.d_beta3_dth - (1.0/(2.0*ut))*(y[3]*y[3]*metric1.d_gamma11_dth + y[4]*y[4]*metric1.d_gamma22_dth + y[5]*y[5]*metric1.d_gamma33_dth);
        equations[5] = 0;
        return equations;
    };
    auto stop_test78 = [] (double x, const std::vector<double>& y) {return (y[2]>4*M_PI || y[0] < 2.01);};
    
    metric1.compute_metric(R, th0);
    u_r = -std::cos(eta)*std::sqrt(metric1.g_11);
    u_th = 0.0;
    u_phi = std::sin(eta) * (std::sqrt(metric1.g_33));
    y0 = {R,th0 ,phi0,u_r,u_th,u_phi};
    rk4.integrate(dydx7, stop_test78, h, x0, y0);
    std::ofstream output_file98("output_test_case_2.csv");
    for (int i = 0; i < rk4.xs.size(); i++) {
        output_file98 << rk4.xs[i] << "," << rk4.result[i][0] <<","<< rk4.result[i][1]<<","<<rk4.result[i][2] <<std::endl;
    };
    output_file98.close();

    std::ofstream output_file99("output_test_case_2xy.csv");
    for (int i = 0; i < rk4.xs.size(); i++) {
        output_file99 << rk4.result[i][0] * std::cos(rk4.result[i][2]) << ", " << rk4.result[i][0] * std::sin(rk4.result[i][2]) << std::endl;
    };
    output_file99.close();
    // second test case
    a=1.0;
     boyer_lindquist_metric metric2(a, M);
     auto dydx2 = [&metric2, &eta] (double x, const std::vector<double> &y) {
        metric2.compute_metric(y[0], y[1]); // y[0] is r, y[1] is theta
        // compute the right hand sides of equations (17)-(22)
        std::vector<double> equations(6);
        double ut = std::sqrt((metric2.gamma11*y[3]*y[3]) + (metric2.gamma22 * y[4] * y[4]) + (metric2.gamma33*y[5]*y[5]))/metric2.alpha;
        equations[0] = metric2.gamma11 * y[3]/ut;
        equations[1] = metric2.gamma22 * y[4]/ut;
        equations[2] = metric2.gamma33 * y[5]/ut - metric2.beta3;
        equations[3] = -metric2.alpha * ut * metric2.d_alpha_dr + y[5]*metric2.d_beta3_dr - (1.0/(2.0*ut))*(y[3]*y[3]*metric2.d_gamma11_dr + y[4]*y[4]*metric2.d_gamma22_dr + y[5]*y[5]*metric2.d_gamma33_dr);
        equations[4] = -metric2.alpha * ut * metric2.d_alpha_dth + y[5]*metric2.d_beta3_dth - (1.0/(2.0*ut))*(y[3]*y[3]*metric2.d_gamma11_dth + y[4]*y[4]*metric2.d_gamma22_dth + y[5]*y[5]*metric2.d_gamma33_dth);
        equations[5] = 0;
        return equations;
    };

    //orbit A
    R = 1.0 + sqrt(2.0);
    u_th = std::sqrt(11.0+8*std::sqrt(2.0));
    u_r =0.0;
    u_phi=0.0;
    th0 = M_PI/2.0;
    phi0=0;
    std::vector<double> y1 = {R,th0 ,phi0,u_r,u_th,u_phi};
    auto stop_test2 = [&R] (double x, const std::vector<double>& y) {return (y[0]>R+R/5.0 || y[0]<R-R/5.0);};
    x0= 0; 
    rk4.integrate(dydx2, stop_test2, h, x0, y1);
    std::ofstream output_file4("output_test_case_2A_error.csv");
    for (int i = 0; i < rk4.xs.size(); i++) {
        output_file4 << rk4.xs[i] << "," << std::abs(rk4.result[i][0]-R)/R <<std::endl;
    };
    output_file4.close();
    std::ofstream output_file3("output_test_case_2A.csv");
    for (int i = 0; i < rk4.xs.size(); i++) {
        double r = rk4.result[i][0];
        double th = rk4.result[i][1];
        double phi = rk4.result[i][2];
        output_file3 << std::sqrt(r*r+a*a)*std::sin(th)*std::cos(phi) << ", " <<std::sqrt(r*r+a*a)*std::sin(th)*sin(phi) << ","<< r*std::cos(th)<<std::endl;
    };
    output_file3.close(); 
    //orbit B
    R = 1.0 + std::sqrt(3.0);
    u_th = std::sqrt(12+8*std::sqrt(3));
    u_phi = -1.0;
    std::vector<double>y2 = {R,th0,phi0,u_r,u_th,u_phi};
    rk4.integrate(dydx2,stop_test2,h,x0,y2);
    std::ofstream output_file7("output_test_case_2B_error.csv");
    for (int i = 0; i < rk4.xs.size(); i++) {
        output_file7 << rk4.xs[i] << "," << std::abs(rk4.result[i][0]-R)/R <<std::endl;
    };
    output_file7.close();
    std::ofstream output_file8("output_test_case_2B.csv");
    for (int i = 0; i < rk4.xs.size(); i++) {
        double r = rk4.result[i][0];
        double th = rk4.result[i][1];
        double phi = rk4.result[i][2];
        output_file8 << std::sqrt(r*r+a*a)*std::sin(th)*std::cos(phi) << ", " <<std::sqrt(r*r+a*a)*std::sin(th)*sin(phi) << ","<< r*std::cos(th)<<std::endl;
    };
    output_file8.close();
    // orbit C 
    R = 3.0; 
    u_th = std::sqrt(27.0); 
    u_phi = -2.0; 
    std::vector<double>y3 = {R,th0,phi0,u_r,u_th,u_phi};
    rk4.integrate(dydx2, stop_test2, h, x0, y3);
    std::ofstream output_file5("output_test_case_2C_error.csv");
    for(int i =0; i<rk4.xs.size(); i++){
        output_file5 << rk4.xs[i] << "," << std::abs(rk4.result[i][0]-R)/R <<std::endl;
    };
    output_file5.close();
    std::ofstream output_file6("output_test_case_2C.csv");
    for(int i = 0; i < rk4.xs.size(); i++){
        double r3 = rk4.result[i][0];
        double th3 = rk4.result[i][1];
        double phi3= rk4.result[i][2];
        output_file6 << std::sqrt(r3*r3+a*a)*std::sin(th3)*std::cos(phi3) << ", " <<std::sqrt(r3*r3+a*a)*std::sin(th3)*sin(phi3) << ","<< r3*std::cos(th3)<<std::endl;
    }; 
    output_file6.close();

    //orbit D

    R = 1.0 + 2.0*std::sqrt(2.0);
    u_th = std::sqrt(-13.0+16.0*std::sqrt(2.0));
    u_phi= -6.0;
    std::vector<double>y4 = {R,th0,phi0,u_r,u_th,u_phi};
    rk4.integrate(dydx2, stop_test2, h, x0, y4);
    std::ofstream output_file9("output_test_case_2D_error.csv");
    for(int i =0; i<rk4.xs.size(); i++){
        output_file9 << rk4.xs[i] << "," << std::abs(rk4.result[i][0]-R)/R <<std::endl;
    };
    output_file9.close();
    std::ofstream output_file10("output_test_case_2D.csv");
    for(int i = 0; i < rk4.xs.size(); i++){
        double r3 = rk4.result[i][0];
        double th3 = rk4.result[i][1];
        double phi3= rk4.result[i][2];
        output_file10 << std::sqrt(r3*r3+a*a)*std::sin(th3)*std::cos(phi3) << ", " <<std::sqrt(r3*r3+a*a)*std::sin(th3)*sin(phi3) << ","<< r3*std::cos(th3)<<std::endl;
    }; 
    output_file10.close();
    //orbit E
    R = 2.0;
    u_th =std::sqrt(16.0);
    u_phi = 1.0;
    std::vector<double>y5 = {R,th0,phi0,u_r,u_th,u_phi};
    rk4.integrate(dydx2, stop_test2, h, x0, y5);
    std::ofstream output_file11("output_test_case_2E_error.csv");
    for(int i =0; i<rk4.xs.size(); i++){
        output_file11 << rk4.xs[i] << "," << std::abs(rk4.result[i][0]-R)/R <<std::endl;
    };
    output_file11.close();
    std::ofstream output_file12("output_test_case_2E.csv");
    for(int i = 0; i < rk4.xs.size(); i++){
        double r3 = rk4.result[i][0];
        double th3 = rk4.result[i][1];
        double phi3= rk4.result[i][2];
        output_file12 << std::sqrt(r3*r3+a*a)*std::sin(th3)*std::cos(phi3) << ", " <<std::sqrt(r3*r3+a*a)*std::sin(th3)*sin(phi3) << ","<< r3*std::cos(th3)<<std::endl;
    }; 
    output_file12.close(); 
    return 0;
} 