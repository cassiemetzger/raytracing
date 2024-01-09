#include<iostream> 
#include<cmath>
#include "rk4_adaptive.h"
#include "EoM.h"
#include <fstream>

int main() {
     double a=0.0; //what alex used 
    double th0= 40 * M_PI/180.0; //observer is located here 
    double Rin = 5.0; //inner radius of accretion disk
    double Rout = 20.0; //outer radius of accretion disk 
    double D = 500.0; //distance to observer (from alex)
    double M = 1.0; //mass of black hole 
    double h = 0.0001; //step size
     double x0 = 0.0;
     double tolerance = 1e-2;//tolerance for how close to pi/2 its getting
     double biggestR = std::sqrt(25*25+25*25+D*D);
     double r_H = M + std::sqrt(M*M +a*a);
     double n = 1.0/8.0;//length of step
     double im_length = 20.0/n;
     double im_width = 50.0/n;
    rk4_adaptive rk4(6, 1e-11, 1e-11);
    //recall metrics to reset variables 
    boyer_lindquist_metric metric3(a, M);
    auto dydx3 = [&metric3] (double x, const std::vector<double> &y) {
        metric3.compute_metric(y[0], y[1]); // y[0] is r, y[1] is theta
        // compute the right hand sides of equations (17)-(22)
        std::vector<double> equations(6);
        double ut = std::sqrt((metric3.gamma11*y[3]*y[3]) + (metric3.gamma22 * y[4] * y[4]) + (metric3.gamma33*y[5]*y[5]))/metric3.alpha;
        equations[0] = metric3.gamma11 * y[3]/ut;
        equations[1] = metric3.gamma22 * y[4]/ut;
        equations[2] = metric3.gamma33 * y[5]/ut - metric3.beta3;
        equations[3] = -metric3.alpha * ut * metric3.d_alpha_dr + y[5]*metric3.d_beta3_dr - (1.0/(2.0*ut))*(y[3]*y[3]*metric3.d_gamma11_dr + y[4]*y[4]*metric3.d_gamma22_dr + y[5]*y[5]*metric3.d_gamma33_dr);
        equations[4] = -metric3.alpha * ut * metric3.d_alpha_dth + y[5]*metric3.d_beta3_dth - (1.0/(2.0*ut))*(y[3]*y[3]*metric3.d_gamma11_dth + y[4]*y[4]*metric3.d_gamma22_dth + y[5]*y[5]*metric3.d_gamma33_dth);
        equations[5] = 0;
        return equations;
    };
    auto stop_test3 = [&Rin,&Rout,&D,&biggestR,&tolerance,&r_H] (double x, const std::vector<double>& y) {return (std::abs(y[1]-M_PI/2.0)<tolerance && y[0]<Rout) || y[0]<1.01*r_H || y[0]>biggestR; };
    std::ofstream out("image_test.csv");
    std::vector<std::vector<double>> pixels(im_length, std::vector<double> (im_width));
    int y_index = 0;
 for(double y = -10; y < 10; y=y+n){
     int x_index = 0;
        for(double x = -25; x < 25; x=x+n){
            double rObs = std::sqrt(D*D + x*x + y*y);
            double beta = x/D;
            double alpha = y/D; 
            double th = th0-alpha;
            metric3.compute_metric(rObs, th);
            double u_r = - std::sqrt(metric3.g_11)*std::cos(beta)*std::cos(alpha); //initial photon momentum in r direction
            double u_th = std::sqrt(metric3.g_22)*std::sin(alpha); //initial photon momentum in theta direction
            double u_phi = std::sqrt(metric3.g_33)*std::sin(beta)*std::cos(alpha); //initial photon momentum in phi direction
            //implementing redshift 
            //from u_upper t = sqrt(g_ij*u_i*u_j)/alpha
            //from u_lower t = -alpha*ut_upper + u_i*beta_i 
            std::vector<double> y6 = {rObs,th,beta,u_r,u_th,u_phi};
            rk4.integrate(dydx3, stop_test3, h, x0, y6);
            double pixel = 0; 
            double final_radius = rk4.result[rk4.result.size()-1.0][0];
            double final_ur = -rk4.result[rk4.result.size()-1.0][3];
            double final_uPhi = -rk4.result[rk4.result.size()-1.0][5];
            double final_uTh = -rk4.result[rk4.result.size()-1.0][4];
            double ohm = 1.0/(a + (std::pow(final_radius, 3.0/2.0)/std::sqrt(M)));
            double final_angle = rk4.result[rk4.result.size()-1.0][1];
            //double ut_upper = std::sqrt(metric3.gamma11 * final_ur * final_ur + metric3.gamma22 * final_uTh * final_uTh + metric3.gamma33 * final_uPhi * final_uPhi)/alpha; 
            //double ut_upper = std::sqrt((metric3.gamma11*u_r*u_r) + (metric3.gamma22 * u_th * u_th) + (metric3.gamma33*u_phi*u_phi))/metric3.alpha;
            //double u_t = -alpha*alpha*ut_upper + final_uPhi*metric3.beta3;
            //double z = (1.0 + ohm*(final_uPhi/u_t))/std::sqrt(-metric3.g_00 - ohm*ohm*metric3.g_33 - 2*ohm*metric3.g_03)- 1.0; 
            //double inverse = 1.0/(1+z); 
            //double I = inverse*inverse*inverse;
            //we assume I_em is 1 for simplicity since it's constant 
            //std::cout << final_radius << std::endl; 
           if(final_radius<Rin){
                 pixel = 0.0; 
            }
            else if(final_radius > Rout){
                 pixel = 0.0;
            }
           else if(std::abs(final_angle-M_PI/2.0)< tolerance){
               double ut_upper = std::sqrt(metric3.gamma11 * final_ur * final_ur + metric3.gamma22 * final_uTh * final_uTh + metric3.gamma33 * final_uPhi * final_uPhi)/metric3.alpha;
               double u_t = -metric3.alpha*metric3.alpha*ut_upper + final_uPhi*metric3.beta3;
               double z = (1.0 + ohm*(final_uPhi/u_t))/std::sqrt(-metric3.g_00 - ohm*ohm*metric3.g_33 - 2*ohm*metric3.g_03)- 1.0; 
               double inverse = 1.0/(1+z); 
                double I = inverse*inverse*inverse;
               pixel = ((std::abs(I)));  
               if (pixel!=pixel){
                    std::cout<<ut_upper<<std::endl;
               }
            }
          pixels[y_index][x_index] = pixel;
          x_index ++;
          //std::cout<<i<<j<<std::endl;

            //out << std::sqrt(r_photon*r_photon+a*a)*std::sin(th)*std::cos(beta) << "," << std::sqrt(r_photon*r_photon+a*a)*std::sin(th)*sin(beta) << "," << r_photon*std::cos(th) << "," << pixel << "\n";

    }
     y_index ++;
    }
    std::ofstream output_file("raytracing_test.csv");
     for(int y_index =0; y_index<im_length;y_index++){
          output_file<<pixels[y_index][0];
          for(int x_index = 1; x_index<im_width;x_index++){
               output_file<<","<<pixels[y_index][x_index];
          }
          output_file<<std::endl;
     }
    //out.close(); 
    return 0;
}