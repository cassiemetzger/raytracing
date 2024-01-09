#include<iostream> 
#include <cmath>
class boyer_lindquist_metric {
public:
boyer_lindquist_metric(double a0, double M0) {
    a = a0; //a is spin parameter, between zero and 1
    M = M0; //mass of black hole
}
void compute_metric(double r, double th) {
    // in general, t is indexed at 0
    // r is indexed at 1
    // theta is indexed at 2
    // phi is indexed at 3
    //weird matrix stuff ask wikipedia not me
    //trig fxns
    double sine = std::sin(th);
    double cosine = std::cos(th);
    //rho,delta,sigma
    double rho_squared = (r*r) + (a*a)*cosine*cosine; 
    double delta = (r*r) + (a*a)-2.0*M*r;
    double sigma = ((r*r)+(a*a))*((r*r)+(a*a))-(a*a)*delta*sine*sine;
    //metric fxns
    g_00 = ((2.0*M*r)/rho_squared) -1.0; //gtt
    g_03 =-((2.0*M*a*r)/rho_squared)* sine*sine; //gtphi
    g_11 = rho_squared/delta; //grr
    g_22 =rho_squared; //gthth
    g_33 =(sigma/rho_squared)*sine*sine; //gphiphi
    //alpha, beta
    alpha = std::sqrt(rho_squared*delta/sigma);
    beta3= (-2*a*M*r)/sigma;
    //inverse metric fxns
    gamma11 = 1.0/g_11;
    gamma22 =1.0/g_22;
    gamma33 = 1.0/g_33;
    //derivatives
    double d_rho_squared_dr = 2.0*r;
    double d_delta_dr = (2.0*r)-(2.0*M);
    double d_sigma_dr = (4.0*r)*(r*r+a*a)-d_delta_dr*a*a*sine*sine;
    double d_rho_squared_dth = -a*a*std::sin(2.0*th);
    double d_delta_dth = 0.0;
    double d_sigma_dth = -a*a*delta*std::sin(2.0*th);
    double z = rho_squared*delta; //intermediate derivative
    double d_z_dr = (d_rho_squared_dr*delta)+(rho_squared*d_delta_dr);
    double d_z_dth = (d_rho_squared_dth*delta)+(rho_squared*d_delta_dth);
    //Use numerical derivatives
    //Mathmatica or wolframalpha
    //Go through the derivatives and check
    d_alpha_dr =(1.0/(2.0*alpha))*(sigma*d_z_dr-z*d_sigma_dr)/(sigma*sigma);
    d_alpha_dth =(1.0/(2.0*alpha))*(sigma*d_z_dth-z*d_sigma_dth)/(sigma*sigma);
    d_beta3_dr = (-2.0*M*a)*((sigma-r*d_sigma_dr)/(sigma*sigma));
    d_beta3_dth = (-2.0*M*a*r)*(-1.0/(sigma*sigma))*(d_sigma_dth);
    d_gamma11_dr = (rho_squared * d_delta_dr - delta * d_rho_squared_dr) / (rho_squared * rho_squared);
    d_gamma11_dth = (rho_squared * d_delta_dth - delta * d_rho_squared_dth) / (rho_squared * rho_squared);
    d_gamma22_dr = (-1.0/(rho_squared * rho_squared)) * d_rho_squared_dr;
    d_gamma22_dth = (-1.0/(rho_squared * rho_squared)) * d_rho_squared_dth;
    d_gamma33_dr = ((sigma*d_rho_squared_dr)-(rho_squared*d_sigma_dr))/(sigma*sigma*sine*sine);
    d_gamma33_dth = ((sigma * sine * sine * d_rho_squared_dth) - (rho_squared*((d_sigma_dth* sine * sine) + (sigma*std::sin(2.0*th)))))/(sigma*sigma*std::pow(sine, 4));
    
}
double a, M;
double alpha, beta3;
double gamma11, gamma22, gamma33; // components of upper gamma^ij
double g_00, g_03, g_11, g_22, g_33; // components of lower g_\mu\nu
double d_alpha_dr, d_beta3_dr, d_gamma11_dr, d_gamma22_dr, d_gamma33_dr;
double d_alpha_dth, d_beta3_dth, d_gamma11_dth, d_gamma22_dth, d_gamma33_dth;
};