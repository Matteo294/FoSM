#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

void RHS(); // right hand side function of the diff. eq. automatically store the results in the matrix f
double calculate_energy();
void euler_step(double dt); // perform an euler integration step with step length dt
void RK2_step(double dt); // perform an RKS predictor-corrector integration step with step length dt
void RK4_step(double dt); // perform a RK4 integration step with step length dt

// Global variables
const double m1=0.5, m2=1.0, l1=2.0, l2=1.0, g=1.0; // physical constants of the problem
const double dt = 0.05; // integration step
vector<double> y(4, 0.0); // state variables phi1, phi2, q1, q2
vector<double> f(4, 0.0); // matrix to store the value of the RHS function
double C, c1, c2, c3, c4; // useful to evaluate RHS
const int nsteps = 2000;

int main(){

    // Initial conditions
    y[0] = 51*M_PI/180.;
    y[1] = -120*M_PI/180.;
    y[2] = 0.;
    y[3] = 0.;

    // Two files to store results: one for plotting (state variables) and one for the animation (cartesian coordinates)
    ofstream datafile;
    datafile.open("data_newIC.csv");
    datafile << "t,phi1,phi2,phi1dot,phi2dot,deltaE" << endl;
    ofstream animfile;
    animfile.open("animation_newIC.csv");
    animfile << "x1,x2,y1,y2" << endl;

    double t=0.0;
    double E0 = calculate_energy();
    double x1, x2, y1, y2;
    for(int i=0; i<nsteps; i++){
        RK4_step(dt);        
        t += dt;
       // Printing to files
        datafile << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << (calculate_energy() - E0)/E0 << endl;
        x1 = l1*sin(y[0]);
        x2 = x1 + l2*sin(y[1]);
        y1 = - l1*cos(y[0]);
        y2 = y1 - l2*cos(y[1]);
        animfile << x1 << "," << x2 << "," << y1 << "," << y2 << endl;
    }

    return 0;
}

void RHS(){
    C = m2 * l1 * l1 * l2 * l2 * (m1 + m2 * pow(sin(y[0] - y[1]), 2));
    c1 = m2 * l2 * l2 / C;
    c2 = -m2 * l1 * l2 * cos(y[0] - y[1]) / C;
    c3 = c2;
    c4 = (m1 + m2) * l1 * l1 / C;
    f[0] = c1 * y[2] + c3 * y[3]; // phi1dot, can be recycled and re-used in the evaluation of qdot
    f[1] = c3 * y[2] + c4 * y[3]; // phi2dot, can be recycled and re-used in the evaluation of qdot
    f[2] = -m2 * l1 * l2 * f[0] * f[1] * sin(y[0] - y[1]) - (m1 + m2) * g * l1 * sin(y[0]);
    f[3] = m2 * l1 * l2  * f[0] * f[1] * sin(y[0] - y[1]) - m2 * g * l2 * sin(y[1]);


}

double calculate_energy(){
    //return 0.5*m1*pow(l1*f[0], 2) + 0.5*m2*(pow(l1*f[0], 2) + pow(l2*f[1], 2) + 2*l1*l2*f[0]*f[1]*cos(y[0]-y[1])) + (m1+m2)*g*l1*(1 - cos(y[0])) + l2*(1 - cos(y[1]));
    return y[2] * f[0] + y[3] * f[1] - (m1 + m2)/2 * l1 * l1 * f[0] * f[0] - m2/2 * l2 * l2 * f[1] * f[1] - m2 * l1 * l2 * f[0] * f[1] * cos(y[0] - y[1]) + (m1 + m2) * g * l1 * (1 - cos(y[0])) + l2 * (1 - cos(y[1]));
}

void euler_step(double dt){
    RHS();
    for(int i=0; i<4; i++){
        y[i] += dt*f[i];
    }
}   

void RK2_step(double dt){
    vector<double> y_old(3, 0.0);
    vector<double> f_old(3, 0.0);
    for(int i=0; i<4; i++){
        f_old[i] = f[i];
        y_old[i] = y[i];
        y[i] += dt*f[i];
    }
    RHS();
    for(int i=0; i<4; i++){
        y[i] = y_old[i] + 0.5*dt*(f_old[i] + f[i]);
    }
}

void RK4_step(double dt){
    RHS();
    vector<double> k1(4, 0.0), k2(4, 0.0), k3(4, 0.0), k4(4, 0.0), aux(4,0.0);
    for(int i = 0; i < 4; i++) {
        k1[i] = dt * f[i];
        aux[i] = y[i];
        y[i] = y[i] + k1[i]/2.;
    }
    RHS();
    for(int i = 0; i < 4; i++) {
        k2[i] = dt * f[i];
        y[i] = aux[i] + k2[i]/2.;
    }
    RHS();
    for(int i = 0; i < 4; i++) {
        k3[i] = dt * f[i];
        y[i] = aux[i] + k3[i];
    }
    RHS();
    for(int i = 0; i < 4; i++) {
        k4[i] = dt * f[i];
        y[i] = aux[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.;
    }
}
