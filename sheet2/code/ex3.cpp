#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

// Declaration of functions. Functions are implemented after the main function
vector<double> RHS(vector<double> y, vector<double> phidot); // Evaluate right hand side of equation ydot = f(y)
vector<double> integrate_euler(double dt, vector<double> y_old, vector<double> phidot);
vector<double> integrate_RK2(double dt, vector<double> y_old, vector<double> phidot); // Makes one integration step: depends on the chosen integration scheme
vector<double> get_phidot(vector<double> y); // velocity as a function of conj momenta
double evaluate_energy(vector<double> y, vector<double> phidot);

double m1=1.1, m2=1.1, l1=1.1, l2=1.1, g=1.1, mu=m1*m2/(m1+m2);

int main(){

    // File to print data for plotting
    ofstream datafile;
    datafile.open("data.csv");
    datafile << "t,phi1,phi2,phi1dot,phi2dot,p1,p2,deltaE" << endl;
    // File for animation
    ofstream animfile;
    animfile.open("animation.csv");
    animfile << "x1,y1,x2,y2" << endl;

    // Initial conditions
    const double phi1_0 = (double) M_PI/8.;
    const double phi2_0 = (double) M_PI/3.;
    const double p1_0 = 0;
    const double p2_0 = 0;

    // Parameters
    const int nsteps = 2000; // Number of integration steps
    const double dt = 0.05; // time step
    const double g = 1; // Gravity

    // State variables
    double t;
    vector<double> y(4, 0);
    vector<double> phidot(2, 0);
    double energy, energy0;
    // Cartesian coordintes
    double x1, x2, y1, y2;

    // Assign initial conditions
    t = 0.0;
    y[0] = phi1_0;
    y[1] = phi2_0;
    y[2] = p1_0;
    y[3] = p2_0;
    phidot = get_phidot(y);

    energy = evaluate_energy(y, phidot);
    energy0 = energy;
    cout << energy0 << endl;
    // Integration
    for(int i=0; i<nsteps; i++){
        phidot = get_phidot(y);
        y = integrate_RK2(dt, y, phidot);
        t += dt;
        energy = evaluate_energy(y, phidot);
        datafile << t << "," << y[0] << "," << y[1] << "," << phidot[0] << "," << phidot[1] << "," << y[2] << "," << y[3] << "," << (energy-energy0)/energy0 << endl;

        // Cartesian coordinates for the animation
        x1 = l1*sin(y[0]);
        y1 = -l1*cos(y[0]);
        x2 = x1 + l2*sin(y[2]);
        y2 = y1 - l2*cos(y[2]);
        animfile << x1 << "," << y1 << "," << x2 << "," << y2 << endl;
    }

    return 0;
}

vector<double> RHS(vector<double> y, vector<double> phidot){
    vector<double> f(4, 0);
    f[0] = phidot[0];
    f[1] = phidot[1];
    f[2] = -m2*l1*l2*phidot[0]*phidot[1]*sin(y[0]-y[1]) - (m1+m2)*g*l1*sin(y[0]);
    f[3] = m2*l1*l2*phidot[0]*phidot[1]*sin(y[0]-y[1]) - m2*g*l2*sin(y[1]);
    return f;
}

vector<double> integrate_euler(double dt, vector<double> y_old, vector<double> phidot){
    vector<double> ynew(4);
    ynew = RHS(y_old, phidot);
    for(int i=0; i<4; i++){
        ynew[i] = y_old[i] + dt*ynew[i];
    }
    return ynew;
}

vector<double> integrate_RK2(double dt, vector<double> y_old, vector<double> phidot){
    vector<double> aux(4); // Auxiliary variable
    vector<double> k1(4);
    vector<double> ynew(4);
    vector<double> phidot_new(4);
    k1 = RHS(y_old, phidot);
    for(int i=0; i<4; i++){
        aux[i] = y_old[i] + 0.5*dt*k1[i];
    }
    phidot_new = get_phidot(aux);
    aux = RHS(aux, phidot_new); // k2
    for(int i=0; i<4;i++){
        ynew[i] = y_old[i] + dt*aux[i];
    }
    return ynew;
}

vector<double> get_phidot(vector<double> y){
    vector<double> phidot(2, 0);
    double fact = m1*l1*(1 + m2/m1*pow(sin(y[0]-y[1]), 2));
    phidot[0] = y[2]/l1/fact - cos(y[0]-y[1])*y[3]/fact/l2;
    phidot[1] = -cos(y[0]-y[1])*y[2]/l2/fact - y[3]*(m1+m2)/m2/l2/fact;
    return phidot;
}

double evaluate_energy(vector<double> y, vector<double> phidot){
    return 0.5*m1*pow(l1*phidot[0], 2) * 0.5*m2*(pow(l1*phidot[0], 2) + pow(l2*phidot[1], 2) + 2*l1*l2*phidot[0]*phidot[1]*cos(y[0]-y[1])) + 
    m1*g*l1*(1-cos(y[0])) + m2*g*(l1*(1-cos(y[0])) + l2*(1-cos(y[1])));
}
