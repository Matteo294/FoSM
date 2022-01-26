#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

vector<double> HLL( vector<double> qL, vector<double> qR);
vector<double> getFlux(vector<double> q);
double get_dt(vector<vector<double>> q);

const double L = 1.0;
const int Nx = 1000;
const double rho0l=1.0, p0l=1.0, u0l=0.0;
const double rho0r=0.125, p0r=0.1, u0r=0.0;
const double adiab_idx=1.4;
const double tsim=0.2;
const double CFL=0.5;
const double dx = (double) L/Nx;

int main(){

    vector<vector<double>> F(Nx+2);
    vector<vector<double>> q(Nx+2);
    vector<double> x(Nx+2, 0.0);

    ofstream datafile;
    double dt;
    double t;
    vector<double> R1(3, 0.0), R2(3, 0.0);
    vector<vector<double>> R(Nx+2);
    for(int i=0; i<=Nx+1; i++){
        R[i].resize(3, 0.0);
    }

    // Create grid
    for(int i=0; i<Nx; i++) x[i+1] = (i+0.5)*dx; // Place point at the middle of the cell
    x[0] = -0.5*dx; x[Nx+1] = L+0.5*dx; // Just to give positions also to ghosts


    // Set initial fluxes and states and PBC
    q[0] = {rho0l, rho0l*u0l, p0l/(adiab_idx-1.)};
    q[Nx+1] = {rho0r, rho0r*u0r, p0r/(adiab_idx-1.)};
    F[0] = getFlux(q[0]);
    F[Nx+1] = getFlux(q[Nx+1]);
    for(int i=1; i<=Nx; i++){
        F[i].resize(3, 0.0);
        q[i].resize(3, 0.0);
        if (x[i] < 0.5) q[i] = q[0];
        else q[i] = q[Nx+1];
        F[i] = getFlux(q[i]);
    }

    // Time evolution
    t = 0.0;
    int n=0;
    while(t <= 0.2){
        
        cout << "t: " << t << " cycle: " << n << endl;

        // Write data to file (no ghosts)
        datafile.open("data/data_" + to_string(n) + ".csv");
        datafile << "x,rho,p,u" << endl;
        for(int i=1; i<=Nx+1; i++) datafile << x[i] << "," << q[i][0] << "," << (double) (adiab_idx-1.0) * (q[i][2] - (double) 0.5*q[i][1]/q[i][0]*q[i][1]) << "," << (double) q[i][1]/q[i][0] << endl;
        datafile.close();

        // Get delta_t through CFL
        dt = get_dt(q);

        // Apply PBC
        // --

        // Cycle through Nx interfaces
        for(int i=1; i<=Nx; i++){
            
            // Evaluate interfaces
            R2 = HLL(q[i], q[i+1]);
            R1 = HLL(q[i-1], q[i]);
            
            // Residuals
            for(int j=0; j<3; j++){
                R[i][j] = (double) (R2[j] - R1[j]) * dt/dx;
            }
        }

        // Update (RK1)
        for(int i=0; i<Nx+2; i++){
            for(int j=0; j<3; j++) q[i][j] -= R[i][j];
        }

        t+=dt;
        n++;
    }

}

vector<double> HLL( vector<double> qL, vector<double> qR){
    vector<double> FL = getFlux(qL);
    vector<double> FR = getFlux(qR);
    double CsL = (double) sqrt(adiab_idx*(adiab_idx-1) * (qL[2] - (double) 0.5*qL[1]/qL[0]*qL[1]) / qL[0]);
    double CsR = (double) sqrt(adiab_idx*(adiab_idx-1) * (qR[2] - (double) 0.5*qR[1]/qR[0]*qR[1]) / qR[0]);
    double SL = min((double) qL[1]/qL[0] - CsL, (double) qR[1]/qR[0] - CsR);
    double SR = max((double) qL[1]/qL[0] + CsL, (double) qR[1]/qR[0] + CsR);
        
    if(SL >= 0) return FL;
    else if (SR <= 0) return FR;
    else {
        vector<double> x(3, 0.0);
        for(int i=0; i<3; i++){
            x[i] = (double) (SR*FL[i] - SL*FR[i] + (double) SL*SR*(qR[i]-qL[i])) / (SR-SL);
        }
        return x;
    }
}

vector<double> getFlux(vector<double> q){
    vector<double> F(3, 0.0);
    double p = (double) (adiab_idx-1) * (q[2] - (double) 0.5*q[1]*q[1]/q[0]);
    double u = (double) q[1]/q[0];
    return vector<double> {q[1], (double) q[1]*u + p, u*(q[2] + p)};
}

double get_dt(vector<vector<double>> q){
    double vmax = 0;
    double v;
    for(auto val: q){
        v = abs(val[1]/val[0]) + sqrt((double) adiab_idx*(adiab_idx-1) * (val[2] - (double) 0.5*val[1]*val[1]/val[0]) / val[0]);
        if (v > vmax) vmax = v;
    }
    return (double) CFL * dx / vmax;
}