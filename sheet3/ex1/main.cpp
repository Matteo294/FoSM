#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

// Unit conversion
const long double Msun=1.99e30, xscale=1.496e11, vscale=2.98e4, tscale=xscale/vscale;
const int nsteps=5e4;
// Constants of the problem
const double Ms=1.0, Mp=1e-3, G=6.67408e-11*Msun/pow(xscale,3)*pow(tscale,2), dt=5e-4;
// Initial conditions
vector<double> xp {1.0, 0.0, 0.0}, vp {0.0, 0.5, 0.0};
vector<double> xs {0.0, 0.0, 0.0}, vs {0.0, 0.0, 0.0};
double E, E0;

// Store things
vector<double> force(3, 0.0); // on the planet

// Norm of a 3d vector
double norm3D(vector<double> x){
    double s=0.0;
    for(int i=0; i<3; i++){
        s += x[i]*x[i];
    }
    return sqrt(s);
}
// Calculate \ration (i.e. right hand side of the diff. equation)
void compute_force(){
    vector<double> diff(3, 0.0);
    for(int i=0; i<3; i++) diff[i] = xs[i] - xp[i];
    double r_cube = pow(norm3D(diff), 3);
    for(int i=0; i<3; i++){
        force[i] = G*Ms*Mp/r_cube*(xs[i]-xp[i]); 
    }
}
// Integrator
void kick_drift_kick(double dt){
    compute_force();
    for(int i=0; i<3; i++) 0.5*dt*force[i]/Mp;
    for(int i=0; i<3; i++){
        vp[i] += 0.5*dt*force[i]/Mp;
        vs[i] -= 0.5*dt*force[i]/Ms;
        xp[i] += vp[i]*dt;
        xs[i] += vs[i]*dt;
    }
    compute_force();
    for(int i=0; i<3; i++){
        vp[i] += 0.5*dt*force[i]/Mp;
        vs[i] -= 0.5*dt*force[i]/Ms;
    }
}
void rk2(double dt){
    vector<double> xp_old(3, 0.0), xs_old(3, 0.0), vp_old(3, 0.0), vs_old(3, 0.0);
    vector<double> k1_xs(3,0), k1_xp(3,0), k1_vs(3,0), k1_vp(3,0);
    compute_force();
    for(int i=0; i<3; i++){
        xp_old[i] = xp[i];
        xs_old[i] = xs[i];
        vp_old[i] = vp[i];
        vs_old[i] = vs[i];
        k1_vs[i] = -dt*force[i]/Ms;
        k1_vp[i] = dt*force[i]/Mp;
        k1_xs[i] = dt*vs[i];
        k1_xp[i] = dt*vp[i];
        xp[i] += k1_xp[i];
        xs[i] += k1_xs[i];
        vp[i] += k1_vp[i];
        vs[i] += k1_vs[i];
    }
    compute_force();
    for(int i=0; i<3; i++){
        xp[i] = xp_old[i] + 0.5*(k1_xp[i] + dt*vp[i]);
        vp[i] = vp_old[i] + 0.5*(k1_vp[i] + dt*force[i]/Mp);
        xs[i] = xs_old[i] + 0.5*(k1_xs[i] + dt*vs[i]);
        vs[i] = vs_old[i] + 0.5*(k1_vs[i] - dt*force[i]/Ms);
    }
}

double calculate_energy(){
    vector<double> diff(3, 0.0);
    for(int i=0; i<3; i++) diff[i] = xs[i] - xp[i];
    return -G*Mp*Ms/norm3D(diff) + 0.5*Mp*pow(norm3D(vp), 2) + 0.5*Ms*pow(norm3D(vs), 2);
}

int main(){

    ofstream outfile;
    outfile.open("data.csv");
    outfile << "x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2,energy,deltaE" << endl;
    outfile << xp[0] << "," << xp[1] << "," << xp[2] << "," << vp[0] << "," << vp[1] << "," << vp[2] << endl;

    E0 = calculate_energy();
    for(int i=0; i<nsteps; i++){
        // Select the integratio method by commeting the other line
        kick_drift_kick(dt);
        //rk2(dt);

        E = calculate_energy();
        outfile << xp[0] << "," << xp[1] << "," << xp[2] << "," << vp[0] << "," << vp[1] << "," << vp[2] << ","
                << xs[0] << "," << xs[1] << "," << xs[2] << "," << vs[0] << "," << vs[1] << "," << vs[2] << "," 
                << E << "," << abs(E-E0)/E0 << endl;
    }

    return 0;
}