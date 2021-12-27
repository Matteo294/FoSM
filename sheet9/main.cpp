#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

// possible speed-up: compute the potential only for the terms involved in the acceleration calculation

double calculate_potential();
double norm3d(vector<double> r);
vector<double> calculate_force(int i);
void init_pos();
void init_vel();
void leapfrog_step();

const double length_fact = 3.4e-10; // m
const double energy_fact = 1.65e-21; // J
const double mass_fact = 6.69e-26; // kg
const int N1d = 10;
const int Npart = N1d*N1d*N1d;
const double h = 1e-3; // numerical derivative step

const double m = 1.0; // m/mass_fact
const double T = 100.0; // K
const double k = 1.381e-23/energy_fact; // K^-1

vector<vector<double>> r(Npart);

int main(){

    init_pos();
    init_vel();

    return 0;
}

double calculate_potential(){
    double V = 0.0;
    vector<double> deltar (3, 0.0);
    double dist = 0.0;
    for(int i=0; i<r.size(); i++){
        for(int j=i+1; j<r.size(); j++){
            if (i != j){
                for(int k=0; k<3; k++) deltar[k] = r[i][k] - r[j][k];
                dist = norm3d(deltar);
                V += pow(length_fact/dist, 12) - pow(length_fact/dist, 6);
            } 
            
        }
    }
    return 4.0*V;
}

double norm3d(vector<double> r){
    return sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
}

vector<double> calculate_force(int i){
    double V1 = calculate_potential();
    double V2;
    vector<double> force(3, 0.0);
    for(int k=0; k<3; k++){
        r[i][k] += h;
        V2 = calculate_potential();
        force[k] = (double) (V2-V1)/h;
        r[i][k] -= h;
    }
    return force;
}

void init_pos(){
    for(int i=0; i<N1d; i++){
        for(int j=0; j<N1d; j++){
            for(int k=0; k<N1d; k++){
                //cout << Npart << " " << i+j*N1d+k*N1d*N1d << endl;
                r[i+j*N1d+k*N1d*N1d] = vector<double> {(double) 5.0*i, (double) 5.0*j, (double) 5.0*k};
            }
        }
    }
}


// Velocity initialized not in the asked way
void init_vel(){
    mt19937 rndgen;
    unsigned long seed = random_device{}();
    rndgen.seed(seed);
    std::normal_distribution<> boltzmann(0.0, sqrt(k*T/m));
    for(int i=0; i<Npart; i++){
            r[i].resize(3, 0.0); // initialize particle i in the origin
            for(int k=0; k<3; k++) r[i][k] = boltzmann(rndgen);
        }
}

void leapfrog_step(){
    
}