#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <omp.h>

#define RUN_PARALLEL 1

using namespace std;

double calculate_kinetic(); // Calculate total kinetic energy
double calculate_potential(); // Calculate total potential energy
double norm3d(vector<double> r); // Calculate norm of a 3D vector
vector<double> calculate_force(int i); // Calculate force on particle i
void init_pos(); // Initialize positions
void init_vel(); // Initialize velocities
void print_pos_to_file(ofstream &posfile); // Print updated positions to file

// Scaling units factors
const double xfact = 3.4e-10; // m (length)
const double efact = 1.65e-21; // J (energy)
const double mfact = 6.69e-26; // kg (mass)
const double vfact = sqrt(efact/mfact); // m/s (speed)

// Constants
const int N1d = 8; // Number of particles for dimension
const int Npart = N1d*N1d*N1d; // Total number of particles
const int Nsteps = (int) 60000; // Simulation steps
const double m = 1.0; // Mass of the particles
const double Temp = 80.0; // K
const double kb = 1.381e-23/efact; // Boltzmann constant (K^-1)
const double dt = 1e-2; // Integration step length
const double K_goal = 1.5*Npart*kb*Temp; // Thermal bath keeps kinetic energy constant

#if RUN_PARALLEL == 1
    const int Nthreads = 16; // Number of parallel threads
#else
    const int Nthreads = 1;
#endif 

vector<vector<double>> r(Npart); // Positions matrix
vector<vector<double>> v(Npart); // Velocities matrix

int main(int argc, char* argv[]){

    

    // Initialize files to store data
    ofstream datafile, posfile_before, posfile_after;
    datafile.open("data.csv");
    datafile << "K,V,E,T,Eerr" << endl;
    posfile_before.open("positions_before.csv");
    posfile_after.open("positions_after.csv");
    posfile_before << "x,y,z" << endl;
    posfile_after << "x,y,z" << endl;

    omp_set_num_threads(Nthreads); // Set number of parallel threads

    //vector<double> force(3, 0.0);
    double T, K, V, E;

    init_pos();
    init_vel(); 

    print_pos_to_file(posfile_before);

    K = calculate_kinetic();
    V = calculate_potential();
    double E0 = K + V;
    E = E0;
    T = K/(1.5*kb*Npart);

    // Integration cycles
    auto start = chrono::high_resolution_clock::now();
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::seconds>(stop-start).count();
    for(int i=0; i<Nsteps; i++){

        // Kick 1
        #pragma omp parallel for shared (v, dt)
        for (int n=0; n<Npart; n++){
            vector<double> force = calculate_force(n);
            for(int j=0; j<3; j++) v[n][j] += 0.5*dt*force[j]; // assuming mass 1
        }
        // Drift 1
        for (int n=0; n<Npart; n++){
            for(int j=0; j<3; j++){
                r[n][j] += dt*v[n][j];
                // Periodic boundary conditions
                //if (r[n][j] > 5*N1d) r[n][j] = r[n][j] - 5*N1d;
                //if (r[n][j] < 0) r[n][j] = 5*N1d + r[n][j];
                if ((r[n][j] > 5*N1d) || (r[n][j] < 0)) {r[n][j] -= dt*v[n][j]; v[n][j] = -v[n][j];}
            }
        }
        // Kick 2
        #pragma omp parallel for shared (v, dt)
        for (int n=0; n<Npart; n++){
            vector<double> force = calculate_force(n);
            for(int j=0; j<3; j++) v[n][j] += 0.5*dt*force[j]; // assuming mass 1
        }

        // Calculate things
        V = calculate_potential();
        K = calculate_kinetic();
        E = K+V;
        T = K/(1.5*kb*Npart);
        datafile << K << "," << V << "," << E << "," << T << "," << (E-E0)/E0 << endl;

        if (i % 1000 == 0){
            stop = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<std::chrono::seconds>(stop-start).count();
            cout << (int) 100*i/Nsteps << "%" << " Cycle time: " << duration << "s Energy: " << E << " Temperature: " << T << endl;
            start = chrono::high_resolution_clock::now();
        }

        // Thermal bath - comment this part if you want to work in the microcanonical ensemble
        if (i % 100 == 0){
            for(int j=0; j<Npart; j++){
                for (int k=0; k<3; k++) v[j][k] *= sqrt(Temp/T);
            }
        }
        // ----------------------------------------------------------------------------------
    }

    print_pos_to_file(posfile_after);


    return 0;
}

double calculate_kinetic(){
    double K = 0.0;
    for (int i=0; i<Npart; i++){
        for (int j=0; j<3; j++) K += 0.5*v[i][j]*v[i][j];
    }
    return K;
}

double calculate_potential(){
    double V = 0.0;
    double dist = 0.0;
    int j, k;
    #pragma omp parallel for reduction(+: V) private(dist) shared(r)
    for(int i=0; i<r.size(); i++){
        double locV = 0.0;
        for(int j=i+1; j<r.size(); j++){
            if (i != j){
                vector<double> deltar (3, 0.0);
                for(int k=0; k<3; k++) deltar[k] = r[i][k] - r[j][k];
                dist = norm3d(deltar);
                if (dist < 10.0) locV += pow(1.0/dist, 12) - pow(1.0/dist, 6);
            } 
            
        }
        V += locV;
    }
    return 4.0*V;
}

double norm3d(vector<double> r){
    return sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
}

vector<double> calculate_force(int i){
    vector<double> force (3, 0.0);
    vector<double> deltar (3, 0.0);
    double dist = 0.0;
    for(int j=0; j<r.size(); j++){
        if (i != j){
            for(int k=0; k<3; k++) deltar[k] = r[i][k] - r[j][k];
            dist = norm3d(deltar);
            for(int k=0; k<3; k++){
                if (dist < 10.0) force[k] += (double) 24.0 * (2*pow(1.0/dist, 13) - pow(1.0/dist, 7)) * deltar[k]/dist;
            }
        } 
        
    }
    return force;
}

void init_pos(){
    for(int i=0; i<N1d; i++){
        for(int j=0; j<N1d; j++){
            for(int k=0; k<N1d; k++){
                r[i+j*N1d+k*N1d*N1d] = vector<double> {(double) 5.0*i, (double) 5.0*j, (double) 5.0*k};
            }
        }
    }
}

void init_vel(){
    mt19937 rndgen[3];
    vector<unsigned long> seed {random_device{}(), random_device{}(), random_device{}()};
    rndgen[0].seed(seed[0]);
    rndgen[1].seed(seed[1]);
    rndgen[2].seed(seed[2]);
    std::normal_distribution<> boltzmann(0.0, sqrt(kb*Temp/m));
    for(int i=0; i<Npart; i++){
            v[i].resize(3, 0.0); // initialize particle i in the origin
            for(int k=0; k<3; k++) v[i][k] = boltzmann(rndgen[k]);
        }
}

void print_pos_to_file(ofstream &posfile){
    for(int i=0; i<Npart; i++){
        for(int j=0; j<2; j++){
            posfile << r[i][j] << ",";
        }
        posfile << r[i][2] << endl;
    }
}