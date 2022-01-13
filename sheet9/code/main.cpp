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
double Temp; int N1d, Nsteps, Npart; // Pass by terminal, see .sh files

// Constants
const double m = 1.0; // Mass of the particles
const double kb = 1.381e-23/efact; // Boltzmann constant (K^-1)
const double dt = 1e-2; // Integration step length

#if RUN_PARALLEL == 1
    const int Nthreads = 16; // Number of parallel threads
#else
    const int Nthreads = 1;
#endif 

vector<vector<double>> r; // Positions matrix
vector<vector<double>> v; // Velocities matrix

int main(int argc, const char* argv[]){

    if (argc != 5) cout << "Please see the .sh files to understand how to run the simulations" << endl;

    const int mode = atoi(argv[1]); // 0 for microcanonical, 1 for canonical (pass by terminal when running)
    Temp = atoi(argv[2]);
    N1d = atoi(argv[3]);
    Nsteps = atoi(argv[4]);
    Npart = N1d*N1d*N1d;

    r.resize(Npart);
    v.resize(Npart);

    if (mode == 0) cout << "Microcanonical ensemble with energy: " << 1.5*kb*Npart*Temp << endl;
    else cout << "Canonical ensemble with temperature: " << Temp << " K" << endl;
    
    double T, K, V, E, avgE, avgT, avgK, avgV;

    // Variables to store time
    auto start = chrono::high_resolution_clock::now();
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::seconds>(stop-start).count();

    // Initialize files to store data
    ofstream datafile, posfile_before, posfile_after;
    datafile.open("data.csv");
    datafile << "K,V,E,T,Eerr" << endl;
    posfile_before.open("positions_before.csv");
    posfile_after.open("positions_after.csv");
    posfile_before << "x,y,z" << endl;
    posfile_after << "x,y,z" << endl; 

    omp_set_num_threads(Nthreads); // Set number of parallel threads

    init_pos();
    init_vel(); 

    print_pos_to_file(posfile_before);

    K = calculate_kinetic();
    V = calculate_potential();
    double E0 = K + V;
    E = E0;
    T = K/(1.5*kb*Npart);

    avgE = 0.0;
    avgK = 0.0;
    avgT = 0.0;
    avgV = 0.0;
    // Integration cycles
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
                while (r[n][j] >= 5*N1d) r[n][j] -= (double) 5*N1d;
                while (r[n][j] < 0) r[n][j] += (double) 5*N1d;
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

        datafile << K/Npart << "," << V/Npart << "," << E/Npart << "," << T << "," << (E-E0)/E0 << endl;

        if (i % 1000 == 0){
            stop = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<std::chrono::seconds>(stop-start).count();
            cout << (int) 100*i/Nsteps << "%" << " Cycle time: " << duration << "s Energy: " << E << " Temperature: " << T << endl;
            start = chrono::high_resolution_clock::now();
        }

        // Thermal bath
        if (mode == 1){
            if (i % 100 == 0){
                for(int j=0; j<Npart; j++){
                    for (int k=0; k<3; k++) v[j][k] *= sqrt(Temp/T);
                }
            }
        }

        // Average energy
        if (i >= Nsteps - 10000) {avgE += E; avgK += K; avgV += V; avgT += T;}
    }

    cout << "Average energy: " << avgE/10000. << " Average kinetic: " << avgK/10000. << " Average potential: " << avgV/10000. << " Average temperature: " << avgT/10000. << endl << endl;

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
            for(int k=0; k<3; k++){
                deltar[k] = r[i][k] - r[j][k];
                if (deltar[k] > 0.5*5.0*N1d) deltar[k] -= 5.0*N1d;
                else if (deltar[k] < -0.5*5.0*N1d) deltar[k] += 5.0*N1d;
            }
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