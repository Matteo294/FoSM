#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

double func1D(double x);
double func(int dim);
void midpoint(int dim_tot, int dim);

const int N = 20000;
const int n = 6;
const int dmax = 10;
const int n_measures = 50;
const double h = (double) 1/n;
double val_midpoint;
vector<double> xvec(dmax, 0.0);

int main(){
    srand(time(NULL));
    auto begin = chrono::high_resolution_clock::now();
    auto end = chrono::high_resolution_clock::now();
    auto deltaT_midpoint = chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
    auto deltaT_MC = chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

    ofstream timesfile, errorsfile;
    timesfile.open("times.csv");
    timesfile << "midpoint,mc" << endl;
    errorsfile.open("errors.csv");
    errorsfile << "midpoint,mc" << endl;

    double x;
    double val_single; // value single run (for MC)
    double val_MC; // value averaged over runs (for MC)
    double f; // auxiliary variable

    for (int d = 1; d <= dmax; d++){

        cout << "###########################################" << endl;
        cout << "Dimension: " << d << endl << endl;
        cout << scientific << setprecision(5);

        for(int i=0; i<dmax; i++) xvec[i] = 0.0; // reset vector before integration

        // Midpoint (multiple measures just for time measurements)
        begin = chrono::high_resolution_clock::now();
        for (int n=0; n<n_measures; n++) {val_midpoint = 0.0; midpoint(d, d);} // integration
        end = chrono::high_resolution_clock::now();
        deltaT_midpoint = chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() / n_measures;
        cout << "Midpoint rule"  << endl;
        cout << "val: " << val_midpoint << endl;
        cout << "time: " <<  deltaT_midpoint << endl << endl;

        // Monte Carlo
        begin = chrono::high_resolution_clock::now();
        val_MC = 0.0; // final value (averaged)
        for(int n=0; n<n_measures; n++){
            val_single = 0.0; // value single run
            for (int i=0; i<N; i++){
                f = 1.0; // value of the function
                for (int dim=0; dim<d; dim++){
                    x = (double) rand()/RAND_MAX;
                    f *= func1D(x); // evaluate for each dimension
                }
                val_single += f;
            }
            val_single = (double) val_single / N;
            val_MC += val_single; // sum values single simulations to calculate average
        }
        end = chrono::high_resolution_clock::now();
        deltaT_MC = chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() / n_measures;
        val_MC = (double) val_MC / n_measures;
        cout << "Monte Carlo: " << endl;
        cout << "val: " << val_MC << endl;
        cout << "time: " << deltaT_MC << endl;
        cout << "###########################################" << endl << endl;

        timesfile << scientific << setprecision(5) << (double) deltaT_midpoint << "," << (double) deltaT_MC << endl;    
        errorsfile << scientific << setprecision(5) << abs(val_midpoint - 1.0) << "," << abs(val_MC - 1.0) << endl;
    }
}

double func1D(double x){
    return 3./2. * (1 - x*x);
}

double func(int dim){
    double f = 1.0;
    for(int i=0; i<dim; i++){
        f *= func1D(xvec[i]);
    }
    return f;
}

// algorithm https://github.com/janosh/numeric-simulations/blob/main/sheet%209/Monte%20Carlo/montecarlo.c
void midpoint(int dim_tot, int dim){
    if (dim == 0) val_midpoint += func(dim_tot)*pow(h, dim_tot);
    else{
        for(int i=0; i<n; i++){
            xvec[dim-1] = (double) i*h + 0.5*h;
            midpoint(dim_tot, dim-1);
        }
    }
}