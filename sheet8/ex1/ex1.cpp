#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <iomanip>
#include <chrono>

using namespace std;

// Performs multiplication between the tridiagonal matrix and the vector
vector<double> tridiagonal_multiplication(vector<vector<double>> M, vector<double> T);

// Build the tridiagonal matrix with Dirichlet boundary conditions
vector<vector<double>> build_tridiagonal(int n);

// Build the matrix of 3 vectors starting from the tridiagonal matrix
vector<vector<double>> build_trivector(vector<vector<double>> A);

// Gauss-Jordan method (forward elimination backward subsitution)
vector<double> gauss_jordan(vector<vector<double>> X, vector<double> b);

// Print 1D vector
void print_vector(vector<double> x);

// Print 2D vector (override previous)
void print_vector(vector<vector<double>> X);

const double D = 0.5; // Difussion coefficient
const double eps = 1.0; // RHS member of the equation
const double T0 = 1.0; // Boundary conditions
const double L = 1.0; // Grid length
const int N = 100; // Grid points
const double h = 1e-2; // Integration step

int main(){

    vector<vector<double>> A = build_tridiagonal(N); // prepare matrix A such that Ax = b
    vector<vector<double>> M = build_trivector(A); // transform tridiagonal matrix A in a matrix M of three vectors containing the diagonals
    // Prepare RHS of the matrix equation
    vector<double> b(N, 0.0);
    for(int i=1; i<N-1; i++) b[i] = -h*h/D*eps;
    b[0] = T0;
    b[N-1] = T0;

    // Some printing
    cout << "A:" << endl;
    print_vector(A);
    cout << "M: "<< endl;
    print_vector(M);

    // Search solution x with the Gauss-Jordan algorithm
    auto begin = chrono::high_resolution_clock::now();
    vector<double> sol = gauss_jordan(A, b);
    auto end = chrono::high_resolution_clock::now();
    auto deltaT_midpoint = chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
    cout << "Solution found. Time taken: " << scientific  << setprecision(3) << deltaT_midpoint*1e-9 << defaultfloat << "s. Solution:" << endl;
    print_vector(sol);

    // Check if solution x is correct by calculating residuals Ax - b
    vector<double> estimated_b = tridiagonal_multiplication(M, sol);
    cout << "Residuals: " << endl;
    for(int i=0; i<N; i++) cout << estimated_b[i] - b[i] << " ";
    cout << endl;

    return 0;
}

// Performs multiplication between the tridiagonal matrix and the vector
vector<double> tridiagonal_multiplication(vector<vector<double>> M, vector<double> T){
    uint n = T.size();
    vector<double> res(n, 0.0);
    for(uint i=1; i<n-1; i++){
        res[i] = M[i][0]*T[i-1] + M[i][1]*T[i] + M[i][2]*T[i+1];
    }
    res[0] = M[0][0]*T[n-1] + M[0][1]*T[0] + M[0][2]*T[1];
    res[n-1] = M[n-1][0]*T[n-2] + M[n-1][1]*T[n-1] + M[n-1][2]*T[0];
    return res;
}

// Build the tridiagonal matrix with Dirichlet boundary conditions
vector<vector<double>> build_tridiagonal(int n){
    vector<vector<double>> A(n); // Creates one column
    for(uint i=1; i<n-1; i++){
        A[i].resize(n, 0.0); // Creates other elements of the row
        // Fill row
        A[i][i-1] = 1.;
        A[i][i] = -2.;
        A[i][i+1] = 1.;
    }
    // Fix first and last rows separately
    A[0].resize(n, 0.0);
    A[0][0] = 1.0;
    A[n-1].resize(n, 0.0);
    A[n-1][n-1] = 1.0;
    return A;
}

// Build the matrix of 3 vectors starting from the tridiagonal matrix
vector<vector<double>> build_trivector(vector<vector<double>> A){
    int n = A[0].size();
    vector<vector<double>> M(n);
    for(uint i=1; i<n-1; i++){
        M[i].resize(3, 0.0);
        M[i][0] = A[i][i-1];
        M[i][1] = A[i][i];
        M[i][2] = A[i][i+1];
    }
    // Fix first and last rows separately
    M[0].resize(3, 0.0);
    M[0][0] = A[0][n-1];
    M[0][1] = A[0][0];
    M[0][2] = A[0][1];
    M[n-1].resize(3, 0.0);
    M[n-1][0] = A[n-1][n-2];
    M[n-1][1] = A[n-1][n-1];
    M[n-1][2] = A[n-1][0];
    return M;
}

// Gauss-Jordan method (forward elimination backward subsitution)
vector<double> gauss_jordan(vector<vector<double>> X, vector<double> b){
    int n = X[0].size();
    double c; // auxiliary variable useful when swapping rows
    vector<double> sol(n, 0.0); // solution vector

    vector<double> indices;
    for(int i=0; i<N; i++){
        indices.push_back(i);
    }

    // Elimination
    for(uint i=1; i<n; i++){
        for(uint j=0; j<i; j++){
            if (X[i][j] != 0){
                c = (double) -X[i][j]/X[j][j];
                for(uint k=j; k<n; k++){
                    X[i][k] += (double) c*X[j][k];
                    //print_vector(X);
                }
                b[i] += (double) c*b[j];
            }
        }
        // If this row is in the wrong position (i.e. element on the diagonal is 0) put it at the end and redo with the next row
        if (X[i][i] == 0){
            swap_ranges( (X.begin()+i)->begin(), (X.begin()+i)->end(), (X.end()-1)->begin() ); // swap this row and the last one
            c = b[i];
            b[i] = b[n-1];
            b[n-1] = c;
            indices[i] = indices[n-1];
            indices[n-1] = indices[i];
            i -= 1; // re-do this row index with the new row
        }
    }

    // Substitution
    for(int i=N-1; i>=0; i--){
        sol[i] = b[i]/X[i][i];
        for(uint j=i+1; j<n; j++){
            sol[i] -= X[i][j]/X[i][i]*sol[j];
        }
    }

    return sol;
}

// Print 1D vector
void print_vector(vector<double> x){
    for(uint i=0; i<x.size(); i++){
        cout << x[i] << " ";
    }
    cout << endl << endl;
}

// Print 2D vector (overrides 1D)
void print_vector(vector<vector<double>> X){
    for(uint i=0; i<X.size(); i++){
        for(uint j=0; j<X[i].size(); j++){
            cout << X[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}