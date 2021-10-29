#include <iostream>
#include <iomanip>

using namespace std;

int main(){

    // Single precision
    float eps = 0.5;
    float old_eps;
    while((1+eps) > 1){
        old_eps = eps;
        eps /= 2.;
    }
    cout << "Single precision machine epsilon: " << old_eps << endl;
    cout << 1+old_eps << endl;

    // Double precision
    double eps2 = 0.5;
    double old_eps2;
    while((1+eps2) > 1){
        old_eps2 = eps2;
        eps2 /= 2.;
    }
    cout << "Double precision machine epsilon: " << old_eps2 << endl;
    cout << 1+old_eps2 << endl;

    // Long double precision
    long double eps3 = 0.5;
    long double old_eps3;
    while((1+eps3) > 1){
        old_eps3 = eps3;
        eps3 /= 2.;
    }
    cout << "Long double precision machine epsilon: " << old_eps3 << endl;
    cout << setprecision(30) << 1+old_eps3 << endl;
    
    return 0;
}