#include <iostream>
#include <math.h>
#include <iomanip>

using namespace std;

double f(double x);
double f2(double x);
long int factorial(long int n);

int main(){

    double x;

    cout << "Point in which to evaluate the function: " << endl;
    cin >> x;
    cout << "Value of the function in x=" << x << " is: " << f(x) << endl;
    cout << "Value of x + exp(-x) is: " << fixed << setprecision(15) << x + exp(-x) << endl;

    cout << "Press x to continue" << endl;
    while (1){
        if ('x' == getchar())
        break;
    }   
    double expmax= -4;
    double expmin = -8;
    double delta_exp = (expmax-expmin) / 100.;
    double exp_x = expmax;
    x = pow(10, exp_x);
    for(int i=0; i<100; i++){
        cout << "x: " << x << "\t f(x): " << f(x) << "\t f2(x): " << f2(x) << setprecision(10) << "\t exp(-x): " << exp(-x)<< endl;
        exp_x -= delta_exp;
        x = pow(10, exp_x);
    }
    return 0;
}

double f(double x){
    return (x + exp(-x) - 1) / pow(x,2);
}

double f2(double x){
    if(x > 1e-5) return f(x);
    else{
        double s = 0.5;
        for(long int i=1; i<29; i++){
            s += (double) pow(-1, i) * pow(x, i) / factorial(i+1);
        }
        return s;
    }
}

long int factorial(long int n){
    if (n == 0) return 1;
    else return n*factorial(n-1);
}