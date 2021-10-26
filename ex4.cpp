#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>

using namespace std;

double f(double x); // Normal function as given in the text of the exercise
double f2(double x); // Modified function to overcome numerical errors
long int factorial(long int n);

int main(){

    // File to print the value of the function to find out where the normal function stops yielding the correct result
    ofstream myfile;
    myfile.open("outputdata.csv");
    myfile << "x,f1(x),f2(x)" << endl;
   
    double x;

    // Evaluate the function in a given point
    cout << "Point in which to evaluate the function: " << endl;
    cin >> x;
    cout << "Value of the function in x=" << x << " is: " << f(x) << endl;

    cout << "Press x to continue" << endl;
    while (1){
        if ('x' == getchar())
        break;
    }  

    // First create a log-spaced array from 10^-4 to 10^-8 (i.e. create a lin space vector of exponents and then take the exponential of those values)
    double expmax= -4;
    double expmin = -8;
    double delta_exp = (expmax-expmin) / 1000.;
    double exp_x = expmax;
    x = pow(10, exp_x);

    // For each point evaluate the two implementations of the functions and compare them
    for(int i=0; i<1000; i++){
        cout << "x: " << x << "\t f(x): " << f(x) << "\t f2(x): " << f2(x) << endl;
        myfile << x << "," << f(x) << "," << f2(x) << endl;
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