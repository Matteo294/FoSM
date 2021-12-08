#include <iostream>
#include <cmath>
#include <fstream>
#include <random>

using namespace std;

int main(int argc, char** argv){

    mt19937_64 prng;
    unsigned long seed = std::random_device{}();
    prng.seed(seed);
    std::uniform_real_distribution<> dis(0.0, 1.0);

    ofstream RANDU_file, gen1_file, gen2_file;
    RANDU_file.open("RANDU_data.csv");
    gen1_file.open("gen1_data.csv");
    gen2_file.open("gen2_data.csv");
    gen1_file << "x,y" << endl;
    gen2_file << "x,y" << endl;
    RANDU_file << "x,y" << endl;
    int ncycles = stoi(argv[1]);

    long int I = (double) 2*time(NULL);
    cout << "Starting value: " << I  << " " << pow(2, 31) << endl;

    for(int i=0; i<ncycles; i++){
        if (i % 10000 == 0) cout << i << endl;
        // Randu
        I = fmod((65539*I) , pow(2, 31));
        RANDU_file << (double) I/pow(2, 31) << ",";
        I = fmod((65539*I) , pow(2, 31));
        RANDU_file << (double) I/pow(2,31) << endl;
        // Gen 1
        gen1_file << dis(prng) << "," << dis(prng) << endl;
        // Gen 2
        gen2_file << (double) rand()/RAND_MAX << "," << (double) rand()/RAND_MAX << endl;
    }
    return 0;

}