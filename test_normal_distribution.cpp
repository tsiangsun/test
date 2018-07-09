// Test C++11 random number generator (Normal distribution)
// compile: g++ -std=c++11 test_gaussian_random_c++11.cpp -o test_gaussian_random_c++11

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <random>
using namespace std;

//CONSTANTS

const int N = 5;
const int LEN = 10;

//DECLARE SUBROUTINES

int test_ran(vector<double>& x, int n, int seed);
int test_ran1(vector<double>& x, int n, mt19937& gen1);

//MAIN PROGRAM
int main (int argc, char *argv[]) {
    
    int id(0); //job number
    stringstream ss;
    string emptystr("");
    string idstr("");
    
    if (argc > 1) {
        ss << argv[1];
        idstr += ss.str();
        ss >> id;
        ss.clear();
    }
    cout << "-------------------------" << endl;
    cout << ">>> Job id # " << id << endl;
    
    int itraj; //index of trajectory
    int i,j,k;
    
    vector<double> x(LEN,0);
    
    mt19937 gen;
    gen.seed(id+15);

    for (itraj = 0; itraj < N; itraj++) {
        cout << " ---> Traj # " << itraj << ":" << endl;
        
        //test_ran(x, LEN, itraj);
        test_ran1(x, LEN, gen);
            
        for (i=0 ; i < LEN; i++) {
            cout << x[i] << endl;
        }
        
    }
    
    

    return 0;
}





// --------------------------------------------------
// ------------------ SUBROUTINES -------------------
// --------------------------------------------------

int test_ran(vector<double>& x, int n, int seed) {
    
    //true random seed based on system time
    //random_device rd;
    //mt19937 gen(rd());
    
    //manual random seed, reproduceble, id number as seed
    mt19937 gen;
    gen.seed(seed);
    normal_distribution<double> normal_dist(0.0,1.0);
    
    for (int i = 0 ; i < n ; i++) {
        x[i] = normal_dist(gen);
    }
    return 0;
}


int test_ran1(vector<double>& x, int n, mt19937& gen1) {
    
    //true random seed based on system time
    //random_device rd;
    //mt19937 gen(rd());
    
    //manual random seed, reproduceble, id number as seed
    //mt19937 gen;
    //gen.seed(myseed);
    normal_distribution<double> normal_dist(0.0,1.0);
    
    for (int i = 0 ; i < n ; i++) {
        x[i] = normal_dist(gen1);
    }
    return 0;
}





