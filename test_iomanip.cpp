// Test iomanip format output
// compile: g++ -std=c++11 test_iomanip.cpp -o test_iomanip

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

const int N = 2;
const int LEN = 21;

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

    cout.precision(6);
    
    for (itraj = 0; itraj < N; itraj++) {
        cout << "---> Traj # " << itraj << ":" << endl;
        
        //test_ran(x, LEN, itraj);
        test_ran1(x, LEN, gen);
            
        for (i=0 ; i < LEN; i++) {
            cout << setw(12) << x[i] << "*";
            
            cout << setw(15) << setprecision(9) << x[i] << endl;
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





