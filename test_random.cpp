#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <random>
//#include "r_1279.h"
using namespace std;

//CONSTANTS

const int N = 5;

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
    
    /*
    long seed;
    seed = seedgen();	// have seedgen compute a random seed
    setr1279(seed);		//seed the genertor
    r1279(); //[0,1] random number
    */
    
    
    int itraj; //index of trajectory
    int i,j,k;
    
    vector<double> x(N,0);
    
    mt19937 gen;
    gen.seed(id+5);

    for (itraj = 0; itraj < N; itraj++) {
        cout << " ---> Traj # " << itraj << ":" << endl;
        
        //test_ran(x, N, itraj);
        test_ran1(x, N, gen);
            
        for (i=0 ; i < N; i++) {
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
    //gen.seed(seed);
    normal_distribution<double> normal_dist(0.0,1.0);
    
    for (int i = 0 ; i < n ; i++) {
        x[i] = normal_dist(gen1);
    }
    return 0;
}





