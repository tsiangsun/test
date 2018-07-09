//  Test cout format output to file
//

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <random>
#include <complex>
using namespace std;


typedef std::complex<double> Complex;
typedef std::vector<vector<Complex> > Complex_Matrix;
typedef std::vector<vector<double> > Real_Matrix;

const int LEN_TRAJ = 500;
const int DOFe = 2;

//MAIN PROGRAM
int main (int argc, char *argv[]) {

    int i,j,k;
    int flag(0);
    int id(0); //job number
    int itraj; //index of trajectory
    
    double a;
    
    stringstream ss;
    string emptystr("");
    string idstr("");
    string bcf_ind("");
    int label;
    
    ifstream infile;
    ofstream outfile;
    
    //BEGIN read in data
    cout << ">>> Converting format" << endl;
    
    Complex_Matrix sigma_D(DOFe*DOFe, vector<Complex>(LEN_TRAJ,0.0));

    
    //loading data
    infile.open((emptystr + "sigma_Kelly_D_init_b10_g0.5_w0.5_s0.dat").c_str());
    for (i = 0 ; i < LEN_TRAJ; i++) {
        infile >> a;   sigma_D[0][i].real(a);
        infile >> a;   sigma_D[0][i].imag(a);
        infile >> a;   sigma_D[1][i].real(a);
        infile >> a;   sigma_D[1][i].imag(a);
        infile >> a;   sigma_D[2][i].real(a);
        infile >> a;   sigma_D[2][i].imag(a);
        infile >> a;   sigma_D[3][i].real(a);
        infile >> a;   sigma_D[3][i].imag(a);
    }
    infile.close();
    infile.clear();
    
    
    //outputing data
    outfile.open((emptystr + "sigma_Kelly_D_init_b10_g0.5_w0.5_s0_reformat.dat").c_str());
    outfile << setprecision(8);
    for (i = 0 ; i < LEN_TRAJ; i++) {
        outfile << setw(19) << sigma_D[0][i].real() << setw(19) <<  sigma_D[0][i].imag()
        << setw(19)<<  sigma_D[1][i].real() << setw(19)<<  sigma_D[1][i].imag()
        << setw(19)<<  sigma_D[2][i].real() << setw(19)<<  sigma_D[2][i].imag()
        << setw(19)<<  sigma_D[3][i].real() << setw(19)<<  sigma_D[3][i].imag()
        << endl;
    }
    outfile.close();
    outfile.clear();


    cout << ">>> Finish. " << endl;
    
    
    
    
    return 0;
    
}



