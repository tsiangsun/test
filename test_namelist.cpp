#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <set>

using namespace std;

int main () {
    int i,j,k;
    
    //variables to initialize
    double eps(0);
    double Gamma_DA(1);
    double beta(1);
    double zeta(1);
    double DT(0.01);
    int Ntraj(100);
    int LEN_TRAJ(1);
    
    
    map<string, double> varmap;
    //set<string> exclude = {",", "="};
    ifstream infile;

    infile.open("test_namelist.inp");
    if (!infile.is_open()) {
        cout << "Error: input file cannot open"<< endl;
        return -1;
    }
    
    i=0;
    string varstr;
    string equ;
    string commentline;
    while (infile >> varstr) {
        //if (varstr[0] == '#') {
        if (varstr.at(0) == '#') {
        	getline(infile, commentline);
        	cout << "Comment: " << varstr << commentline << endl;
        	continue;
        }
        i++;
        infile >> equ;
        if (equ != "=" ) cout << "missing = sign" << endl;
        infile >> varmap[varstr];
        cout << varstr << ": " << varmap[varstr] << endl;
    }
    
    cout << ">> number of entries = " << varmap.size() << endl;
    
    if (varmap.find("eps") != varmap.end()) eps = varmap["eps"];
    if (varmap.find("Gamma_DA") != varmap.end()) Gamma_DA = varmap["Gamma_DA"];
    if (varmap.find("beta") != varmap.end()) beta = varmap["beta"];
    if (varmap.find("zeta") != varmap.end()) zeta = varmap["zeta"];
    if (varmap.find("DT") != varmap.end()) DT = varmap["DT"];
    if (varmap.find("Ntraj") != varmap.end()) Ntraj = varmap["Ntraj"];
    if (varmap.find("LEN_TRAJ") != varmap.end()) LEN_TRAJ = varmap["LEN_TRAJ"];
    
    cout << " eps = " << eps << endl;
    cout << " Gamma_DA = " << Gamma_DA << endl;
    cout << " beta = " << beta << endl;
    cout << " zeta = " << zeta << endl;
    cout << " DT = " << DT << endl;
    cout << " Ntraj = " << Ntraj << endl;
    cout << " LEN_TRAJ = " << LEN_TRAJ << endl;
    
    return 0;
}

