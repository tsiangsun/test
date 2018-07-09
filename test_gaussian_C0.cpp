//  test_gaussian_C0.cpp
//
//  " Test time-domain vs freq-domain C0 approximation using Gaussian distribution "
//
//  Compile: g++ -std=c++11 test_gaussian_C0.cpp -o test_gaussian_C0
//
//  Created by Xiang Sun on 9/28/16.
//
// 1/(2 pi) \int du P(u) \int dt e^{i*u*t} = Prob(u=0) = <delta(u)>
// where, P(u) is gaussian distribution with shifted mean value


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <random>
#include <complex>
using namespace std;

const int LEN_TRAJ_EQ = 100000;//total number of data
const int LEN_TRAJ = 5000;
const int LEN = 1024; //512;//number of t choices for FFT
const double DT = 0.05;
const double pi=3.14159265358979324;
const int MAXBIN = 200;
double T0 = -DT*(LEN*0.5);//-DeltaT*LEN/2+DeltaT/2;
double sigma_u = 1; //sigma of normal distribution
double mean_u = 3.7;  //mean value of normal distribution
double u_min = mean_u - 5 * sigma_u;
double u_max = mean_u + 5 * sigma_u;

void FFT(int dir, int m, double *x, double *y); //Fast Fourier Transform, 2^m data
double Integrate(double *data, int n, double dx);
double Integrate_from(double *data, int sp, int n, double dx);
double Sum(double *data, int n);

int main (int argc, char *argv[]) {
    
    int i,j,k;
    int flag(0);
    int id(0); //job number
    int itraj; //index of trajectory
    
    stringstream ss;
    string emptystr("");
    string filename;
    string idstr("");
    string nameapp("");

    double *u = new double [LEN_TRAJ_EQ]; //all input data
    double uav = 0; //average value of data
    
    ofstream outfile;
    outfile.open((emptystr+ "Data_Gauss_input" + nameapp + ".dat").c_str());
    
    
    cout << "------- Prob[u=0] for Gaussian distribution mean=" << mean_u << ", sigma=" << sigma_u << " --------"<< endl;
    cout << "    Total number of input data = " << LEN_TRAJ_EQ << endl;
    
    //1. generate Gaussian distributed data
    
    // using normal distribution random number generator c++11
    /*
     random_device rd;
     mt19937 gen(rd());
     */
    //mt19937 gen;
    //gen.seed(seed);
    
    random_device rd;
    mt19937 gen; // for random generator
    //gen.seed(id*7+153); //the same seed for debug
    gen.seed(rd());
    
    normal_distribution<double> normal_dist(0.0,1.0);
    
    for (i = 0 ; i < LEN_TRAJ_EQ ; i++) {
        //Gaussian sampling
        u[i] = normal_dist(gen) * sigma_u + mean_u;
        uav += u[i];
        outfile << u[i] << endl;
    }
    uav /= LEN_TRAJ_EQ;
    outfile.close();
    outfile.clear();
    
    
    //Method [0]: Analytical prob[u=0]
    double prob_u0;
    prob_u0 = exp(-(0-mean_u)*(0-mean_u)/2/sigma_u/sigma_u) / sqrt(2*pi) / sigma_u;
    
    cout << "[0] --> Analytical Prob[u=0] = " << prob_u0 << endl;
    
    
    
    //2. set up discrete histogram for prob[u]
    double prob[MAXBIN];
    int bin;
    int extreme_count(0);
    int hist_count(0);
    double min_val = u_min;
    double dr = (u_max - u_min)/MAXBIN; //bin size
    for (bin = 0; bin < MAXBIN; bin++) prob[bin] = 0.0;

    // ************** Method [1]: probablity of u=0 *****************
    for (i = 0; i < LEN_TRAJ_EQ; i++) {
        bin = static_cast<int> ((u[i] - min_val)/dr);
        if (bin >=0 && bin < MAXBIN) {
            prob[bin] += 1; //histogramming
            hist_count++;
        }
        else {
            extreme_count++;
            cout << "  1 extreme points. "<< endl;
        }
    }
    for (bin = 0; bin < MAXBIN; bin++) {//normalization of prob
        prob[bin] /= hist_count*dr;
    }
    
    //bin index corresponding to u=0
    bin = static_cast<int> ((0 - min_val)/dr);
    cout << "[1] --> Numerical  Prob[u=0] = " << prob[bin]  << endl;
    
    outfile.open((emptystr + "Prob(u)" + nameapp + ".dat").c_str());
    for (bin = 0; bin < MAXBIN; bin++)  outfile << prob[bin] << endl;
    outfile.close();
    outfile.clear();
    
    
    // ************** Method [2]: time domain direct integral **************
    double *Ct_re = new double [LEN_TRAJ];
    double *Ct_im = new double [LEN_TRAJ];
    double integral_re, integral_im;
    int count(0);
    double k_C0(0);
    double I(0);
    
    //initialize C(t)=exp(i*u*t)
    for (j = 0; j < LEN_TRAJ ; j++) {
        Ct_re[j] = 0;
        Ct_im[j] = 0;
    }
    
    for (i = 0; i < LEN_TRAJ_EQ - LEN_TRAJ; i++) {
        integral_re = 0;
        count = 0;
        for (j = 0; j < LEN_TRAJ ; j++) {
            Ct_re[j] += cos(integral_re);
            Ct_im[j] += sin(integral_re);
            count++;
            integral_re += DT * u[i]; //C-0
        }
    }
    for (j = 0; j < LEN_TRAJ ; j++) {
        Ct_re[j] /= count;
        Ct_im[j] /= count;
    }
    
    outfile.open((emptystr+ "C(t)_re_C0_" + nameapp + ".dat").c_str());
    for (j = 0 ; j < LEN_TRAJ ; j++) {
        outfile << Ct_re[j] << endl;
    }
    outfile.close();
    outfile.clear();
    
    I = 0;
    for (j = 0; j < LEN_TRAJ ; j++) {//integrate 0 to infinity
        I += DT * Ct_re[j];
    }
    
    k_C0 =  I /pi;
    
    cout << "[2] --> 1/2pi*2Re int_0^infty dt <exp(i*u*t)> = " << k_C0  << endl;
    
    
    
    
    
    // ************** Method [3]: time domain via FFT **************
    int mm(0), nn(1); // nn = 2^mm is number of (complex) data to FFT
    double t;
    double shift = T0 / DT; //discrete time origin shift
    while (nn < LEN ) {
        mm++;
        nn *= 2;
    } //nn is the first 2^m that larger LEN
    
    double *corr1 = new double [nn];
    double *corr2 = new double [nn];
    
    double df= 1.0/LEN/DT;
    double domega = df * 2 * pi;
    
    double *corr1_orig = new double [nn]; //shifted origin to T0
    double *corr2_orig = new double [nn];
    
    for (i = 0; i < nn; i++) corr1[i] = corr2[i] = 0; //zero padding
    count = 0;
    
    for (j = 0; j< LEN_TRAJ_EQ; j++) {
        for (k = 0; k < LEN; k++) {
            t = T0 + DT * k;
            integral_re = (u[j]-uav) * t; //C-0 approx. int dt delta u(t)
            corr1[k] += cos(integral_re);
            corr2[k] += sin(integral_re);
        }
        count++;
    }
    
    for (k = 0; k < LEN; k++) {
        corr1[k] /=count;
        corr2[k] /=count;
    }
    
    outfile.open((emptystr + "Ct_re" + nameapp + ".dat").c_str());
    for (i=0; i < nn; i++) outfile << corr1[i]*LEN*DT*0.5/pi << endl;
    outfile.close();
    outfile.clear();
    
    outfile.open((emptystr + "Ct_im" + nameapp + ".dat").c_str());
    for (i=0; i < nn; i++) outfile << corr2[i]*LEN*DT*0.5/pi << endl;
    outfile.close();
    outfile.clear();
    
    FFT(-1, mm, corr1, corr2);//inverse FT
    
    
    for(i = 0; i < nn; i++) { //shift time origin
        corr1_orig[i] = corr1[i] * cos(2*pi*i*shift/nn) - corr2[i] * sin(-2*pi*i*shift/nn);
        corr2_orig[i] = corr2[i] * cos(2*pi*i*shift/nn) + corr1[i] * sin(-2*pi*i*shift/nn);
    }
    
    outfile.open((emptystr + "FFT_Ct_re" + nameapp + ".dat").c_str());
    for (i=nn/2; i < nn; i++) outfile << corr1_orig[i]*LEN*DT*0.5/pi << endl;
    for (i=0; i < nn/2; i++) outfile << corr1_orig[i]*LEN*DT*0.5/pi << endl;
    outfile.close();
    outfile.clear();
    
    outfile.open((emptystr + "FFT_Ct_im" + nameapp + ".dat").c_str());
    for (i=nn/2; i < nn; i++) outfile << corr2_orig[i]*LEN*DT*0.5/pi << endl;
    for (i=0; i < nn/2; i++) outfile << corr2_orig[i]*LEN*DT*0.5/pi << endl;
    outfile.close();
    outfile.clear();
    
    //k = static_cast<int> (uav/domega+0.5);//index for <u>
    k = static_cast<int> (uav/domega);
    double FFTav;
    FFTav = 0.5*(corr1_orig[k]+corr1_orig[k+1]);//average of the left and right FFT result
    
    cout << "[3] --> 1/2pi*DT*LEN * FFT[ <exp(i*du*t)> ](<u>) = " << FFTav*LEN*DT*0.5/pi  << endl;
    cout << "    LEN for FFT = " << LEN << endl;
    cout << "    delta omega = " << domega << endl;
    cout << "    omega_max   = " << LEN/2 * domega << endl;
    cout << "    <u> = " << uav << "   -> index = " << k << endl;
    
    
    return 0;
}





// --------------------------------------------------
// ------------------ SUBROUTINES -------------------
// --------------------------------------------------

double Integrate(double *data, int n, double dx){
    double I =0;
    I += (data[0]+data[n-1])/2;
    for (int i=1; i< n-1; i++) {
        I += data[i];
    }
    I *= dx;
    return I;
}

double Integrate_from(double *data, int sp, int n, double dx){
    double I(0);
    I += (data[sp]+data[n-1])/2;
    for (int i=sp+1; i< n-1; i++) {
        I += data[i];
    }
    I *= dx;
    return I;
}

double Sum(double *data, int n){
    double I = 0;
    for (int i=0; i< n; i++) {
        I += data[i];
    }
    return I;
}


void FFT(int dir, int m, double *x, double *y)
{
  //This code computes an in-place complex-to-complex FFT Written by Paul Bourke
  //x and y are the real and imaginary arrays of N=2^m points.
  // dir =  1 gives forward transform
  // dir = -1 gives reverse transform
  // Formula: forward
  //           N-1
  //           ---
  //        1  \           - i 2 pi k n / N
  //X(n) = ---  >   x(k) e                   = forward transform
  //        1  /                                   n = 0...N-1
  //           ---
  //           k=0
  //
  // Formula: reverse
  //           N-1
  //           ---
  //        1  \           i 2 pi k n / N
  //x(k) = ---  >   X(n) e                  = reverse transform
  //        N  /                                  k = 0..N-1
  //           ---
  //           n=0

    int n,i,i1,j,k,i2,l,l1,l2;
    double c1,c2,tx,ty,t1,t2,u1,u2,z;
    
    // Calculate the number of points
    n = 1;
    for (i=0;i<m;i++)
        n *= 2;
    
    // Do the bit reversal
    i2 = n >> 1; //i2 = (010 ...0)_2,second highest bit of n=(100 ...0)_2
    j = 0; //reversely bit accumulater from the second highest bit, i2.
    for (i=0;i<n-1;i++) {
        if (i < j) {
            tx = x[i]; //swap(i,j)
            ty = y[i];
            x[i] = x[j];
            y[i] = y[j];
            x[j] = tx;
            y[j] = ty;
        }
        //to find the highest non-one bit, k, from the second highest bit
        k = i2;
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k; //add 1 reversly
    }
    
    // Compute the Radix-2 FFT: Cooley-Tukey Algorithm
    c1 = -1.0; // c1+i*c2 = -1 = c^(i 2Pi/2) = W_2, def W_N^j = e^(i 2j*Pi/N)
    c2 = 0.0;
    l2 = 1;
    for (l=0;l<m;l++) {
        l1 = l2;
        l2 <<= 1;
        u1 = 1.0;
        u2 = 0.0;
        for (j=0;j<l1;j++) {
            for (i=j;i<n;i+=l2) {
                //Butterfly calculation of x,y[i] and x,y[i1]:
                //t1+i*t2 =(u1+i*u2)(x[i1]+i*y[i2]) where u1+i*u2=W_N^j=e^(i 2j*Pi/N)
                i1 = i + l1;
                t1 = u1 * x[i1] - u2 * y[i1];
                t2 = u1 * y[i1] + u2 * x[i1];
                x[i1] = x[i] - t1;
                y[i1] = y[i] - t2;
                x[i] += t1;
                y[i] += t2;
            }
            // i1+i*u2 *= c1+i*c2, or W_N
            z =  u1 * c1 - u2 * c2;
            u2 = u1 * c2 + u2 * c1;
            u1 = z;
        }
        //c1+i*c2 = sqrt(c1+i*c2) eg. W_2 --> W_4 ...
        c2 = sqrt((1.0 - c1) / 2.0);
        if (dir == 1)
            c2 = -c2;
        c1 = sqrt((1.0 + c1) / 2.0);
    }
    
    // times STEPS*DeltaT forward FFT (time --> freq)
    /*if (dir == 1) {
     for (i=0; i<n; i++) {
     x[i] *= 1;//DeltaT;
     y[i] *= 1;//DeltaT;
     }
     }*/
    
    // Scaling for inverse transform
    
    if (dir == -1) {
        for (i=0;i<n;i++) {
            x[i] /= n;
            y[i] /= n;
        }
    }
    
    /*
     //for symmetrical FT,
     double sqn;
     sqn = sqrt(n);
     for (i=0;i<n;i++) {
     x[i] /= sqn;
     y[i] /= sqn;
     }
     */
    
    return;
}


