//  Testing 2D FFT code, Complex to Complex
//  src file:  test_2dFFT.cpp
//  Compile: g++ -o test_2dFFT  test_2dFFT.cpp
//  Created by Xiang Sun on 7/6/18.

#include <iostream>
#include <vector>
#include <complex>
using namespace std;

typedef std::complex<double>  Complex;
typedef std::vector<vector<Complex> > Complex_Matrix;
typedef std::vector<vector<double> > Real_Matrix;

const int NX=64;
const int NY=64;
const Complex I(0,1);
const double pi = std::acos(-1);

int FFT2D(Complex_Matrix& c,int nx,int ny,int dir);
int FFT(int dir,int m,double *x,double *y);
bool Powerof2(int n,int *m,int *twopm);

int main(int argc, char *argv[]) {
    int i, j ,k;
    double ax(800), ay(1200), ex;
    Complex a;
    a.real(1.2);
    a.imag(5);
    //cout << a << endl;
    a = Complex(-5, 9);
    //cout << a << endl;
    
    //Time shifting property: if F(t) -> F(t-t0), Maximum at positive t
    //FT[R(t)] = FT[R(t-t0)] * exp(2*pi*i*shiftX/NX + 2*pi*j*shiftY/NX)
    //(i=0,...NX, j=0,...,NY)
    
    double shiftX = NX/2; //time origin shift (X axis),
    double shiftY = NY/2; //time origin shift (Y axis)

    Complex_Matrix Mat(NX,vector<Complex>(NY,0.0));
    for (i=0; i < NX; i++ ) {
        for (j=0; j < NY; j++) {
            ex = exp(-ax*(i-NX/2)/NX*(i-NX/2)/NX - ay*(j-NY/2)/NY*(j-NY/2)/NY);
            Mat[i][j].real(ex);
        }
    }
    FFT2D(Mat, NX, NY, 1);
    
    //fix time shifting
    Complex_Matrix Mat_orig(NX,vector<Complex>(NY,0.0));
    for (i=0; i < NX; i++) {
        for (j=0; j < NY; j++)
            Mat_orig[i][j] = Mat[i][j] * exp(I* Complex(2*pi*i*shiftX/NX) + I*Complex(2*pi*j*shiftY/NY));
    }
    
    for (i=0; i < NX; i++ ) {
        for (j=0; j < NY; j++) {
            cout<< Mat_orig[i][j]<<" ";
        }
        cout << endl;
    }
    
    return 0;
}




/*-------------------------------------------------------------------------
 Perform a 2D FFT inplace given a complex 2D array
 The direction dir, 1 for forward, -1 for reverse
 The size of the array (nx,ny)
 Return false if there are memory problems or
 the dimensions are not powers of 2
 */
int FFT2D(Complex_Matrix& c,int nx,int ny,int dir)
{
    int i,j;
    int m,twopm;
    //double *real,*imag;
    
    /* Transform the rows */
    //real = (double *)malloc(nx * sizeof(double));
    //imag = (double *)malloc(nx * sizeof(double));
    double *real = new double [nx];
    double *imag = new double [nx];
    if (real == NULL || imag == NULL)
        return(-1);
    if (!Powerof2(nx,&m,&twopm) || twopm != nx)
        return(-1);
    for (j=0;j<ny;j++) {
        for (i=0;i<nx;i++) {
            real[i] = c[i][j].real();
            imag[i] = c[i][j].imag();
        }
        FFT(dir,m,real,imag);
        for (i=0;i<nx;i++) {
            c[i][j].real(real[i]);
            c[i][j].imag(imag[i]);
        }
    }
    //free(real);
    //free(imag);
    delete [] real;
    delete [] imag;
    
    /* Transform the columns */
    //real = (double *)malloc(ny * sizeof(double));
    //imag = (double *)malloc(ny * sizeof(double));
    //double *real = new double [ny];
    //double *imag = new double [ny];
    real = new double [ny];
    imag = new double [ny];
    if (real == NULL || imag == NULL)
        return(-1);
    if (!Powerof2(ny,&m,&twopm) || twopm != ny)
        return(-1);
    for (i=0;i<nx;i++) {
        for (j=0;j<ny;j++) {
            real[j] = c[i][j].real();
            imag[j] = c[i][j].imag();
        }
        FFT(dir,m,real,imag);
        for (j=0;j<ny;j++) {
            c[i][j].real(real[j]);
            c[i][j].imag(imag[j]);
        }
    }
    //free(real);
    //free(imag);
    delete [] real;
    delete [] imag;
    
    return(0);
}

/*-------------------------------------------------------------------------
 This computes an in-place complex-to-complex FFT
 x and y are the real and imaginary arrays of 2^m points.
 dir =  1 gives forward transform
 dir = -1 gives reverse transform
 
 Formula: forward
             N-1
             ---
         1   \          - j k 2 pi n / N
 X(n) = ---   >   x(k) e                    = forward transform
         N   /                                n=0..N-1
             ---
             k=0
 
 Formula: reverse
             N-1
             ---
             \          j k 2 pi n / N
 X(n) =       >   x(k) e                    = forward transform
             /                                n=0..N-1
             ---
             k=0
 ******/
int FFT(int dir,int m,double *x,double *y)
{
    long nn,i,i1,j,k,i2,l,l1,l2;
    double c1,c2,tx,ty,t1,t2,u1,u2,z;
    
    /* Calculate the number of points */
    nn = 1;
    for (i=0;i<m;i++)
        nn *= 2;
    
    /* Do the bit reversal */
    i2 = nn >> 1;
    j = 0;
    for (i=0;i<nn-1;i++) {
        if (i < j) {
            tx = x[i];
            ty = y[i];
            x[i] = x[j];
            y[i] = y[j];
            x[j] = tx;
            y[j] = ty;
        }
        k = i2;
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k;
    }
    
    /* Compute the FFT */
    c1 = -1.0;
    c2 = 0.0;
    l2 = 1;
    for (l=0;l<m;l++) {
        l1 = l2;
        l2 <<= 1;
        u1 = 1.0;
        u2 = 0.0;
        for (j=0;j<l1;j++) {
            for (i=j;i<nn;i+=l2) {
                i1 = i + l1;
                t1 = u1 * x[i1] - u2 * y[i1];
                t2 = u1 * y[i1] + u2 * x[i1];
                x[i1] = x[i] - t1;
                y[i1] = y[i] - t2;
                x[i] += t1;
                y[i] += t2;
            }
            z =  u1 * c1 - u2 * c2;
            u2 = u1 * c2 + u2 * c1;
            u1 = z;
        }
        c2 = sqrt((1.0 - c1) / 2.0);
        if (dir == 1)
            c2 = -c2;
        c1 = sqrt((1.0 + c1) / 2.0);
    }
    
    /* Scaling for forward transform */
    if (dir == 1) {
        for (i=0;i<nn;i++) {
            x[i] /= (double)nn;
            y[i] /= (double)nn;
        }
    }
    
    return(0);
}

/*-------------------------------------------------------------------------
 Calculate the closest but lower power of two of a number
 twopm = 2**m <= n
 Return TRUE if 2**m == n
 */
bool Powerof2(int n,int *m,int *twopm)
{
    if (n <= 1) {
        *m = 0;
        *twopm = 1;
        return(false);
    }
    
    *m = 1;
    *twopm = 2;
    do {
        (*m)++;
        (*twopm) *= 2;
    } while (2*(*twopm) <= n);
    
    if (*twopm != n)
        return(false);
    else
        return(true);
}



