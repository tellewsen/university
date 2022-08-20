//   This is a simple program which tests the trapezoidal rule, Simpsons' rule,
//   and Gaussian quadrature using Legendre and Laguerre polynomials
//   It integrates the simple function x* exp(-x) for the interval
//   x \in [0,infty). The exact result is 1. For Legendre based quadrature a
//   tangent mapping is also used.

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <lib.h>
#include <armadillo>

#define EPS 3.0e-14
#define MAXIT 10
#define   ZERO       1.0E-10
using namespace std;
using namespace arma;
//     Here we define various functions called by the main program

//Our function
double int_function(double x1, double y1, double z1, double x2, double y2, double z2);
//Our function in spherical coordinates
double int_function2(double x1, double x2, double costheta1, double costheta2, double phi1, double phi2);
//Gaussian quadrature
void gauss_laguerre(double *, double *, int, double);
void gauleg(double, double, double *, double *, int);
double gammln(double);
//Our Function for brute force Monte Carlo integration
double brute_force_MC(double *);
//Our Function in spherical coords for imp samp Monte Carlo integration
double imp_samp_MC(double *x);


//   Main function begins here
int main()
{
    int N; //Iterations
    cout << "Read in the number of integration points" << endl;
    cin >> N;

    //Start of Gauss-Legendre method

    double a, b;//Integration limits
    cout << "Read in integration limits" << endl;
    cin >> a >> b;

    clock_t start, finish;
    start = clock();

    double *x = new double [N];
    double *w = new double [N];
    //Set up the mesh points and weights
    gauleg(a,b,x,w,N);

    //Initialize the sum
    double int_gauss = 0.;

    //Six-double loops
    //Note here that you get n^6 calculations total here, so this adds up VERY quickly to enormous numbers

    for (int i=0;i<N;i++){
        for (int j = 0;j<N;j++){
            for (int k = 0;k<N;k++){
                for (int l = 0;l<N;l++){
                    for (int m = 0;m<N;m++){
                        for (int n = 0;n<N;n++){
                            int_gauss+=w[i]*w[j]*w[k]*w[l]*w[m]*w[n]*int_function(x[i],x[j],x[k],x[l],x[m],x[n]);
                        }
                    }
                }
            }
        }
    }
    finish = clock();
    cout << "Gauss-Legendre time = "<< ((float)(finish-start)/CLOCKS_PER_SEC) <<" sec" << endl;
    //End of Gauss-Legendre

    //Start of Gauss-Laguerre intergration
    clock_t start1, finish1;
    start1 = clock();

    double alf = 2.;

    //Initialize vectors
    double *r = new double[N];
    double *theta = new double[N];
    double *phi   = new double[N];
    //Init weights
    double *wr= new double[N];
    double *wphi  = new double[N];
    double *wtheta  = new double[N];

    //Calculate mesh points and weights
    gauss_laguerre(r,wr,N,alf);
    gauleg(0., M_PI,theta,wtheta,N);
    gauleg(0.,2.*M_PI,phi,wphi,N);

    //Initialize sum
    double int_gausslag = 0.;

    //Six-double loop to calculate integral
    for (int i=1;i<=N;i++){
        for (int j = 1;j<=N;j++){
            for (int k = 0;k<N;k++){
                for (int l = 0;l<N;l++){
                    for (int m = 0;m<N;m++){
                        for (int n = 0;n<N;n++){
                            int_gausslag+=wr[i]*wr[j]*wtheta[k]*wtheta[l]*wphi[m]*wphi[n]*
                                    int_function2(r[i],r[j],theta[k],theta[l],phi[m],phi[n]);                            }
                    }
                }
            }
        }
    }
    double alpha = 2.;
    int_gausslag /= pow(2.*alpha,5);

    finish1 = clock();
    cout << "Gauss-Laguerre time = "<< ((float)(finish1-start1)/CLOCKS_PER_SEC) <<" sec" << endl;
    //End of gauss-laguerre computation




    //    //Monte Carlo without importance sampling
    long idum=-1 ;

    cout << "Read in the number of Monte-Carlo samples" << endl;
    cin >> N;



    clock_t start2, finish2;
    start2 = clock();

    double xmc[6], fx;
    double int_mc = 0.;
    double variance = 0.;
    double sum_sigma= 0. ;
    double length = 3.; // we fix the max size of the box to L=3
    double jacobidet = pow((2*length),6);
    // evaluate the integral without importance sampling
    for ( int i = 1; i <= N; i++){
        // x[] contains the random numbers for all dimensions
        for (int j = 0; j< 6; j++) {
            xmc[j]=-length+2*length*ran0(&idum);
        }
        fx=brute_force_MC(xmc);
        int_mc += fx;
        sum_sigma += fx*fx;
    }
    int_mc = int_mc/((double) N );
    sum_sigma = sum_sigma/((double) N );
    variance=sum_sigma-int_mc*int_mc;

    finish2 = clock();
    cout << "MC time = "<< ((float)(finish2-start2)/CLOCKS_PER_SEC) <<" sec" << endl;


    //Monte Carlo with importance sampling
    clock_t start3, finish3;
    start3 = clock();

    double xis[6];
    double int_mcis = 0.;
    double varianceis = 0.;
    double sum_sigmais= 0. ;
    double jacobidetis = 4.*pow(M_PI,4.)/16.;
    double fxis = 0;

    for ( int i = 0; i <= N; i++){
        // x[] contains the random numbers for all dimensions
        for (int j = 0; j<= 1; j++) {
            xis[j]= -log(1-ran0(&idum))/4.;//Exponential distribution
            xis[j+2] = M_PI*ran0(&idum);//Uniform [0,2pi]
            xis[j+4] = 2.*M_PI*ran0(&idum);//Uniform [0,pi]
        }

        fxis=imp_samp_MC(xis);
        int_mcis += fxis;
        sum_sigmais += fxis*fxis;
    }

    int_mcis    = int_mcis/((double) N );
    sum_sigmais = sum_sigmais/((double) N );
    varianceis  = sum_sigmais- int_mcis*int_mcis;
    finish3 = clock();
    cout << "MC imp samp time = "<< ((float)(finish3-start3)/CLOCKS_PER_SEC) <<" sec" << endl;
    //End of MC with importance sampling

    //Final output
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << "Gaussian-Legendre quad = "<< setw(20) << setprecision(15)  << int_gauss << endl;
    cout << "Gaussian-Laguerre quad = " << setw(20) << setprecision(15) << int_gausslag << endl;
    cout << "Monte carlo result=      " << setw(10) << setprecision(8) << jacobidet*int_mc<< endl;
    cout << "Sigma=                   " << setw(10) << setprecision(8) << jacobidet*sqrt(variance/((double) N)) << endl;
    cout << "Mc imp samp result=      " << setw(10) << setprecision(8) << jacobidetis*int_mcis<< endl;
    cout << "Sigma imp samp=          " << setw(10) << setprecision(8) << jacobidetis*sqrt(varianceis/((double) N)) << endl;
    cout << "Exact solution =         " << 5.*M_PI*M_PI/16./16. << endl;

    //Plotting e^-x
    vec t(101);
    vec function(101);
    for(int i=0; i < 101 ; i++){
        t[i] = i/10.;
        function[i] = exp(-2.*t[i]);
    }

    //Print output to file
    ofstream out ("output.txt");
    for(int i = 0; i <=100; i++){
        out << t[i] << "\t" << function[i] << endl;
    }
    out.close();
    cout<< "Output written to file output.txt" << endl;
    return 0;
}  // end of main program

// this function defines the integrand to integrate
double brute_force_MC(double *x)
{
    double alpha = 2.;
    // evaluate the different terms of the exponential
    double exp1=-2*alpha*sqrt(x[0]*x[0] + x[2]*x[2] + x[4]*x[4]);
    double exp2=-2*alpha*sqrt(x[1]*x[1] + x[3]*x[3] + x[5]*x[5]);
    double deno=sqrt(pow((x[0]-x[1]),2)+pow((x[2]-x[3]),2)+pow((x[4]-x[5]),2));

    if(deno <pow(10.,-6.)) {
        return 0;}
    else
        return exp(exp1+exp2)/deno;


} // end function for the integrand
double imp_samp_MC(double *r)
{
    // evaluate the different terms of the exponential
    double cosbeta = cos(r[2]) * cos(r[3]) + sin(r[2]) * sin(r[3]) * cos(r[4]-r[5]);
    double deno=sqrt(r[0]*r[0] + r[1]*r[1] - 2.*r[0]*r[1]*cosbeta);

    //Ignore terms with complex deno or very small deno
    if(deno <pow(10.,-6.) || isnan(deno) ){
        return 0;}
    else
        return r[0]*r[0]*r[1]*r[1]/deno* sin(r[2]) * sin(r[3]);
}

//  this function defines the function to integrate
double int_function(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double alpha = 2.;
    // evaluate the different terms of the exponential
    double exp1=-2*alpha*sqrt(x1*x1+y1*y1+z1*z1);
    double exp2=-2*alpha*sqrt(x2*x2+y2*y2+z2*z2);
    double deno=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
    if(deno <pow(10.,-6.)) {
        return 0;}
    else
        return exp(exp1+exp2)/deno;
} // end of function to evaluate

double int_function2(double r1, double r2, double theta1, double theta2, double phi1, double phi2)
{
    double cosbeta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
    // evaluate the different terms of the exponential
    double deno=sqrt(r1*r1 + r2*r2 -2.*r1*r2*cosbeta);
    if(deno <pow(10.,-6.) || isnan(deno) ) {
        return 0;}
    else
        return 1./deno * sin(theta1)*sin(theta2);
} // end of function to evaluate

void gauss_laguerre(double *x, double *w, int n, double alf)
{
    int i,its,j;
    double ai;
    double p1,p2,p3,pp,z,z1;

    for (i=1;i<=n;i++) {
        if (i == 1) {
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        } else if (i == 2) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        } else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                  (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=1;its<=MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
        x[i]=z;
        w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    }
}
// end function gauss_laguerre

double gammln( double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
                          24.01409824083091,-1.231739572450155,
                          0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}
// end function gammln
