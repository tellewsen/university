#include <iostream>
#include <armadillo>
#include "time.h"

using namespace arma;
using namespace std;



int main()
{
    //Define arrays
    int n = 1000;
    vec a(n+1);
    a.fill(-1);

    vec b(n+2);
    b.fill(2);

    vec c(n+1);
    c.fill(-1);

    vec b_tilde(n+2);

    double h = 1./(n+1);
    //cout << h<< endl;

    //Define and fill x vector
    vec x(n+2);
    for(int i=0; i <= n+1 ; i++){
        x[i] = i*h;
    }

    //Define and fill u vector with exact solution
    vec u(n+2);
    for(int i = 0; i <= n+1 ; i++){
        u[i] = 1 -(1-exp(-10))*x[i]-exp(-10*x[i]);
    }

    //Fill b_tilde with the source function at x values in x vector
    for(int i = 0; i <= n+1 ; i++){
        b_tilde[i] = h*h*100*exp(-10*x[i]);
    }
    vec b_tilde_orig = b_tilde;


    //Initialize tridiagonal timer
    clock_t start, finish;
    start = clock();


    //Forward substituion
    for(int i=0 ; i <= n-1 ; i++){
        //a[i]   = 0;
        b[i+1] = b[i+1]*b[i] + c[i];
        c[i+1] = c[i+1]*b[i];
        b_tilde[i+1] = b_tilde[i+1]*b[i] +b_tilde[i];
    }
    //a[n] = 0;
    b[n+1] = b[n+1]*b[n] + c[n];

    //Backward substition
    for(int i=n ; i >=0 ; i--){
        b_tilde[i] = b_tilde[i] - b_tilde[i+1]*c[i]/b[i+1];
        //c[i] = 0;
    }
    //Find v
    vec v(n+2);
    for(int i=0; i<=n+1; i++){
        v[i] = b_tilde[i]/b[i];
    }


    //Terminate tridiagonal timer
    finish = clock();
    cout << "Tridiagonal time = " << ( (float)(finish - start)/CLOCKS_PER_SEC)<< endl;



    //Compute error in v
    vec epsilon(n+2);
    for(int i = 1 ; i<=n+1;i++){
        epsilon[i] = log10(abs((v[i] -u[i])/u[i]));
    }
    //cout <<epsilon<<endl;


    //Compute max error with current n
    float maxerror;

    for(int i= 0 ; i<n+1;i++){
        if(maxerror < epsilon[i]){
            maxerror = epsilon[i];
        }
    }
    cout << "Maxerror = "<< maxerror << endl;


    //Solve using Gaussian elimination in armadillo
    mat A = zeros<mat>(n+2,n+2);
    A.diag() +=2;
    A.diag(-1) -= 1;
    A.diag(1)  -= 1;
    //cout << A<< endl;
    vec b_tilde_gauss = b_tilde_orig;
    vec v_gauss;


    v_gauss = solve(A,b_tilde_gauss);
    //cout << v_gauss << endl;


    //Solve using LU decompositin in armadillo
    vec b_tilde_LU = b_tilde_orig;
    vec v_LU;
    mat L, U, P;

    //Initialize LU decomposition timer
    clock_t start1, finish1;
    start1 = clock();

    lu(L,U, P, A);

    //Terminate LU decomposition timer
    finish1 = clock();
    cout <<"LU time = " << ( (float)(finish1 - start1)/CLOCKS_PER_SEC)<< endl;

    v_LU = inv(U)*inv(L)*b_tilde_LU;


    //Compute error in v_LU
    vec epsilon_LU(n+2);
    for(int i = 1 ; i<=n+1;i++){
        epsilon_LU[i] = log10(abs((v_LU[i] -u[i])/u[i]));
    }

    //Compute max error with current n
    float maxerror_LU;

    for(int i= 0 ; i<n+1;i++){
        if(maxerror_LU < epsilon_LU[i]){
            maxerror_LU = epsilon_LU[i];
        }
    }
    cout << "Maxerror LU = "<< maxerror_LU << endl;

    //cout << L << endl;
    //cout << U << endl;
    //cout << P << endl;
    //cout <<v_LU<< endl;


    //Printout FLOPS
    cout << "FLOPS_tri = " << 8.*n << endl;
    cout << "FLOPS_gauss = " << 3./2*n*n*n + n*n<< endl;
    cout << "FLOPS_LU = " << 1.*n*n*n << endl;




    //Print stuff for testing
    //cout << "a = "<< endl<< a << endl;
    //cout << "b = " << endl << b << endl;
    //cout << "c = "<< endl << c << endl;
    //cout << "b_tilde = " << endl << b_tilde << endl;
    //cout << "u = " << endl << u << endl;
    //cout << "v = " << endl << v << endl;
    //cout << "x = " << endl << x << endl;

    //Print output to file
    ofstream out ("output.txt");
    for(int i = 0; i <=n+1; i++){
        out<< x[i]<< "\t" << u[i] <<"\t" <<v[i] <<"\t" << v_gauss[i]<< "\t" << v_LU[i] << endl;
    }
    out.close();
    cout<< "Output written to file output.txt"<<endl;
    return 0;

}
