#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>
#include "time.h"

using namespace arma;
using namespace  std;

// The maxoffdiag function, using Armadillo
double maxoffdiag ( mat A, int * k, int * l, int n ){
    double max = 0.0;
    for ( int i = 0; i < n; i++ ) {
        for ( int j = i + 1; j < n; j++ ) {
            if ( fabs(A(i,j)) > max ) {
                max = fabs(A(i,j));
                *l = i;
                *k = j;
            }
        }
    }
    return max;
}
//Jacobi rotations
// Function to find the values of cos and sin
mat rotate ( mat A, mat R, int k, int l, int n )
{
    double s, c;
    if ( A(k,l) != 0.0 ) {
        double t, tau;
        tau = (A(l,l) - A(k,k))/(2*A(k,l));
        if ( tau > 0 ) {
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        } else {
            t = -1.0/( -tau + sqrt(1.0 + tau*tau));
        }
        c = 1/sqrt(1+t*t);
        s = c*t;
    } else {
        c = 1.0;
        s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    // changing the matrix elements with indices k and l
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0; // hard-coding of the zeros
    A(l,k) = 0.0;
    // and then we change the remaining elements
    for ( int i = 0; i < n; i++ ) {
        if ( i != k && i != l ) {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }
        // Finally, we compute the new eigenvectors
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
    return A;
}


//Main program
int main(){
    //Define limits and steps
    double rho_max = 5;
    double rho_min = 0;
    int n = 1000; //number of steps
    double h = (rho_max - rho_min)/n;
    double omega = 5;

    //Printing analytical solutions
//    for(int m=0; m<=3; m++){
//        cout << "Analytical solution= " << 3.*pow((omega/2.),(2./3)) + sqrt(3)*omega*(2*m +1) << endl;
//    }



    //Define arrays
    vec V(n);
    vec rho(n);
    vec d(n);
    double e = -1./(h*h);

    //Fill arrays
    for(int i=0; i<=n; i++){
        rho[i] = rho_min + (i+1)*h;
        //V[i]   = rho[i]*rho[i];
        V[i] = omega*omega*rho[i]*rho[i] + 1/(rho[i]);
        d[i]   = 2./(h*h) + V[i];
    }
    //Fill matrix with values
    mat A = zeros<mat>(n,n);
    A(0,1) = e; //two first elements
    A(0,0) = d[0];//same
    for(int i= 1; i < n-1 ; i++){
        A(i,i-1)   = e;
        A(i,i)     = d[i];
        A(i,i+1)   = e;
    }
    A(n-1,n-2) = e;//Two last elements
    A(n-1,n-1) = d[n-1];//same
    mat B = A;

    //cout << A<< endl; //check that it works

    // Setting up the eigenvector
    mat R = zeros<mat>(n,n);
    for ( int i = 0; i < n; i++ ) {
        for ( int j = 0; j < n; j++ ) {
            if ( i == j ) {
                R(i,j) = 1.0;
            } else {
                R(i,j) = 0.0;
            }
        }
    }

//Start of our mathod
    //Doing the Jacobi rotations with our algorithm
//    int k, l;
//    double epsilon = 1.0e-10;
//    double max_number_iterations = (double) n * (double) n * (double) n;
//    int iterations = 0;
//    double max_offdiag = maxoffdiag ( A, &k, &l, n );

//    clock_t start, finish;
//    start = clock();

//    while ( fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations ) {
//        A = rotate ( A, R, k, l, n );
//        max_offdiag = maxoffdiag ( A, &k, &l, n );
//        iterations++;
//    }

//    finish = clock();
//    cout << "Jacobi time = "<< ((float)(finish-start)/CLOCKS_PER_SEC) <<" sec" << endl;

//    cout << "Number of iterations: " << iterations << "\n";
//    vec Eigenvalues(n);
//    for(int i=0;i<n;i++){
//        Eigenvalues[i] = A(i,i);
//    }
//    Eigenvalues = sort(Eigenvalues,"ascend");
//    for(int i= 0;i<5; i++){
//        cout << Eigenvalues[i]<< endl;
//    }
    //cout << R << endl;
    //Test of iterations dependent on n
    //cout << "3/5* n^2 = " << (float)7./4*n*n << endl;
//End of our method

//Start of Armadillo method
    clock_t starta,finisha;
    starta = clock();
    vec eigval;
    mat eigvec;
    //cout << B << endl; this line is used for testing

    //vec arm = eig_sym(B);
    eig_sym(eigval,eigvec,B);

    finisha = clock();
    cout << "Armadillo time = "<< ((float)(finisha-starta)/CLOCKS_PER_SEC)<<" sec"  << endl;

    for(int i= 0;i<5; i++){
        cout <<eigval[i]<< endl;
    }
    int notortho = 0;
    double test;
    //Check of orthogonality
    for(int i=0; i<=3; i++){
        //cout << eigvec(i,0)*eigvec(i+1,0)<< endl;
        test = eigvec(i,0)*eigvec(i+1,0);
        if(test > 1e-4){
            notortho += 1;
        }
    }
    if(notortho == 0 ){
        cout << "All the eigenvectors are orthogonal!" << endl;
    }
    else{
        cout <<"All the eigenvectors are  NOT orthogonal!" << endl;
    }
    //cout <<eigval<< endl;    
//End of Armadillo method
    //Print output to file
    ofstream out ("output.txt");
    for(int i = 0; i < n; i++){
        out<< rho[i]<< "\t" << eigvec(i,0) << endl;
    }
    out.close();
    cout<< "Output written to file output.txt"<<endl;

    return 0;
}
