#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <mpi.h>
#include <lib.h>
#include <p5funcs.h>
#include <armadillo>
#include <random>
#include <time.h>

using namespace std;
using namespace arma;

ofstream ofile;
ofstream ofile2;
ofstream ofile3;
ofstream ofile4;

int main(int argc, char* argv[])
{
    //Numerical constants
    double M_sun = 1.989*1e30;//kg
    double AU  = 1.496e11;
    double ly  = 9.4607e15;
    //outputfilenames
    char *outfilename,*outfilename2,*outfilename3,*outfilename4;
    //RNG seed
    long idum = time(NULL);

    // MPI initializations
    int my_rank, numprocs;
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    // Read in output file, abort if there are too few command-line arguments
    if (my_rank == 0 && argc <= 2) {
        cout << "Bad Usage: " << argv[0] <<
                " read output file" << endl;
        exit(1);
    }
    if (my_rank == 0 && argc > 2) {
        outfilename=argv[1];
        outfilename2=argv[2];
        outfilename3=argv[3];
        outfilename4=argv[4];

        ofile.open(outfilename);
        ofile2.open(outfilename2);
        ofile3.open(outfilename3);
        ofile4.open(outfilename4);
    }

    //Decide which method to use.
    int method = 1;


    //Initialize number of planets variable
    int n_planets;

    //Create mass vector contatining all planets


    /*Set timestep depending on method and create the
     * data vector corresponding to this time step*/
    double final_time;
    double step_size;
    int time_steps;
    double epsilon;
    if(method==0){
        n_planets = 2;
        final_time = 3600.*24.*365.*20;
        step_size = 3600.*24.*2;
        time_steps = final_time/step_size;
    }
    if(method==1){
        n_planets = 100;
        final_time = 5;
        step_size = 0.001;
        time_steps = final_time/step_size;
        epsilon = 0.1;
    }
    vec masses(n_planets);
    masses.fill(0);

    cube data_vector(n_planets,6,time_steps);
    data_vector.fill(0);

    if(method==0){
        /*Sun-earth system with sun at 0 and earth at 1AU on x axis*/
        data_vector(1,4,0) = 29.78e3;//m/s
        data_vector(1,0,0) = 1*AU;//m
        masses(0) = 1.9898e30;//kg
        masses(1) = 5.972e24; //kg
    }
    double R_0;
    if(method==1){
        //Initialize positions of all planets at random inside sphere of radius R_0
        R_0 = 20; //Radius of sphere we want objects inside
        random_sphere(data_vector,R_0,idum,n_planets);

        /*Set masses of all objects with gaussian distribution
        around 10solar masses, with standard deviation 1solar mass*/
        default_random_engine generator;
        normal_distribution<double> distribution(10,1);
        for(int i=0;i<n_planets;i++){
            masses[i] = distribution(generator)/5.;
        }
    }

    cout <<"time steps: "<< time_steps << endl;


    /*Define G depending on if one is using dimensionless variables
     * or not*/
    double G;
    if(method==0){
        /*Method 0 is the regular G when using SI units*/
        G = 6.67e-11;//m3kg-1s-2
    }
    if(method==1){
        G = M_PI*M_PI*pow(R_0,3)/(8*sum(masses));
    }



    //Initialize velocity of all planets
    /*Random velocities in normal distributon*/
    //    default_random_engine generator1;
    //    normal_distribution<double> distribution1(1e5,1e2);
    //    for(int i=0;i<n_planets;i++){
    //        for(int j =0;j<3;j++){
    //            data_vector(i,3+j,0) = distribution1(generator1);
    //        }
    //    }


    //set masses manually for two body problem sun earth

    //Use the Runge Kutta 4 method to evolve system through time
    //RK4_solve(data_vector,step_size,time_steps,n_planets,masses,G,epsilon);
    //Print positions of all planets at all times to file
    //output(data_vector,n_planets,ofile,time_steps);

    /*  Solve system using Velocity-Verlet method*/
    verlet(data_vector, step_size, time_steps,n_planets,masses,G,epsilon);
    //Print positions of all planets at all times to file
    output(data_vector,n_planets,ofile2,time_steps);

    //Print masses of all objects in simulation
    outputmass(masses,ofile3);

    /*Solve system using Euler's method*/
    //Euler_solve(data_vector,step_size,time_steps,n_planets,masses,G,epsilon);



    /*Calculate and write kinetic,potential, and total energy to file*/
    mat kinetic(n_planets,time_steps);
    kinetic.fill(0);
    mat potential(n_planets,time_steps);
    potential.fill(0);
    Energy(kinetic,potential,data_vector, n_planets, masses, time_steps, G);

    mat energies(n_planets,time_steps);
    energies.fill(0);
    energies = kinetic+potential;
    outputenergy(kinetic,potential,n_planets,time_steps,ofile4);

    return 0;
}

