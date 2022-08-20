#ifndef P5FUNCS_H
#define P5FUNCS_H


#include <armadillo>
using namespace arma;
#endif // P5FUNCS_H


/*Generates random positions inside sphere of radius R_0
for n planets*/
void random_sphere(cube &position_vector,double R_0,long &idum,int n_planets);

// prints to file the results of the calculations for all planets for all times
void output(cube position_vector, int n_planets, std::ofstream &ofile, int time_steps);
/*Prints mass of all planets*/
void outputmass(vec masses,std::ofstream &ofile2);

/*Calculates the gravitational force between all objects in the
data cube and puts that in the forces vector*/
void accel_calculate(cube data_vector, vec masses, int n_planets, mat &accel, int t, double G, double epsilon);
//The RK4 solver
void RK4_solve(cube &data_vector, double step_size, int time_steps, int n_planets, vec masses,double G,double epsilon);
void Euler_solve(cube &data_vector, double step_size, int time_steps, int n_planets, vec masses,double G,double epsilon);
void verlet(cube &data_vector, double step_size, int time_steps, int n_planets, vec masses, double G, double epsilon);
void Energy(mat &kinetic,mat &potential,cube data_vector, int n_planets, vec masses, int time_steps, double G);
void outputenergy(mat kinetic,mat potential,int n_planets,int time_steps,std::ofstream &ofile4);
