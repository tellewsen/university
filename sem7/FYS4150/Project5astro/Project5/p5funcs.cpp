#include "p5funcs.h"

#include <iostream>
#include <lib.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <armadillo>

void random_sphere(cube &data_vector,double R_0,long &idum,int n_planets){
    double u;
    double v;
    double w;
    double c =2*M_PI;
    double r;
    double theta;
    double phi;

    for(int i=0;i<n_planets;i++){
        u = ran0(&idum);
        v = ran0(&idum);
        w = ran0(&idum);
        r= R_0*pow(u,1./3);
        theta= acos(1.-2.*v);
        phi= c*w;
        data_vector(i,0,0) = r*sin(theta)*cos(phi);;
        data_vector(i,1,0) = r*sin(theta)*sin(phi);
        data_vector(i,2,0) = r*cos(theta);
    }
}
void output(cube position_vector, int n_planets,ofstream &ofile,int time_steps){
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for(int j=0;j<time_steps;j++){
        for(int i=0;i<n_planets;i++){
            ofile << setw(15) << setprecision(8) << position_vector(i,0,j);
            ofile << setw(15) << setprecision(8) << position_vector(i,1,j);
            ofile << setw(15) << setprecision(8) << position_vector(i,2,j) <<endl;
        }
    }
}
void outputmass(vec masses,std::ofstream &ofile2){
    ofile2 << setiosflags(ios::showpoint | ios::uppercase);
    int ifinal = masses.size();
    for(int i=0; i< ifinal; i++){
        ofile2 << setw(15) << setprecision(8) << masses[i]<<endl;
    }
}

void accel_calculate(cube data_vector,vec masses,int n_planets,mat &accel,int t,double G,double epsilon){
    double r2=0;
    double force=0;
    vec forces(3);
    forces.fill(0);

    /*One should find a way to capitalize on using newtons third law to minimize
     * computation time for this calculation. The current way wastes a lot of
     * computation time calculating quantities already calculated*/
    for(int i=0;i<n_planets;i++){
        for(int j=0;j<n_planets;j++){
            if(i == j) continue;
            r2 =
                    (data_vector(j,0,t)-data_vector(i,0,t))*
                    (data_vector(j,0,t)-data_vector(i,0,t))

                    + (data_vector(j,1,t)-data_vector(i,1,t))*
                    (data_vector(j,1,t)-data_vector(i,1,t))

                    + (data_vector(j,2,t)-data_vector(i,2,t))*
                    (data_vector(j,2,t)-data_vector(i,2,t))
                    ;

            force = -G*masses(j)/((r2+epsilon*epsilon)*sqrt(r2));

            forces(0) += force*(data_vector(i,0,t)-data_vector(j,0,t));
            forces(1) += force*(data_vector(i,1,t)-data_vector(j,1,t));
            forces(2) += force*(data_vector(i,2,t)-data_vector(j,2,t));

        }
        /*Note here that we dont divide with the mass of the object since
        that quantity is not multiplied with for the forces. This is done
        to save computation time.*/
        accel(i,0) = forces(0);
        accel(i,1) = forces(1);
        accel(i,2) = forces(2);
        force = 0;
        forces.fill(0);
        //cout <<accel << endl;

    }

}
void verlet(cube &data_vector, double step_size, int time_steps, int n_planets, vec masses,double G,double epsilon){
    mat accel(n_planets, 3);
    mat start_data(n_planets,6);

    accel.fill(0);
    start_data.fill(0);
    clock_t start,finish;
    start = clock();
    //Evolve through time
    for(int t=0;t< time_steps-1;t++){
        //Calculate inital acceleration
        accel_calculate(data_vector,masses,n_planets,accel,t,G,epsilon);

        for(int i=0; i<n_planets; i++){
            for(int j=0; j <3; j++){
                start_data(i,j) = data_vector(i,j,t);
                start_data(i,j+3) = data_vector(i,j+3,t);
            }

            for(int j = 0; j<3; j++){
                data_vector(i,j+3,t) = start_data(i,j+3) + 0.5*step_size*accel(i,j);
                data_vector(i,j,t+1) = start_data(i,j) + step_size*start_data(i, j+3) + 0.5*step_size*step_size*accel(i,j);
            }
        }

        accel_calculate(data_vector,masses,n_planets,accel,t,G,epsilon);
        for(int i=0; i<n_planets; i++){
            for(int j=0; j<3; j++){
                data_vector(i,j+3,t) = data_vector(i, j+3,t) + 0.5*accel(i,j)*step_size;
            }
        }
        //Update new position
        for(int i=0; i<n_planets; i++){
            for(int j=0; j<3; j++){
                data_vector(i,j+3,t+1) = data_vector(i,j+3,t);
            }
        }
        //Return start values
        for(int i=0; i<n_planets; i++){
            for(int j=0; j<3; j++){
                data_vector(i,j,t) = start_data(i,j);
                data_vector(i,j+3,t) = start_data(i,j+3);
            }
        }
    }
    finish = clock();
    cout << "Total Velocity-Verlet time = "<< ((float)(finish-start)/CLOCKS_PER_SEC)<<" sec"  << endl;
    cout << "Velocity-Verlet time/step = "<< ((float)(finish-start)/CLOCKS_PER_SEC)/time_steps<<" sec"  << endl;
}
void RK4_solve(cube &data_vector,double step_size,int time_steps, int n_planets,vec masses,double G,double epsilon){

    //Define all the different variables used for the method
    mat k1(n_planets,6);
    mat k2(n_planets,6);
    mat k3(n_planets,6);
    mat k4(n_planets,6);
    mat start_data(n_planets,6);
    mat accel(n_planets,3);

    k1.fill(0);
    k2.fill(0);
    k3.fill(0);
    k4.fill(0);
    start_data.fill(0);
    accel.fill(0);

    //Evolve through time
    clock_t start,finish;
    start = clock();
    for(int t=0;t< time_steps-1;t++){

        //Calculate inital acceleration
        accel_calculate(data_vector,masses,n_planets,accel,t,G,epsilon);

        for(int i=0;i< n_planets;i++){
            //Set start data
            for(int j=0;j<6;j++){
                start_data(i,j) = data_vector(i,j,t);

                //Calculate k1 for each dim
                for(int j=0;j<3;j++){
                    k1(i,j)= step_size*data_vector(i,3+j,t);
                    k1(i,3+j)= step_size*accel(i,j);
                }
            }
            //Save temp vel
            for(int j=0;j<6;j++){
                data_vector(i,j,t)= start_data(i,j) +k1(i,j)/2.;
            }
        }

        //end k1 step

        //Calculate new acceleration at this temporary position
        accel_calculate(data_vector,masses,n_planets,accel,t,G,epsilon);
        for(int i=0;i<n_planets;i++){
            //Calculate k2 for each dim
            k2(i,0)= step_size*data_vector(i,3,t);
            k2(i,1)= step_size*data_vector(i,4,t);
            k2(i,2)= step_size*data_vector(i,5,t);
            k2(i,3)= step_size*accel(i,0);
            k2(i,4)= step_size*accel(i,1);
            k2(i,5)= step_size*accel(i,2);

            //Save temp pos & temp vel
            data_vector(i,0,t)= start_data(i,0) +k2(i,0)/2.;
            data_vector(i,1,t)= start_data(i,1) +k2(i,1)/2.;
            data_vector(i,2,t)= start_data(i,2) +k2(i,2)/2.;
            data_vector(i,3,t)= start_data(i,3) +k2(i,3)/2.;
            data_vector(i,4,t)= start_data(i,4) +k2(i,4)/2.;
            data_vector(i,5,t)= start_data(i,5) +k2(i,5)/2.;

        }

        //Calculate new acceleration again
        accel_calculate(data_vector,masses,n_planets,accel,t,G,epsilon);
        for(int i=0;i<n_planets;i++){
            //Calculate k3 for each dim
            k3(i,0)= step_size*data_vector(i,3,t);
            k3(i,1)= step_size*data_vector(i,4,t);
            k3(i,2)= step_size*data_vector(i,5,t);
            k3(i,3)= step_size*accel(i,0);
            k3(i,4)= step_size*accel(i,1);
            k3(i,5)= step_size*accel(i,2);

            //Save temp vel and pos
            data_vector(i,0,t)= start_data(i,0) +k3(i,0);
            data_vector(i,1,t)= start_data(i,1) +k3(i,1);
            data_vector(i,2,t)= start_data(i,2) +k3(i,2);
            data_vector(i,3,t)= start_data(i,3) +k3(i,3);
            data_vector(i,4,t)= start_data(i,4) +k3(i,4);
            data_vector(i,5,t)= start_data(i,5) +k3(i,5);
        }

        //Calculate new acceleration again
        accel_calculate(data_vector,masses,n_planets,accel,t,G,epsilon);
        for(int i=0;i<n_planets;i++){
            //Calculate k4 for each dim
            k4(i,0)= step_size*data_vector(i,3,t);
            k4(i,1)= step_size*data_vector(i,4,t);
            k4(i,2)= step_size*data_vector(i,5,t);
            k4(i,3)= step_size*accel(i,0);
            k4(i,4)= step_size*accel(i,1);
            k4(i,5)= step_size*accel(i,2);

            //Return start position to data_vector
            data_vector(i,0,t)= start_data(i,0);
            data_vector(i,1,t)= start_data(i,1);
            data_vector(i,2,t)= start_data(i,2);
            data_vector(i,3,t)= start_data(i,3);
            data_vector(i,4,t)= start_data(i,4);
            data_vector(i,5,t)= start_data(i,5);
        }


        for(int i=0;i<n_planets;i++){
            //Update new position
            for(int j=0;j<6;j++){
                data_vector(i,j,t+1) = data_vector(i,j,t) + 1/6.*(k1(i,j) + 2.*k2(i,j) + 2.*k3(i,j) + k4(i,j));    }
        }
    }
    finish = clock();
    cout << "Total RK4 time = "<< ((float)(finish-start)/CLOCKS_PER_SEC)<<" sec"  << endl;
    cout << "RK4 time/step = "<< ((float)(finish-start)/CLOCKS_PER_SEC)/time_steps<<" sec"  << endl;
}

void Euler_solve(cube &data_vector,double step_size,int time_steps, int n_planets,vec masses,double G,double epsilon){

    //Define all the different variables used for the method
    vec k1(6);
    mat accel(n_planets,3);

    k1.fill(0);
    accel.fill(0);

    //Evolve through time
    for(int t=0;t< time_steps-1;t++){

        //Calculate inital acceleration
        accel_calculate(data_vector,masses,n_planets,accel,t,G,epsilon);

        for(int i=0;i< n_planets;i++){

            //Calculate k1 for each dim
            k1(3)= step_size*accel(i,0);
            k1(4)= step_size*accel(i,1);
            k1(5)= step_size*accel(i,2);
            for(int j=3;j<6;j++){
                data_vector(i,j,t+1) = data_vector(i,j,t) + k1(j);
            }

            k1(0)= step_size*data_vector(i,3,t+1);
            k1(1)= step_size*data_vector(i,4,t+1);
            k1(2)= step_size*data_vector(i,5,t+1);

            //Update new position
            for(int j=0;j<3;j++){
                data_vector(i,j,t+1) = data_vector(i,j,t) + k1(j);
            }
        }
    }
}
void Energy(mat &kinetic,mat &potential,cube data_vector, int n_planets,vec masses,int time_steps,double G){
    double r=0;
    double v_sqrd=0;

    for(int t=0;t<time_steps;t++){
        for(int i=0;i<n_planets;i++){
            r = sqrt(data_vector(i,0,t)*data_vector(i,0,t)+
                     data_vector(i,1,t)*data_vector(i,1,t)+
                     data_vector(i,2,t)*data_vector(i,2,t));

            v_sqrd = data_vector(i,3,t)*data_vector(i,3,t)+
                    data_vector(i,4,t)*data_vector(i,4,t)+
                    data_vector(i,5,t)*data_vector(i,5,t);

            kinetic(i,t) = 0.5*masses(i)*v_sqrd;
            for(int j=0;j<n_planets;j++){
                if(i==j) continue;
                r = sqrt(
                        (data_vector(j,0,t)-data_vector(i,0,t))*
                        (data_vector(j,0,t)-data_vector(i,0,t))

                        + (data_vector(j,1,t)-data_vector(i,1,t))*
                        (data_vector(j,1,t)-data_vector(i,1,t))

                        + (data_vector(j,2,t)-data_vector(i,2,t))*
                        (data_vector(j,2,t)-data_vector(i,2,t))
                        );
                potential(i,t) += -masses(i)*masses(j)*G/r;
            }
        }
    }
}

void outputenergy(mat kinetic,mat potential,int n_planets,int time_steps,std::ofstream &ofile4){
    ofile4 << setiosflags(ios::showpoint | ios::uppercase);
    for(int t=0;t<time_steps;t++){
        for(int i=0; i< n_planets; i++){
            ofile4 << setw(15) << setprecision(8) << kinetic(i,t);
            ofile4 << setw(15) << setprecision(8) << potential(i,t) <<endl;
        }
    }
}
