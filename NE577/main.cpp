
/*
 * File:   2phase.cpp
 * Author: nsaini3
 *
 * Created on November 1, 2017
 */

#define PI 3.1415926535897
#include <iostream>
# include <stdlib.h>
# include <math.h>
# include <fstream>
#include <cmath>
#include <time.h>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>
#include <fenv.h>
#include <stdio.h>
#include <limits.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <time.h>
#include <cctype>
#include <cstring>
#include <omp.h>


using namespace std;

#include "common.h"



//Include user generated header files. Must come after global variable decalaration
#include "grid.h"
#include "control.h"
#include "output.h"
#include "read_write.h"
#include "bound_cond.h"
#include "pressure_solver.h"
#include "variable_pressure.h"
#include "heavy_delta.h"
#include "initial_conditions.h"
#include "surface_tension.h"
#include "body_force.h"
#include "rhs.h"
#include "functions.h"
#include "rhs_bub.h"
#include "bub_advect.h"
#include "re_distance.h"
#include "hyperbolic.h"
#include "fast_march.h"
#include "direct_redist.h"
#include "calc_vf.h"


int main()
{
    ///Catches mathematical exceptions
    //feenableexcept(FE_INVALID | FE_OVERFLOW |FE_DIVBYZERO);
 
 
    mkdir((getexepath()+"/output").c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir((getexepath()+"/laststep").c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    time_t t1,t2;
    t1 = time(0);
    //********Read control file******/

    omp_set_num_threads(4);

    control();


    xelem=xelem+2; //Include 2 ghost cells
    yelem=yelem+2; //Include 2 ghost cells
    xnode=xelem+1; //Include 2 ghost nodes
    ynode=yelem+1; //Include 2 ghost nodes

    ///*******Resize the vectors and initialize data structures*******//
    x.resize(xnode);
    y.resize(xnode);
    for(int i=0; i < xnode; i++)
    {
        x[i].resize(ynode,0.0);
        y[i].resize(ynode,0.0);
    }

    xc.resize(xelem);
    yc.resize(xelem);
    vol.resize(xelem);
    area.resize(xelem);
    for(int i=0; i<xelem; i++)
    {
        xc[i].resize(yelem,0.0);
        yc[i].resize(yelem,0.0);
        vol[i].resize(yelem,0.0);
        area[i].resize(yelem);
        for(int j=0; j<yelem; j++)
        {
            area[i][j].resize(2);
            for(int k=0; k<2; k++)
            {
                area[i][j][k].resize(2,0.0);
            }
        }
    }

    /*****Read grid and populate element and node vectors******/
    gridread();


    /****Initialize solution vectors*****/
    elemsclr sclr;

    sclr.p.resize(xelem);
    sclr.u.resize(xelem);
    sclr.v.resize(xelem);
    sclr.phi.resize(xelem);
    sclr.rho.resize(xelem);
    sclr.mu.resize(xelem);
    for(int i=0; i<xelem; i++)
    {
        sclr.p[i].resize(yelem);
        sclr.u[i].resize(yelem);
        sclr.v[i].resize(yelem);
        sclr.phi[i].resize(yelem);
        sclr.rho[i].resize(yelem);
        sclr.mu[i].resize(yelem);
        for(int j=0; j<yelem; j++)
        {
            sclr.p[i][j].resize(zelem,0.0);
            sclr.u[i][j].resize(zelem,0.0);
            sclr.v[i][j].resize(zelem,0.0);
            sclr.phi[i][j].resize(zelem,0.0);
            sclr.rho[i][j].resize(zelem,rhof);
            sclr.mu[i][j].resize(zelem,muf);
        }
    }

    initialize(sclr);

    /*Read from file is startstep != 0*/
    read(sclr);

    if(case_tog == 1)
    {
    /*Velocity field for vortex*/
        for(int i=0; i<xelem; i++)
        {
            for(int j=0; j<yelem; j++)
            {
                sclr.u[i][j][0] = -pow(sin(PI*x[i][j]),2.0) * sin(2.0*PI*y[i][j]);
                sclr.v[i][j][0] = pow(sin(PI*y[i][j]),2.0) * sin(2.0*PI*x[i][j]);
            }
        }
    }
    else if(case_tog == 4)
    {
    /*Velocity field for zalesak*/
        for(int i=0; i<xelem; i++)
        {
            for(int j=0; j<yelem; j++)
            {
                sclr.u[i][j][0] = PI*(50.0-y[i][j])/314.0;
                sclr.v[i][j][0] = PI*(x[i][j]-50.0)/314.0;
            }
        }
    }

    //fast_march(sclr);


    vector<double> ires(3,0.0);
    bool exitflag = false;
    int iter;
    int print_count=0;
    double deltat=advect_deltat;
    double init_vf=0.0;

    ofstream out;
    out.open("sim_out.txt",ios::trunc);


    for(iter=startstep; iter<itermax; iter++)
    {

        if(flow_solve == 1)
        {
            vector< vector< vector<double> > > utemp(xelem, vector<vector <double> >(yelem, vector <double> (zelem,0.0)));
            vector< vector< vector<double> > > vtemp(xelem, vector<vector <double> >(yelem, vector <double> (zelem,0.0)));
            for(int i=1;i<xelem-1;i++)
            {
                for(int j=1;j<yelem-1;j++)
                {
                    utemp[i][j][0]=sclr.u[i][j][0];
                    vtemp[i][j][0]=sclr.v[i][j][0];
                }
            }
            vector< vector<double> > rhsx(xelem, vector<double> (yelem,0.0));
            vector< vector<double> > rhsy(xelem, vector<double> (yelem,0.0));
            rhscalc(sclr, rhsx, rhsy, iter, exitflag);

            vector< vector< vector<double> > > ustar(xelem, vector<vector <double> >(yelem, vector <double> (zelem,0.0)));
            vector< vector< vector<double> > > vstar(xelem, vector<vector <double> >(yelem, vector <double> (zelem,0.0)));
            //Predictor Step
            #pragma omp parallel for schedule(dynamic)
            for(int i=1; i<xelem-1; i++)
            {
                for(int j=1; j<yelem-1; j++)
                {
                    ustar[i][j][0] = sclr.u[i][j][0]+deltat*rhsx[i][j];
                    vstar[i][j][0] = sclr.v[i][j][0]+deltat*rhsy[i][j];
                }
            }

            vel_BC(ustar, vstar);

            ///****Calculate contribution from source term - Surface tension force
            /**Note that surface tension force is calculated at the centre of cell at i,j (where p and phi are stored)*/
            vector< vector< vector<double> > > st_forcex(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));
            vector< vector< vector<double> > > st_forcey(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));


            surface(sclr,st_forcex, st_forcey);
            body(sclr,st_forcex,st_forcey);

            if(p_solver == 1)
            {
                pressure(ustar,vstar, sclr.p, deltat);

            }
            else if(p_solver == 2)
            {
                variable_pressure(ustar, vstar, sclr.p, deltat, sclr.rho, st_forcex, st_forcey);
            }

            //Projection Step
            #pragma omp parallel for schedule(dynamic)
            for(int i=1; i<xelem-1; i++)
            {
                for(int j=1; j<yelem-1; j++)
                {

                    sclr.u[i][j][0] = ustar[i][j][0] - deltat*((2.0/(sclr.rho[i][j][0]+sclr.rho[i+1][j][0]))*((sclr.p[i+1][j][0]-sclr.p[i][j][0])/area[i][j][1][1] + 0.5*(st_forcex[i+1][j][0]+st_forcex[i][j][0])));
                    sclr.v[i][j][0] = vstar[i][j][0] - deltat*((2.0/(sclr.rho[i][j][0]+sclr.rho[i][j+1][0]))*((sclr.p[i][j+1][0]-sclr.p[i][j][0])/area[i][j][0][0] + 0.5*(st_forcey[i][j+1][0]+st_forcey[i][j][0])));

                }
            }

            vel_BC(sclr.u, sclr.v);

            cout<<"Step: "<<iter+1;
            out<<"Step: "<<iter+1;
            if(exitflag == false && sol_type == 0)
            {
                monitor_res(ires, exitflag, iter, sclr,utemp,vtemp);
            }
            if(exitflag == true && sol_type == 0)
            {
                cout<<"Flow solution converged"<<endl;
                break;
            }
        }
        if(flow_solve == 0)
        {
            cout<<"Step: "<<iter+1;
            out<<"Step: "<<iter+1;
        }
        /****Bubble Advection and re-distance solvers*****/
        if(advect_solve == 1)
        {
            bub_advect(sclr, iter, deltat);
            //re_distance(sclr);
            if(redist_method == 1)
            {
                hyperbolic(sclr);
            }
            else if(redist_method == 2)
            {
                fast_march(sclr);
                //cout<<"here2"<<endl;
            }

            else if(redist_method == 3)
            {
                direct_redist(sclr);
            }

        }

        print_count++;
        if(print_count == 1)
        {
            output_vtk(sclr,iter);//,st_forcex, st_forcey);
        }
        if(print_count == print_gap)
        {
            print_count = 0;
        }


        /*****Determine time step based on CFL*****/

        if(time_control == 1)
        {
            double cfl;
            timestep_calc(sclr, deltat, cfl);
            cout<<" CFL number: "<<cfl<<" time step: "<<deltat;
            out<<" CFL number: "<<cfl<<" time step: "<<deltat;
        }


        /****Determine Void Fraction*****/
        double vf = 0.0;
        double err = 0.0;
        calc_vf(sclr.phi, init_vf, vf, err);
        cout<<"Void fraction: "<<vf<<"%"<<" Error in vf: "<<err<<"%"<<endl<<endl;
        out<<"Void fraction: "<<vf<<"%"<<" Error in vf: "<<err<<"%"<<endl<<endl;

        cout<<endl;
        out<<endl;
    }



    output_vtk(sclr,iter);
    write(sclr);


    t2 = time(0);
    double seconds = difftime(t2,t1);
    cout<<"Total run time: "<<seconds<<" secs"<<endl;
    out<<"Total run time: "<<seconds<<" secs"<<endl;
    out.close();

}

