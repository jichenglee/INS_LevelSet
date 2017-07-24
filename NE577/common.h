/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   common.h
 * Author: nsaini3
 *
 * Created on October 27, 2016, 1:44 PM
 */

#ifndef COMMON_H
#define COMMON_H

///Global Variable declaration (so that we do not have to pass around information between functions)
double nu;
double cfl;
double tol;
int itermax;
int xelem; //Total elem in x
int yelem; //Total elem in y
int zelem; //Total elem in z
int xnode; //Total nodes in x
int ynode; //Total nodes in y
int znode; //Total nodes in z 
double xlen;
double ylen;
double zlen;
vector< vector<double> > x; //x xoord of nodes
vector< vector<double> > y; //y coord of nodes
    
vector< vector<double> > xc; //x coord of Centroid of cell
vector< vector<double> > yc; //y coord of Centroid of cell
    
vector< vector<double> > vol; //Stores vol of cell
    
vector< vector< vector< vector<double> > > > area; //Stores area of face of parallelogram
    
///Variables for bubble
double rb_in;
double xb_in;
double yb_in;
int advect_steps;
double advect_deltat;
int solnread;
int bub_conv_scheme;
double rhof;
double rhog;
double muf;
double mug;
double epsilon;
double sf_coeff;
double relax;
double ptol;
double re_time;
double re_loops;
int print_gap;
int startstep;
double gx;
double gy;

string inttostr (int n)
{
    stringstream ss;
    ss<<n;
    return ss.str();
}

std::string getexepath()
{
  char result[ MAX_PATH ];
  return std::string( result, GetModuleFileName( NULL, result, MAX_PATH ) );
}


/*****Some Simulation Control variables****/
int sf_toggle;
int flow_solve;
int p_solver;
int x_bound;
int y_bound;
int advect_solve;
int sol_type;
int vf_control;
int time_control;
double max_cfl;
int redist_method;
int case_tog;




struct elemsclr
{
    vector< vector< vector<double> > > p;
    vector< vector< vector<double> > > u;
    vector< vector< vector<double> > > v;
    vector< vector< vector<double> > > phi;
    vector< vector< vector<double> > > rho;
    vector< vector< vector<double> > > mu;
};
#endif /* COMMON_H */

