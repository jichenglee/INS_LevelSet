/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   variable_pressure.h
 * Author: nsaini3
 *
 * Created on November 7, 2016, 4:12 PM
 */

#ifndef VARIABLE_PRESSURE_H
#define VARIABLE_PRESSURE_H

#include "gauss_siedel.h"

void variable_pressure(vector< vector< vector<double> > > &ustar, vector< vector< vector<double> > > &vstar, vector< vector< vector<double> > > &p, double deltat, vector< vector< vector<double> > > &rho, vector< vector< vector<double> > > &stx, vector< vector< vector<double> > > &sty)
{
    /****Calculate the RHS of matrix****/
    vector< vector<double> > b(xelem, vector<double> (yelem,0.0));
    
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            double hx = area[i][j][1][1];
            double hy = area[i][j][0][0];
            double volume = vol[i][j];
            b[i][j] = (volume/deltat)*((ustar[i][j][0] - ustar[i-1][j][0])/hx + (vstar[i][j][0]-vstar[i][j-1][0])/hy);            
        }
    }
    
    /**Now subtract contribution due to surface tension force***/
    /**First calculate the force at the faces***/
    double stx_face[xelem][yelem];
    double sty_face[xelem][yelem];
    for(int i=0; i<xelem-1; i++)
    {
        for(int j=0; j<yelem-1; j++)
        {
            stx_face[i][j] = 0.5*(stx[i][j][0] + stx[i+1][j][0]);
            sty_face[i][j] = 0.5*(sty[i][j][0] + sty[i][j+1][0]);
        }
    }
    /**Now calculate and add the contribution from force**/
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            double volume = vol[i][j];
            b[i][j] -= (1.0/rho[i][j][0]) * volume*((stx_face[i][j]-stx_face[i-1][j])/area[i][j][1][1] + (sty_face[i][j]-sty_face[i][j-1])/area[i][j][0][0]);
        }
    }
    
    /*Create a matrix to store the rho values at all 4 faces of CV*/
    double elem_rho[xelem][yelem][4];
    for(int i=1; i < xelem-1; i++)
    {
        for(int j=1; j< yelem-1; j++)
        {
            elem_rho[i][j][0] = 0.5*(rho[i+1][j][0] + rho[i][j][0]);
            elem_rho[i][j][1] = 0.5*(rho[i][j+1][0] + rho[i][j][0]);
            elem_rho[i][j][2] = 0.5*(rho[i][j][0] + rho[i-1][j][0]);
            elem_rho[i][j][3] = 0.5*(rho[i][j][0] + rho[i][j-1][0]);
        }
    }
    
    vector< vector< vector<double> > > a(xelem, (vector< vector<double> >(yelem, vector<double>(5,0.0))));
    
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            double hx = area[i][j][1][1];
            double hy = area[i][j][0][0];
            double den = hx*hy*elem_rho[i][j][0]*elem_rho[i][j][1]*elem_rho[i][j][2]*elem_rho[i][j][3];
            a[i][j][0] = hy*hy*elem_rho[i][j][0]*elem_rho[i][j][1]*elem_rho[i][j][3]/den;
            a[i][j][1] = hx*hx*elem_rho[i][j][0]*elem_rho[i][j][1]*elem_rho[i][j][2]/den;
            a[i][j][2] = -(hy*hy*(elem_rho[i][j][1]*elem_rho[i][j][2]*elem_rho[i][j][3] + elem_rho[i][j][0]*elem_rho[i][j][1]*elem_rho[i][j][3]) + hx*hx*(elem_rho[i][j][0]*elem_rho[i][j][2]*elem_rho[i][j][3] + elem_rho[i][j][0]*elem_rho[i][j][1]*elem_rho[i][j][2]) )/den;
            a[i][j][3] = hx*hx*elem_rho[i][j][0]*elem_rho[i][j][2]*elem_rho[i][j][3]/den;
            a[i][j][4] = hy*hy*elem_rho[i][j][1]*elem_rho[i][j][2]*elem_rho[i][j][3]/den;
        }
    }
    
    gs_solver(a,b,p);
}

#endif /* VARIABLE_PRESSURE_H */

