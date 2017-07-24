/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   body_force.h
 * Author: nsaini3
 *
 * Created on November 24, 2016, 4:58 PM
 */

#ifndef BODY_FORCE_H
#define BODY_FORCE_H

void body(elemsclr &sclr, vector< vector< vector<double> > > &st_forcex, vector< vector< vector<double> > > &st_forcey)
{
     /**Compute eps based on grid size*/
    double eps = epsilon*max(xlen/(xelem-2), ylen/(yelem-2));
    
    /**Compute Heavyside function**/
    vector< vector< vector<double> > > H(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));
    
    heavy_func(H, sclr.phi, eps);
    double line = 2.5*2.0*rb_in;
    
    for(int i=0; i<xelem; i++)
    {
        for(int j=0; j<yelem; j++)
        {
            st_forcex[i][j][0] += (min(sclr.phi[i][j][0],0.0)/sclr.phi[i][j][0])*sclr.rho[i][j][0]*vol[i][j]*gx;
            st_forcey[i][j][0] += (min(sclr.phi[i][j][0],0.0)/sclr.phi[i][j][0])*sclr.rho[i][j][0]*vol[i][j]*gy;
            
        }
    }
}

#endif /* BODY_FORCE_H */

