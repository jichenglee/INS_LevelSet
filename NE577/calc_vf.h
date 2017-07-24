/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   calc_vf.h
 * Author: nsaini3
 *
 * Created on November 23, 2016, 1:11 PM
 */

#ifndef CALC_VF_H
#define CALC_VF_H

void calc_vf(vector< vector< vector< double > > > &phi, double &init_vf, double &vf, double &err)
{
    double eps=epsilon*max(xlen/(xelem-2), ylen/(yelem-2));
    vector< vector< vector<double> > > H(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));
    
    heavy_func(H,  phi, eps);
    
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            vf += (1.0 - H[i][j][0])*area[i][j][0][0]*area[i][j][1][1];
        }
    }
    
    double totvol = xlen*ylen;
    vf = vf*100.0/totvol;
    
    if(init_vf == 0.0)
    {
        init_vf = vf;
    }
    
    err = (vf - init_vf)*100.0/init_vf;
    
    
}

#endif /* CALC_VF_H */

