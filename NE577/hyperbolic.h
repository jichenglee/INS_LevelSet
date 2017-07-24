/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   re_distance2.h
 * Author: nsaini3
 *
 * Created on November 8, 2016, 12:17 AM
 */

#ifndef HYPERBOLIC_H
#define HYPERBOLIC_H


void hyperbolic(elemsclr &sclr)
{
    /***Store phi values in a separate matrix***/
    vector< vector<vector<double> > > phi2(xelem, vector<vector<double> > (yelem,vector<double> (zelem,0.0)));
    for(int i=0; i< xelem; i++)
    {
        for(int j=0; j< yelem; j++)
        {
            phi2[i][j][0] = sclr.phi[i][j][0];
        }
    }
    
   
    /**Compute eps based on grid size*/
    double eps = epsilon*max(xlen/(xelem-2), ylen/(yelem-2));
    
     double deltat=re_time;
    
    double ires=0.0;  
    
    bool exitflag = false;
    
    /****Heavyside and delta functions for volume constraint***/
    vector< vector< vector<double> > > H(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));
    vector< vector< vector<double> > > delta(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));
    vector< vector< vector<double> > > grad_phi(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));
    heavy_func(H,sclr.phi,eps);
    delta_func(delta,sclr.phi,eps);
    grad_func(grad_phi, sclr.phi);
    
    
    
    for(int iter=0; iter < re_loops; iter++)
    {
        vector< vector< vector<double> > > lambda(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));
        
        
        vector< vector<vector<double> > > temp_phi2(xelem, vector<vector<double> > (yelem,vector<double> (zelem,0.0)));
        for(int i=1;i<xelem-1;i++)
        {
            for(int j=1;j<yelem-1;j++)
            {
                temp_phi2[i][j][0]=phi2[i][j][0];
            }
        }
        
                
        /*****Now onto calculating fluxes******/
        vector< vector<double> > rhs(xelem, vector<double> (yelem,0.0));
        
        rhs_redist2(rhs, phi2, sclr.phi);
        
        
        
        vector< vector<vector<double> > > phistar(xelem, vector<vector<double> > (yelem,vector<double> (zelem,0.0)));
        
        for(int i=1; i<xelem-1; i++) 
        {
            for(int j=1; j<yelem-1; j++)
            {
                phistar[i][j][0] = phi2[i][j][0] + deltat * (rhs[i][j]);
            }
        }
        
        
        if(vf_control == 1)
        {
            vol_contraint(phistar, sclr.phi, grad_phi, delta, deltat);
        }
        
        //bothscalarBC(phistar);
        level_setBC(phistar);
                
        vector< vector<double> > rhstar(xelem, vector<double> (yelem,0.0));
        //Calculate the star fluxes
        rhs_redist2(rhstar, phistar, sclr.phi);
        
        for(int i=1; i<xelem-1; i++)
        {
            for(int j=1; j<yelem-1; j++)
            {
                phi2[i][j][0] = phistar[i][j][0] + 0.5*deltat *(rhs[i][j] + rhstar[i][j]);
                //cout<<" "<<phi[i][j][0]<<endl;
            }
        }
        
        //level_setBC(phi2);
        if(vf_control == 1)
        {
            vol_contraint(phi2, sclr.phi, grad_phi, delta, deltat);
        }
        level_setBC(phi2);
        
        //bothscalarBC(phistar);
        /****Apply volume constraint****/
        
        if(exitflag == false)
        {
            monitor_res_redist(ires, exitflag, iter,  phi2,  temp_phi2);
        }
        else
        {
            break;
        }
        
        
    }
    
    
    /*Reassign values*/
    for(int i=0; i< xelem; i++)
    {
        for(int j=0; j< yelem; j++)
        {
            sclr.phi[i][j][0] = phi2[i][j][0];
        }
    }
}

#endif /* RE_DISTANCE2_H */

