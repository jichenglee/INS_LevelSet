/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   re_distance.h
 * Author: nsaini3
 *
 * Created on November 7, 2016, 8:24 PM
 */

#ifndef RE_DISTANCE_H
#define RE_DISTANCE_H

#include "rhs_bub_redist.h"
void monitor_res_redist(double &ires, bool &exitflag, int iter, vector< vector<vector<double> > > phi, vector< vector<vector<double> > > phitemp)
{
    double res=0.0;
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            res=res + pow(phi[i][j][0]-phitemp[i][j][0],2.0)*vol[i][j];
        }
    }
    
    res=sqrt(res);
    
    
    
    if(iter == 0)
    {
            ires=res;
    }
    else
    {
        //cout<<"Step: "<<iter<<" phi residual: "<<res/ires<<endl;//<<" V vel residual: "<<res[1]/ires[1]<<endl;
        
        if(res/ires < tol)
        {
            exitflag=true;
        }
    }
    
}



void re_distance(elemsclr &sclr)
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
    
    for(int iter=0; iter < re_loops; iter++)
    {
        /**Compute Heavyside function**/
        vector< vector< vector<double> > > H(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));
        //heavy(H,phi2,eps);
        
        vector< vector< vector<double> > > signnew(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));
        for(int i=0; i<xelem; i++)
        {
            for(int j=0; j<yelem; j++)
            {
                signnew[i][j][0] = 2.0*(H[i][j][0] -0.5);
            }
        }
        
        vector< vector<vector<double> > > temp_phi2(xelem, vector<vector<double> > (yelem,vector<double> (zelem,0.0)));
        for(int i=1;i<xelem-1;i++)
        {
            for(int j=1;j<yelem-1;j++)
            {
                temp_phi2[i][j][0]=phi2[i][j][0];
            }
        }
        
        //cout<<iter<<endl;
        /*Calculate convection velocites*/
        vector< vector< vector<double> > > grad_phix(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0))); 
        vector< vector< vector<double> > > grad_phiy(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));
        
        vector< vector< vector<double> > > phiRface(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));
        vector< vector< vector<double> > > phiTface(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));
        for(int i=0; i<xelem-1; i++)
        {
            for(int j=0; j<yelem-1; j++)
            {
                phiRface[i][j][0] = 0.5*(phi2[i+1][j][0] + phi2[i][j][0]);
                phiTface[i][j][0] = 0.5*(phi2[i][j+1][0] + phi2[i][j][0]);
            }
        }

        for(int j=1; j<yelem-1; j++)
        {
            for(int i=1; i<xelem-1; i++)
            {
                grad_phix[i][j][0] = (phiRface[i][j][0] - phiRface[i-1][j][0])/area[i][j][1][1];
                grad_phiy[i][j][0] = (phiRface[i][j][0] - phiRface[i][j-1][0])/area[i][j][0][0];
            }
        }
        periodicBC(grad_phix);
        periodicBC(grad_phiy);
        gradBC(grad_phix);
        gradBC(grad_phiy);
        
        vector< vector<vector<double> > > ucen(xelem, vector<vector<double> > (yelem,vector<double> (zelem,0.0))); //Stored at cell centers
        vector< vector<vector<double> > > vcen(xelem, vector<vector<double> > (yelem,vector<double> (zelem,0.0)));
        
        for(int j=1; j<yelem-1; j++)
        {
            for(int i=1; i<xelem-1; i++)
            {
                double mag_phi = sqrt(pow(grad_phix[i][j][0],2.0) + pow(grad_phiy[i][j][0],2.0));
                ucen[i][j][0] = signnew[i][j][0]*grad_phix[i][j][0]/mag_phi;
                vcen[i][j][0] = signnew[i][j][0]*grad_phiy[i][j][0]/mag_phi;
                //cout<<grad_phiy[i][j][0]<<" ";
            }
            //cout<<endl;
        }
        //exit(0);
        
        periodicBC(ucen);
        periodicBC(vcen);
        zerogradBC(ucen); //Note velocity at the wall should not be zero
        zerogradBC(vcen); //Note velocity at the wall should not be zero
        
        
        
        
        /*****Now onto calculating fluxes******/
        vector< vector<double> > rhsx(xelem, vector<double> (yelem,0.0));
        vector< vector<double> > rhsy(xelem, vector<double> (yelem,0.0));
        
        rhs_bub(rhsx, rhsy, ucen, vcen, phi2);
        
        vector< vector<vector<double> > > phistar(xelem, vector<vector<double> > (yelem,vector<double> (zelem,0.0)));
        for(int i=1; i<xelem-1; i++) 
        {
            for(int j=1; j<yelem-1; j++)
            {
                phistar[i][j][0] = phi2[i][j][0] + deltat * (signnew[i][j][0] + rhsx[i][j] + rhsy[i][j]);
            }
        }
        periodicBC(phistar);
        zerogradBC(phistar);
        //bothscalarBC(phistar);
                
        vector< vector<double> > rhstarx(xelem, vector<double> (yelem,0.0));
        vector< vector<double> > rhstary(xelem, vector<double> (yelem,0.0));
        //Calculate the star fluxes
        rhs_bub(rhstarx, rhstary, ucen, vcen, phistar);
        
        vector< vector< vector<double> > > Hstar(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));
//        /heavy(Hstar,phistar,eps);
        
        vector< vector< vector<double> > > signnew2(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));
        for(int i=0; i<xelem; i++)
        {
            for(int j=0; j<yelem; j++)
            {
                signnew2[i][j][0] = 2.0*(Hstar[i][j][0] -0.5);
            }
        }
        
        
        for(int i=1; i<xelem-1; i++)
        {
            for(int j=1; j<yelem-1; j++)
            {
                phi2[i][j][0] = phistar[i][j][0] + 0.5*deltat *(signnew[i][j][0] + signnew2[i][j][0] + rhsx[i][j] + rhstarx[i][j] + rhsy[i][j] + rhstary[i][j]);
                //cout<<" "<<phi[i][j][0]<<endl;
            }
        }
        
        periodicBC(phi2);
        zerogradBC(phi2);
        
        if(exitflag == false)
        {
            monitor_res_redist(ires, exitflag, iter,  phi2,  temp_phi2);
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

#endif /* RE_DISTANCE_H */

