/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   heavy_delta.h
 * Author: nsaini3
 *
 * Created on November 11, 2016, 7:58 PM
 */

#ifndef HEAVY_DELTA_H
#define HEAVY_DELTA_H


void find_density_visc(vector< vector< vector<double> > > &H, vector< vector< vector<double> > > &rho, vector< vector< vector<double> > > &mu)
{
    for(int i=0; i<xelem; i++)
    {
        for(int j=0; j<yelem; j++)
        {
            //cout<<i<<" "<<rhog + (rhof -rhog)*H[i][j][0]<<" "<<H[i][j][0]<< endl;
            rho[i][j][0] = rhog + (rhof -rhog)*H[i][j][0];
            mu[i][j][0] = mug + (muf -mug)*H[i][j][0];
        }
    }
}

void heavy_func(vector< vector< vector<double> > > &H, vector< vector< vector<double> > > &phi, double eps)
{
    for(int i=1; i<xelem-1; i++)
        {
            for(int j=1; j<yelem-1; j++)
            {
                if(phi[i][j][0] < -eps)
                {
                    H[i][j][0] = 0.0;
                }
                else if(phi[i][j][0] > eps)
                {
                    H[i][j][0] = 1.0;
                }
                else
                {
                    //H[i][j][0] = 1.0/(1.0 + exp(-3.0*sclr.phi[i][j][0]/eps));
                    H[i][j][0] = 0.5 + phi[i][j][0]/(2.0*eps) + (1.0/(2.0*PI)*sin(PI * phi[i][j][0]/eps));
                }
                //cout<<H[i][j][0]<<" ";
            }
            //cout<<endl;
        }
        //exit(0);
        level_setBC(H);
}

void delta_func(vector< vector< vector<double> > > &delta, vector< vector< vector<double> > > &phi, double eps)
{
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            if(fabs(phi[i][j][0]) > eps)
            {
                delta[i][j][0] = 0.0;
            }
            else
            {
                //delta[i][j][0] = 3.0*exp(-3.0*sclr.phi[i][j][0]/eps)/(1.0 * pow(1.0 + exp(-3.0*sclr.phi[i][j][0]/eps),2.0));
                delta[i][j][0] = (1.0/2.0*eps) * (1.0 + cos(PI * phi[i][j][0]/eps));
            }
            //cout<<delta[i][j][0]<<" ";
        }
        //cout<<endl;
    }
    //exit(0);
}

void grad_func(vector< vector< vector<double> > > &grad_phi, vector< vector< vector<double> > > &phi)
{
    vector< vector< vector<double> > > grad_phix(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0))); 
     vector< vector< vector<double> > > grad_phiy(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0))); 
     
     vector< vector< vector<double> > > phiRface(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));
     vector< vector< vector<double> > > phiTface(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));
     for(int i=0; i<xelem-1; i++)
     {
         for(int j=0; j<yelem-1; j++)
         {
             phiRface[i][j][0] = 0.5*(phi[i+1][j][0] + phi[i][j][0]);
             phiTface[i][j][0] = 0.5*(phi[i][j+1][0] + phi[i][j][0]);
         }
     }
     
     for(int j=1; j<yelem-1; j++)
     {
         for(int i=1; i<xelem-1; i++)
         {
             grad_phix[i][j][0] = (phiRface[i][j][0] - phiRface[i-1][j][0])/area[i][j][1][1];
             grad_phiy[i][j][0] = (phiRface[i][j][0] - phiRface[i][j-1][0])/area[i][j][0][0];
             //if(delta[i][j][0] != 0.0){
             //cout<<grad_phix[i][j][0]<<" ";
             //}
         }
         //cout<<endl;
     }
     //exit(0);
     /*Need to impose BC for grad_phix and grad_phiy*/
     grad_level_setBC(grad_phix);
     grad_level_setBC(grad_phiy);
     
     for(int i=0; i<xelem; i++)
     {
         for(int j=0; j<yelem; j++)
         {
             grad_phi[i][j][0] = sqrt(pow(grad_phix[i][j][0],2.0) + pow(grad_phiy[i][j][0],2.0));
         }
     }
}


void vol_contraint(vector< vector< vector<double> > > &phi2, vector< vector< vector<double> > > &phi, vector< vector< vector<double> > > &grad_phi, vector< vector< vector<double> > > &delta, double deltat)
{
    for(int i=0; i<xelem; i++)
    {
        for(int j=0; j<yelem; j++)
        {
            if(delta[i][j][0] != 0.0)
            {
                phi2[i][j][0] = phi[i][j][0];
            }
        }
    }
}

#endif /* HEAVY_DELTA_H */

