/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   rhs.h
 * Author: nsaini3
 *
 * Created on September 27, 2016, 4:03 PM
 */

#ifndef RHS_H
#define RHS_H

double flux(double velcen, double velplus, double velplus2, double velminus, int dir)
{
    int uinc=0;
    int vinc=0;
    int winc=0;
    if(dir==1)
    {
        uinc=1;
    }
    else if(dir==2)
    {
        vinc=1;
    }
    else if(dir == 3)
    {
        winc=1;
    }
    
    double velface1 =  
}


void rhscalc(vector< vector< vector<double> > > vel, vector< vector<double> > &rhs, int iter, int dir, vector< vector< vector<double> > > othervel)
{
    int uinc=0;
    int vinc=0;
    int winc=0;
    if(dir==1)
    {
        uinc=1;
    }
    else if(dir==2)
    {
        vinc=1;
    }
    else if(dir == 3)
    {
        winc=1;
    }
    //***Calculate contribution from advection
    //Note for index i,j the CV under consideration is the CV between i,j and i+1,j
    vector< vector<double> > adv(xelem, vector<double>(yelem,0.0));
    for(int i=vinc; i<xelem-1; i++)
    {
        for(int j=uinc; j<yelem-1; j++)
        {
            //Calculate the velocities at the edge centers of CV
            double velcen, velplus, velplus2, velminus;
            
            
            if(i == 0)
            {
                velcen = vel[i][j][0];
                velplus = vel[i+uinc][j+vinc][0+winc];
                velplus2 =  vel[i+2*uinc][j+2*vinc][0+2*winc];
                velminus = vel[xelem-3][j][0];
            }
            else if (j ==0)
            {
                velcen = 0.0;
                velplus = vel[i+uinc][j+vinc][0+winc];
                velplus2 =  vel[i+2*uinc][j+2*vinc][0+2*winc];
                velminus = 0.0;
            }
            else if (i == xelem-2 && dir == 1)
            {
                velcen = vel[i][j][0];
                velplus = vel[i+uinc][j+vinc][0+winc];
                velplus2 =  vel[1][j][0];
                velminus = vel[i-uinc][j-vinc][0-winc];
            }
            else if (j == yelem-2 && dir == 2)
            {
                velcen = 0.0;
                velplus = 0.0;
                velplus2 = 0.0;
                velminus = vel[i-uinc][j-vinc][0-winc];
            }
            else
            {
                velcen = vel[i][j][0];
                velplus = vel[i+uinc][j+vinc][0+winc];
                velplus2 =  vel[i+2*uinc][j+2*vinc][0+2*winc];
                velminus = vel[i-uinc][j-vinc][0-winc];
            }
                
            double othervel1, othervel2;
            if(dir ==1)
            {
                othervel1 = othervel[i+1][j];
                othervel2 = othervel[i+1][j-1];
            }
            else if(dir == 2)
            {
                othervel1 = othervel[i][j];
                othervel2 = othervel[i+1][j-1];
            }
            double vel1=0.5*(vel[i+uinc][j+vinc][0+winc]+vel[i][j][0]);
            double vel2=0.5*(vel[i][j][0]+vel[i-uinc][j-vinc][0-winc]);
            
            //double v1=0.5*(v[i+1][j]+v[i][j]);
            //double v2=0.5*(v[i+1][j-1]+v[i][j-1]);            
            adv[i][j]=area[i][j][0][0]*(vel1*vel1-vel2*vel2)/vol[i][j];   
            //cout<<area[i][j][0][0]<<" "<<vol[i][j]<<" "<<advx[i][j]<<" "<<area[i][j][1][1]<<endl;
            //exit(0);
        }
    }
    
    //***Calculate contribution from diffusion
    vector< vector<double> > diff(xelem, vector<double>(yelem,0.0));
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            //Simple taylor series approximation of second diff is used
            diff[i][j] = (((vel[i+uinc][j+vinc][0+winc]-2.0*vel[i][j][0]+vel[i-uinc][j-vinc][0-winc])/pow(area[i][j][1][1],2.0)) + ((vel[i+vinc][j+uinc][0]-2.0*vel[i][j][0]+vel[i-vinc][j-uinc][0])/pow(area[i][j][0][0],2.0)))*nu;
        }
    }
    
    for(int j=1; j<yelem-1; j++)
    {
        for(int i=1; i<xelem-1; i++)
        {
            rhs[i][j]=diff[i][j]-adv[i][j];
            if(iter == 1)
            {
                //cout<<diffx[i][j]<<" ";
            }
            
        }
        if(iter == 1){
        //cout<<endl;
        }
    }
    
}

#endif /* RHS_H */

