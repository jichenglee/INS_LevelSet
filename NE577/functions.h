/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   functions.h
 * Author: nsaini3
 *
 * Created on September 27, 2016, 7:44 PM
 */

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

void calcp(elemsclr &sclr)
{
    for(int i=0; i<xelem; i++)
    {
        for(int j=0; j<yelem; j++)
        {
            for(int k=0; k<zelem; k++)
            {
                sclr.p[i][j][k] = (-4800*xc[i][j] + 96)*nu;
            }
        }
        //cout<<endl;
    }
    
    /*for(int i=1; i < xelem-1; i++)
    {
        sclr.p[i][0][0]= sclr.p[i][1][0];
        sclr.p[i][yelem-1][0] = sclr.p[i][yelem-2][0];
    }*/
    
    /*for(int j=1; j < yelem-1; j++)
    {
        sclr.p[0][j][0] = sclr.p[xelem-2][j][0];
        sclr.p[xelem-1][j][0] = sclr.p [1][j][0];
    }*/
}

void monitor_res(vector<double> &ires, bool &exitflag, int iter, elemsclr sclr, vector< vector<vector<double> > > utemp,  vector< vector<vector<double> > > vtemp)
{
    vector<double> res(3,0.0);
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            res[0]=res[0] + pow(sclr.u[i][j][0]-utemp[i][j][0],2.0)*vol[i][j];
            res[1]=res[1] + pow(sclr.v[i][j][0]-vtemp[i][j][0],2.0)*vol[i][j];
        }
    }
    
    for(int i=0; i<3; i++)
    {
        res[i]=sqrt(res[i]);
    }
    
    
    if(iter == 0)
    {
        for(int i=0; i<3; i++)
        {
            ires[i]=res[i];
            //cout<<ires[i]<<endl;
        }
    }
    else
    {
        cout<<" U vel residual: "<<res[0]/ires[0]<<" V vel residual: "<<res[1]/ires[1];
        
        if(res[0]/ires[0] < tol && res[1]/ires[1] <  tol)
        {
            exitflag=true;
        }
    }
    
}


void timestep_calc(elemsclr &sclr, double &deltat, double &cfl)
{
    /*****Find the existing maximum cfl of the domain***/
    cfl=0.0;
    int reqi, reqj;
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            double surf_int = area[i][j][0][0]* (fabs(sclr.u[i][j][0]) + fabs(sclr.u[i-1][j][0]));
            surf_int = surf_int + area[i][j][1][1]* (fabs(sclr.v[i][j][0]) + fabs(sclr.v[i][j-1][0]));
            
            /*double surf_int = area[i][j][0][0]* (fabs(sclr.u[i][j][0] - sclr.u[i-1][j][0]));
            surf_int = surf_int + area[i][j][1][1]* (fabs(sclr.v[i][j][0] - sclr.v[i][j-1][0]));*/
            
            double temp_cfl = 0.0;
            if(cfl == 0.0)
            {
                cfl = surf_int*deltat/(area[i][j][0][0]*area[i][j][1][1]);
                reqi = i;
                reqj = j;
            }
            else
            {
                double temp_cfl = surf_int*deltat/(area[i][j][0][0]*area[i][j][1][1]);
                if(temp_cfl > cfl)
                {
                    cfl = temp_cfl;
                    reqi = i;
                    reqj = j;
                }
            }
        }
    }
    
    
    
    if(cfl < max_cfl)
    {
        cfl = cfl + 0.01*(max_cfl - cfl);
    }
    else
    {
        cfl = max_cfl;
    }
    
    double surf_int = area[reqi][reqj][0][0]* (fabs(sclr.u[reqi][reqj][0]) + fabs(sclr.u[reqi-1][reqj][0]));
    surf_int = surf_int + area[reqi][reqj][1][1]* (fabs(sclr.v[reqi][reqj][0]) + fabs(sclr.v[reqi][reqj-1][0]));
    
    /*double surf_int = area[reqi][reqj][0][0]* (fabs(sclr.u[reqi][reqj][0] - sclr.u[reqi-1][reqj][0]));
    surf_int = surf_int + area[reqi][reqj][1][1]* (fabs(sclr.v[reqi][reqj][0] - sclr.v[reqi][reqj-1][0]));*/
    
    deltat = cfl * area[reqi][reqj][0][0]*area[reqi][reqj][1][1]/surf_int;
    
    //deltat = advect_deltat;
    
}
#endif /* FUNCTIONS_H */

