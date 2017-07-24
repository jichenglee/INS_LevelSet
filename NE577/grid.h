/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   grid.h
 * Author: nsaini3
 *
 * Created on September 26, 2016, 4:34 PM
 */

#ifndef GRID_H
#define GRID_H

void gridread()
{
    double deltax=xlen/(xelem-2);
    double deltay=ylen/(yelem-2);
    
    //Assign coordinates
    for(int i=1; i<xnode-1; i++)
    {
        for(int j=1; j<ynode-1; j++)
        {
            x[i][j] = (i-1)*deltax;
            y[i][j] = (j-1)*deltay;
        }
    }
    
    //Generate ghost nodes
    //Add cells on both sides of x
    for(int j=1; j < ynode-1; j++)
    {
        x[0][j]=2*x[1][j]-x[2][j];
        y[0][j]=2*y[1][j]-y[2][j];
        x[xnode-1][j] = 2*x[xnode-2][j]-x[xnode-3][j];
        y[xnode-1][j] = 2*y[xnode-2][j]-y[xnode-3][j];
    }
    
    for(int i=1; i < xnode-1; i++)
    {
        x[i][0]=2*x[i][1]-x[i][2];
        y[i][0]=2*y[i][1]-y[i][2];
        x[i][ynode-1] = 2*x[i][ynode-2]-x[i][ynode-3];
        y[i][ynode-1] = 2*y[i][ynode-2]-y[i][ynode-3];
    }
    //Not required but assign proper coordinates to corner ghost nodes
    x[0][0] = 2*x[0][1]-x[0][2];
    y[0][0] = 2*y[0][1]-y[0][2];
    x[0][ynode-1] = 2*x[0][ynode-2]-x[0][ynode-3];
    y[0][ynode-1] = 2*y[0][ynode-2]-y[0][ynode-3];
    x[xnode-1][0] = 2*x[xnode-2][0]-x[xnode-3][0];
    y[xnode-1][0] = 2*y[xnode-2][0]-y[xnode-3][0];
    x[xnode-1][ynode-1] =  2*x[xnode-2][ynode-1]-x[xnode-3][ynode-1];
    y[xnode-1][ynode-1] =  2*y[xnode-2][ynode-1]-y[xnode-3][ynode-1];
    /*for(int i=0;i<xnode;i++)
    {
        for(int j=0;j<ynode;j++)
        {
            cout<<y[i][j]<<" ";
        }
        cout<<endl;
    }*/
    //Area components of faces of parallelogram
    for(int i=0; i<xelem; i++)
    {
        for(int j=0; j<yelem; j++)
        {
            area[i][j][0][0] = y[i+1][j+1] - y[i+1][j];
            area[i][j][0][1] = -(x[i+1][j+1] - x[i+1][j]);
            area[i][j][1][0] = y[i][j+1] - y[i+1][j+1];
            area[i][j][1][1] = -(x[i][j+1]-x[i+1][j+1]);
        }
    }
    
    for(int i=0; i<xelem; i++)
    {
        for(int j=0; j<yelem; j++)
        {
            //Heron's Formula for area of triangle (area of cell)
            double a=sqrt(pow(x[i][j]-x[i+1][j],2)+pow(y[i][j]-y[i+1][j],2));
            double b=sqrt(pow(x[i+1][j+1]-x[i+1][j],2)+pow(y[i+1][j+1]-y[i+1][j],2));
            double c=sqrt(pow(x[i][j]-x[i+1][j+1],2)+pow(y[i][j]-y[i+1][j+1],2));
            double s=(a+b+c)/2;
            vol[i][j]=sqrt(s*(s-a)*(s-b)*(s-c));

            double a1=sqrt(pow(x[i][j]-x[i][j+1],2)+pow(y[i][j]-y[i][j+1],2));
            double b1=sqrt(pow(x[i][j+1]-x[i+1][j+1],2)+pow(y[i][j+1]-y[i+1][j+1],2));
            double s1=(a1+b1+c)/2;
            vol[i][j]=vol[i][j]+sqrt(s1*(s1-a1)*(s1-b1)*(s1-c));

            //////Centroid of cell///////
            xc[i][j]=0.25*(x[i][j]+x[i+1][j]+x[i+1][j+1]+x[i][j+1]);
            yc[i][j]=0.25*(y[i][j]+y[i+1][j]+y[i+1][j+1]+y[i][j+1]);
        }
    }
}

#endif /* GRID_H */

