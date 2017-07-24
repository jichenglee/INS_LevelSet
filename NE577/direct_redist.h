/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   direct_redist.h
 * Author: nsaini3
 *
 * Created on November 26, 2016, 5:21 PM
 */

#ifndef DIRECT_REDIST_H
#define DIRECT_REDIST_H

void direct_redist(elemsclr &sclr)
{
    vector< vector<double> > phi2(xelem, vector<double> (yelem,0.0));
    vector< vector<int> > tag(xelem, vector<int> (yelem,3));
    vector< vector<int> > int_cells;
    
    
    for(int i=0; i<xelem; i++)
    {
        for(int j=0; j<yelem; j++)
        {
            phi2[i][j] = sclr.phi[i][j][0];
        }
    }
    /****Calculate interface nodes*****/
    interface_nodes(sclr.phi, int_cells);
    
    /**Tag interface nodes as 1 and -1*/
    
    for(int k=0; k<int_cells.size(); k++)
    {
        int i = int_cells[k][0];
        int j = int_cells[k][1];
        
        if(sclr.phi[i][j][0] > 0.0)
        {
            tag[i][j] = 1;
        }
        else
        {
            tag[i][j] = -1;
        }
        
    }
    
    //cout<<"Interface cells: "<<int_cells.size()<<endl;
    
    
    
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            if(tag[i][j] != 1 && tag[i][j] != -1)
            {
                double multiplier;
                if(sclr.phi[i][j][0] > 0.0)
                {
                    multiplier = 1.0;                    
                }
                else if(sclr.phi[i][j][0] < 0.0)
                {
                    multiplier = -1.0;
                }
                
                if(multiplier > 0.0)
                {
                    double dist=10000.0;
                    double temp = 0.0;
                    for(int k=0; k < int_cells.size(); k++)
                    {
                        int x_in = int_cells[k][0];
                        int y_in = int_cells[k][1];
                        if(tag[x_in][y_in] == 1)
                        {
                            double x1 = xc[x_in][y_in];
                            double y1 = yc[x_in][y_in];
                            
                            double x2 = xc[i][j];
                            double y2 = yc[i][j];
                            
                            temp =sqrt(pow(x1-x2,2.0) + pow(y1-y2,2.0));
                            
                            if(temp < dist)
                            {
                                dist = temp;
                            }
                            
                        }
                    }                    
                    phi2[i][j] = dist;
                    
                }
                
                
                else
                {
                    double dist=10000.0;
                    double temp = 0.0;
                    for(int k=0; k < int_cells.size(); k++)
                    {
                        int x_in = int_cells[k][0];
                        int y_in = int_cells[k][1];
                        if(tag[x_in][y_in] == -1)
                        {
                            double x1 = xc[x_in][y_in];
                            double y1 = yc[x_in][y_in];
                            
                            double x2 = xc[i][j];
                            double y2 = yc[i][j];
                            
                            temp =sqrt(pow(x1-x2,2.0) + pow(y1-y2,2.0));
                            
                            if(temp < dist)
                            {
                                dist = temp;
                            }
                            
                        }
                        
                    }
                    phi2[i][j] = -dist;
                    
                }
                
            }
        }
    }
    
    
    
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            sclr.phi[i][j][0] = phi2[i][j];
        }
    }
    level_setBC(sclr.phi);
    
    
}

#endif /* DIRECT_REDIST_H */

