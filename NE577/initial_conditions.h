/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   initial_conditions.h
 * Author: nsaini3
 *
 * Created on November 24, 2016, 4:20 PM
 */

#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H

void initialize_phi(vector< vector<vector<double> > > &phi)
{
    for (int i=0; i<xelem; i++)
    {
        for (int j=0; j<yelem; j++)
        {
            phi[i][j][0] = sqrt(pow(xb_in - xc[i][j],2) + pow(yb_in - yc[i][j],2)) - rb_in;
        }
    }


    /******Zalesak Problem*****/
    if(case_tog == 4)
    {
        for (int i=0; i<xelem; i++)
        {
            for (int j=0; j<yelem; j++)
            {
                double a,b,c,d;
                a=(0.5-0.03)*100.0;
                b=(0.5+0.03)*100.0;
                c=(0.75-0.15)*100.0;
                d=(0.75+0.15-0.06)*100.0;
                if(xc[i][j]-a == 0.0 || xc[i][j]-b == 0.0 || yc[i][j]-c == 0.0 || yc[i][j]-d == 0.0)
                {
                    phi[i][j][0] =0.0;
                }
                else
                {

                    double dist[4];
                    dist[0] = sqrt(pow(xc[i][j]-a,2.0) + pow(yc[i][j]-c,2.0));
                    dist[1] = sqrt(pow(xc[i][j]-b,2.0) + pow(yc[i][j]-c,2.0));
                    dist[2] = sqrt(pow(xc[i][j]-b,2.0) + pow(yc[i][j]-d,2.0));
                    dist[3] = sqrt(pow(xc[i][j]-a,2.0) + pow(yc[i][j]-d,2.0));

                    int min = 0;
                    for (int i=1; i<4; i++)
                    {
                        if(dist[min] > dist[i])
                        {
                            min=i;
                        }
                    }

                    double base, height;
                    if(min == 0)
                    {
                        base = a - xc[i][j];
                        height = c - yc[i][j];
                    }
                    else if(min == 1)
                    {
                        base = b - xc[i][j];
                        height = c - yc[i][j];
                    }
                    else if(min == 2)
                    {
                        base = b - xc[i][j];
                        height = d - yc[i][j];
                    }
                    else if(min == 3)
                    {
                        base = a - xc[i][j];
                        height = d - yc[i][j];
                    }

                    double angle;
                        angle = atan2(height,base);
                        if(min == 0)
                        {
                            if(angle < PI/4.0 && angle > -3.0*PI/4.0)
                            {
                                phi[i][j][0] = base;//(height);
                            }
                            else
                            {
                                phi[i][j][0] = height;//(base);
                            }
                        }

                        else if(min == 2)
                        {
                            if(angle < PI/4.0 && angle > -3.0*PI/4.0)
                            {
                                phi[i][j][0] = -height;//(base);
                            }
                            else
                            {
                                phi[i][j][0] = -base;//(height);
                            }
                        }

                        else if(min == 1)
                        {
                            if(angle < 3.0*PI/4.0 && angle > -PI/4.0)
                            {
                                phi[i][j][0] = height;//(base);
                            }
                            else
                            {
                                phi[i][j][0] = -base;//(height);
                            }
                        }

                        else if (min == 3)
                        {
                            if(angle <= 3.0*PI/4.0 && angle >= -PI/4.0)
                            {
                                phi[i][j][0] = base;//(height);
                            }
                            else
                            {
                                phi[i][j][0] = -height;//(base);
                            }
                        }
                    //}

                }

            }
        }

        for(int i=0; i<xelem; i++)
        {
            for(int j=0; j<yelem; j++)
            {
                phi[i][j][0] = -phi[i][j][0];
            }
        }


        double temp_phi[xelem][yelem];
        for(int i=0; i<xelem; i++)
        {
            for(int j=0; j<yelem; j++)
            {
                temp_phi[i][j] = phi[i][j][0];
            }
        }

        for(int i=0; i<xelem; i++)
        {
            for(int j=0; j<yelem; j++)
            {
                double cir = sqrt(pow(xc[i][j]-50.0,2.0) + pow(yc[i][j]-75.0,2.0)) - 15.0;
                if(phi[i][j][0] < 0.0)
                {
                    phi[i][j][0] = max(cir,phi[i][j][0]);
                }
            }
        }
    }
}


void initialize(elemsclr &sclr)
{
    //calcp(sclr); //For the simple advection case
    //pressureBC(sclr.p);

    /*This initializes a single bubble*/
    initialize_phi(sclr.phi);
    level_setBC(sclr.phi);
    double eps = epsilon*max(xlen/(xelem-2), ylen/(yelem-2));

    vector< vector< vector<double> > > H(xelem, vector< vector<double> >(yelem, vector<double>(zelem,0.0)));

    heavy_func(H, sclr.phi, eps);


    /*****Bubble breakup case****/
    if (case_tog == 3)
    {
    double line = 2.5*2.0*rb_in;
    double temp_phi[xelem][yelem];
    for(int i=0; i<xelem; i++)
    {
        for(int j=0; j<yelem; j++)
        {
            temp_phi[i][j] = line - yc[i][j];
        }
    }

    for(int i=0; i<xelem; i++)
    {
        for(int j=0; j<yelem; j++)
        {
            if(fabs(temp_phi[i][j]) < fabs(sclr.phi[i][j][0]))
            {
                sclr.phi[i][j][0] = temp_phi[i][j];
            }
        }
    }
    level_setBC(sclr.phi);
    }

    /****Pressure*****/
    /*for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            if(yc[i][j] < line)
            {
                sclr.p[i][j][0] = rhof * fabs(gy)*fabs(line-yc[i][j]);
            }
            else
            {
                sclr.p[i][j][0] = 0.0;
            }
        }
    }
    double p_at_cen = rhof * fabs(gy)*fabs(line-yb_in);*/



    /*for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            if(H[i][j][0] < 1.0 && yc[i][j] < line)
            {
                sclr.p[i][j][0] = (1.0 - H[i][j][0])*p_at_cen;
            }
        }
    }
    /*********End of Bubble breakup initial conditions************/

    /***********************************************/

    /****Initializes density and viscosity based on phi****/
    find_density_visc(H, sclr.rho, sclr.mu);
    level_setBC(sclr.rho);
    level_setBC(sclr.mu);




}
#endif /* INITIAL_CONDITIONS_H */

