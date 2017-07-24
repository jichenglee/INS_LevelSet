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

void flux(double u1R, double u1L, double u1T, double u1B, double u2R, double u2L, double u2T, double u2B, int dir, double &flux1, double &flux2)
{
    if(dir ==1)
    {
         flux1 = pow(u1R + u1L,2.0)/4.0;
         flux2 = (u2R +u2L)*0.5 * (u2T + u2B)*0.5;
    }
    else if (dir == 2)
    {
        flux1 = (u1R +u1L)*0.5 * (u1T + u1B)*0.5;
        flux2 =  pow(u2T + u2B,2.0)/4.0;
    }
}




void quick(double u1R, double u1L, double u1T, double u1B, double u2R, double u2L, double u2T, double u2B, double u11, double u22, int dir, double &flux1, double &flux2)
{
    if(dir==1)
    {
        if(u1L+u1R > 0.0)
        {
            flux1 = 0.125*(-u11 + 6*u1L + 3*u1R);
        }
        else
        {
            flux1 = 0.125*(-u22 + 6*u1R + 3*u1L);
        }

        flux2 = (u2R +u2L)*0.5 * (u2T + u2B)*0.5;
    }
    else if(dir ==2)
    {
        flux1 = (u1R +u1L)*0.5 * (u1T + u1B)*0.5;
        if(u2B+u2T > 0.0)
        {
            flux2 = 0.125*(-u11 + 6*u2B + 3*u2T);
        }
        else
        {
            flux2 = 0.125*(-u22 + 6*u2T + 3*u2B);
        }
    }
}

void rhscalc(elemsclr &sclr, vector< vector<double> > &rhsx, vector< vector<double> > &rhsy, int iter, bool exitflag)
{





    //***Calculate contribution from advection
    //Note for index i,j the CV under consideration is the CV between i,j and i+1,j
    vector< vector<double> > advx(xelem, vector<double>(yelem,0.0));
    vector< vector<double> > advy(xelem, vector<double>(yelem,0.0));

    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            //Calculate the velocities at the edge centers of CV
            double u1R, u1L;
            double u1T, u1B;
            double u2R, u2L;
            double u2T, u2B;

            double u1LL, u1RR;

            u1R = sclr.u[i+1][j][0];
            u1L = sclr.u[i][j][0];

            u1T = sclr.v[i+1][j][0];
            u1B = sclr.v[i+1][j-1][0];

            u2R = sclr.v[i+1][j][0];
            u2L = sclr.v[i][j][0];

            u2T = sclr.u[i][j+1][0];
            u2B = sclr.u[i][j][0];

            if(i==0)
            {
                u1LL = sclr.u[xelem-3][j][0];
                u1RR = sclr.u[i+2][j][0];
            }
            else if(i==xelem-2)
            {
                u1LL = sclr.u[i-1][j][0];
                u1RR = sclr.u[2][j][0];
            }
            else
            {
                u1LL = sclr.u[i-1][j][0];
                u1RR = sclr.u[i+2][j][0];
            }


            double flux1, flux2;
            flux(u1R, u1L, u1T, u1B, u2R, u2L, u2T, u2B, 1, flux1, flux2);
            //quick(u1R, u1L, u1T, u1B, u2R, u2L, u2T, u2B, u1LL, u1RR, 1, flux1, flux2);
            advx[i][j]=(area[i][j][0][0]*flux1 + area[i][j][1][1]*flux2)/vol[i][j];

        }
    }

    #pragma omp parallel for schedule(dynamic)
    for (int i=1; i< xelem-1 ;i++)
    {
        for (int j=0; j< yelem-1; j++)
        {
            double v1R, v1L;
            double v1T, v1B;
            double v2R, v2L;
            double v2T, v2B;

            double v2TT, v2BB;

            v1R = sclr.v[i+1][j][0];
            v1L = sclr.v[i][j][0];

            v1T = sclr.u[i][j+1][0];
            v1B = sclr.u[i][j][0];

            v2R = sclr.u[i][j+1][0];
            v2L = sclr.u[i-1][j+1][0];

            v2T = sclr.v[i][j+1][0];
            v2B = sclr.v[i][j][0];

            if(j==0)
            {
                v2TT = sclr.v[i][j+2][0];
                v2BB = 0.0;
            }
            else if(j==yelem-2)
            {
                v2TT = 0.0;
                v2BB = sclr.v[i][j-1][0];
            }
            else
            {
                v2TT = sclr.v[i][j+2][0];
                v2BB = sclr.v[i][j-1][0];
            }


            double flux1, flux2;
            flux(v1R, v1L, v1T, v1B, v2R, v2L, v2T, v2B, 2, flux1, flux2);
            //quick(v1R, v1L, v1T, v1B, v2R, v2L, v2T, v2B, v2BB, v2TT, 2, flux1, flux2);
            advy[i][j]=(area[i][j][0][0]*flux1 + area[i][j][1][1]*flux2)/vol[i][j];
        }
    }


    //***Calculate contribution from diffusion
    vector< vector<double> > diffx(xelem, vector<double>(yelem,0.0));
    vector< vector<double> > diffy(xelem, vector<double>(yelem,0.0));
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            //Simple taylor series approximation of second diff is used
            double u, uR, uL, uT, uB;
            if (i == 0)
            {
                u = sclr.u[i][j][0];
                uR = sclr.u[i+1][j][0];
                if(x_bound == 1)
                {
                    uL = -sclr.u[i+1][j][0];
                }
                else if(x_bound == 2)
                {
                    uL = sclr.u[i+1][j][0];
                }
                else if(x_bound == 3)
                {
                    uL = sclr.u[xelem-3][j][0];
                }

                uT = sclr.u[i][j+1][0];
                uB = sclr.u[i][j-1][0];
            }
            else
            {
                u = sclr.u[i][j][0];
                uR = sclr.u[i+1][j][0];
                uL = sclr.u[i-1][j][0];
                uT = sclr.u[i][j+1][0];
                uB = sclr.u[i][j-1][0];
            }

            diffx[i][j] = ((sclr.mu[i+1][j][0] + sclr.mu[i][j][0])/(sclr.rho[i+1][j][0] + sclr.rho[i][j][0]))*((uR + uL - 2.0*u)/pow(area[i][j][1][1],2.0) + (uT + uB - 2.0*u)/pow(area[i][j][0][0],2.0));
        }
    }

    #pragma omp parallel for schedule(dynamic)
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=0; j<yelem-1; j++)
        {
            //Simple taylor series approximation of second diff is used
            double v, vR, vL, vT, vB;
            if (j == 0)
            {
                v = sclr.v[i][j][0];
                vR = sclr.v[i+1][j][0];
                vL = sclr.v[i-1][j][0];
                vT = sclr.v[i][j+1][0];
                if(y_bound == 1)
                {
                    vB = -vT;
                    //vB = 0.0;
                }
                else if(y_bound == 2)
                {
                    vB = vT;
                }
                else if(y_bound == 3)
                {
                    vB = sclr.v[i][yelem-3][0];
                }

            }
            else
            {
                v = sclr.v[i][j][0];
                vR = sclr.v[i+1][j][0];
                vL = sclr.v[i-1][j][0];
                vT = sclr.v[i][j+1][0];
                vB = sclr.v[i][j-1][0];
            }
            diffy[i][j] = ((sclr.mu[i][j+1][0] + sclr.mu[i][j][0])/(sclr.rho[i][j+1][0] + sclr.rho[i][j][0]))*((vR + vL - 2.0*v)/pow(area[i][j][1][1],2.0) + (vT + vB - 2.0*v)/pow(area[i][j][0][0],2.0));
        }
    }





    #pragma omp parallel for schedule(dynamic)
    for(int j=1; j<yelem-1; j++)
    {
        for(int i=1; i<xelem-1; i++)
        {
            rhsx[i][j]=diffx[i][j]-(advx[i][j]-advx[i-1][j]);
            rhsy[i][j]=diffy[i][j]-(advy[i][j]-advy[i][j-1]);

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

