/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   rhs_bub.h
 * Author: nsaini3
 *
 * Created on October 20, 2016, 2:45 PM
 */

#ifndef RHS_BUB_H
#define RHS_BUB_H

double minmod(double a, double b)
{
    double result;

    if(fabs(a) < fabs(b))
    {
        result = a;
    }
    else
    {
        result = b;
    }

    return result;
}

double phi_at_edge(double uR, double uL, double phiR, double phiL, double phi2R, double phi2L)
{
    double phi_face;

    double a = 0.5*(uR + uL); // wave speed at cell interface

    if(a > 0)
    {
        double dxplus = phiR - phiL;
        double dxminus = phiL - phi2L;

        phi_face = phiL + 0.5* minmod(dxplus, dxminus);
    }
    else
    {
        double dxplus = phi2R - phiR;
        double dxminus = phiR - phiL;

        phi_face = phiR - 0.5* minmod(dxplus, dxminus);
    }

    return phi_face;
}

double weno(double a, double b, double c, double d)
{
    double IS0 = 13.0*pow(a-b,2.0) + 3.0*pow(a-3.0*b,2.0);
    double IS1 = 13.0*pow(b-c,2.0) + 3.0*pow(b+c,2.0);
    double IS2 = 13.0*pow(c-d,2.0) + 3.0*pow(3.0*c-d,2.0);

    double eps=10e-6;

    double alpha0 = 1.0/pow(eps + IS0,2.0);
    double alpha1 = 6.0/pow(eps + IS1,2.0);
    double alpha2 = 3.0/pow(eps + IS2,2.0);

    double w0 = alpha0/(alpha0+alpha1+alpha2);
    double w2 = alpha2/(alpha0+alpha1+alpha2);

    double result = (1.0/3.0)*w0*(a-2*b+c) + (1.0/6.0)*(w2-0.5)*(b-2.0*c+d);

    return result;
}

double delplus(double phi[6], int index)
{
    double result = phi[index+1] - phi[index];
    return result;
}

double delminus(double phi[6], int index)
{
    double result = phi[index] - phi[index-1];
    return result;
}

double hj_weno(double vel1, double vel2, double phicen, double phi1, double phi11, double phi111, double phi2, double phi22, double phi222, double len)
{
    double phi[7];

    phi[0]=phi111;
    phi[1]=phi11;
    phi[2]=phi1;
    phi[3]=phicen;
    phi[4]=phi2;
    phi[5]=phi22;
    phi[6]=phi222;

    /*phi[3] = 0.5*(phicen + phi2);
    phi[2] = 0.5*(phi1 + phicen);
    phi[1] = 0.5*(phi11 + phi1);
    phi[0] = 0.5*(phi111 + phi11);
    phi[4] = 0.5*(phi2 + phi22);
    phi[5] = 0.5*(phi22 + phi222);*/


    double i=3;
     double a,b,c,d;
    double phiplus = (1.0/12.0)*(-delplus(phi,i-2)/len + 7*delplus(phi,i-1)/len + 7*delplus(phi,i)/len - delplus(phi,i+1)/len);

    a = (delminus(phi,i+3) - delminus(phi,i+2))/len;
    b = (delminus(phi,i+2) - delminus(phi,i+1))/len;
    c = (delminus(phi,i+1) - delminus(phi,i))/len;
    a = (delminus(phi,i) - delminus(phi,i-1))/len;
    phiplus = phiplus - weno(a,b,c,d);

    double phiminus = (1.0/12.0)*(-delplus(phi,i-2)/len + 7*delplus(phi,i-1)/len + 7*delplus(phi,i)/len - delplus(phi,i+1)/len);
    a = (delminus(phi,i-1) - delminus(phi,i-2))/len;
    b = (delminus(phi,i) - delminus(phi,i-1))/len;
    c = (delminus(phi,i+1) - delminus(phi,i))/len;
    a = (delminus(phi,i+2) - delminus(phi,i+1))/len;
    phiminus = phiminus - weno(a,b,c,d);

    double phi_face=phiminus;

    double wavespeed = vel2;

    if(wavespeed >= 0.0)
    {
        phi_face=phiminus;
    }
    else
    {
        phi_face=phiplus;
    }

    return phi_face;
}


void rhs_bub(vector< vector<double> > &rhsx, vector< vector<double> > &rhsy, vector< vector<vector<double> > > ucen, vector< vector<vector<double> > > vcen, vector< vector<vector<double> > > phi)
{
    /*First of, find the value of phi at cell edges and store it*/
    vector< vector<double> > gradx(xelem, vector<double> (yelem,0.0));
    vector< vector<double> > grady(xelem, vector<double> (yelem,0.0));
    #pragma omp parallel for schedule(dynamic)
    for (int i=0; i < xelem-1; i++)
    {
        for(int j=1; j < yelem-1; j++)
        {
           //cout<<i<<" "<<j<<endl;
            double phicen, phiL, phiLL, phiLLL, phiR, phiRR, phiRRR;
            if (i == 0)
            {
                phicen = phi[i][j][0];
                if(x_bound == 1 || x_bound == 2)
                {
                    phiL = phicen;
                }
                else if(x_bound == 3)
                {
                    phiL=phi[xelem-3][j][0];
                }

                phiLL=phi[xelem-4][j][0];
                phiLLL=phi[xelem-5][j][0];
                phiR=phi[i+1][j][0];
                phiRR=phi[i+2][j][0];
                phiRRR=phi[i+3][j][0];

            }

            else if (i == xelem-2)
            {
                phicen = phi[i][j][0];
                phiL=phi[i-1][j][0];
                phiLL=phi[i-2][j][0];
                phiLLL=phi[i-3][j][0];
                phiR=phi[i+1][j][0];
                if(x_bound == 1 || x_bound == 2)
                {
                    phiRR=phiR;
                }
                else if(x_bound == 3)
                {
                    phiRR=phi[2][j][0];
                }

                phiRRR=phi[3][j][0];
            }

            else if(i == 1)
            {
                phicen = phi[i][j][0];
                phiL=phi[i-1][j][0];
                phiLL=phi[xelem-3][j][0];
                phiLLL=phi[xelem-4][j][0];
                phiR=phi[i+1][j][0];
                phiRR=phi[i+2][j][0];
                phiRRR=phi[i+3][j][0];
            }
            else if(i == 2)
            {
                phicen = phi[i][j][0];
                phiL=phi[i-1][j][0];
                phiLL=phi[i-2][j][0];
                phiLLL=phi[xelem-3][j][0];
                phiR=phi[i+1][j][0];
                phiRR=phi[i+2][j][0];
                phiRRR=phi[i+3][j][0];
            }

            else if (i == xelem-3)
            {
                phicen = phi[i][j][0];
                phiL=phi[i-1][j][0];
                phiLL=phi[i-2][j][0];
                phiLLL=phi[i-3][j][0];
                phiR=phi[i+1][j][0];
                phiRR=phi[i+2][j][0];
                phiRRR=phi[2][j][0];
            }

            else
            {
                phicen = phi[i][j][0];
                phiL=phi[i-1][j][0];
                phiLL=phi[i-2][j][0];
                phiLLL=phi[i-3][j][0];
                phiR=phi[i+1][j][0];
                phiRR=phi[i+2][j][0];
                phiRRR=phi[i+3][j][0];
            }


            if(bub_conv_scheme == 1)
            {
                gradx[i][j] = phi_at_edge(ucen[i+1][j][0], ucen[i][j][0], phiR, phicen, phiRR, phiL);
            }
            else if(bub_conv_scheme == 2)
            {
                gradx[i][j] = hj_weno(ucen[i+1][j][0], ucen[i][j][0], phicen, phiL, phiLL, phiLLL, phiR, phiRR, phiRRR, area[i][j][1][1]);
            }
        }


    }

    #pragma omp parallel for schedule(dynamic)
    for (int i=1; i < xelem-1; i++)
    {
        for(int j=0; j < yelem-1; j++)
        {
           //cout<<i<<" "<<j<<endl;
            double phicen, phiT, phiTT, phiTTT, phiB, phiBB, phiBBB;
            if (j == 0) ///gradient must be 0
            {
                phicen = phi[i][j][0];
                if(y_bound ==1 || y_bound ==2)
                {
                    phiB = phicen;
                }
                else if(y_bound == 3)
                {
                    phiB = phi[i][yelem-3][0];
                }
                phiB=phicen;
                phiBB=phicen;
                phiBBB=phicen;
                phiT=phi[i][j+1][0];
                phiTT=phi[i][j+2][0];
                phiTTT=phi[i][j+3][0];

            }

            else if (j == yelem-2)
            {
                phicen = phi[i][j][0];
                phiB=phi[i][j-1][0];
                phiBB=phi[i][j-2][0];
                phiBBB=phi[i][j-3][0];
                phiT=phi[i][j+1][0];
                if(y_bound == 1 || y_bound == 2)
                {
                    phiTT = phiT;
                }
                else if(y_bound == 3)
                {
                    phiTT = phi[i][2][0];
                }
                phiTT=phicen;
                phiTTT=phicen;
            }

            else if(j == 1)
            {
                phicen = phi[i][j][0];
                phiB=phicen;
                phiBB=phicen;
                phiBBB=phicen;
                phiT=phi[i][j+1][0];
                phiTT=phi[i][j+2][0];
                phiTTT=phi[i][j+3][0];
            }
            else if(j == 2)
            {
                phicen = phi[i][j][0];
                phiB=phi[i][j-1][0];
                phiBB=phiB;
                phiBBB=phiB;
                phiT=phi[i][j+1][0];
                phiTT=phi[i][j+2][0];
                phiTTT=phi[i][j+3][0];
            }

            else if (j == yelem-3)
            {
                phicen = phi[i][j][0];
                phiB=phi[i][j-1][0];
                phiBB=phi[i][j-2][0];
                phiBBB=phi[i][j-3][0];
                phiT=phi[i][j+1][0];
                phiTT=phiT;
                phiTTT=phiT;
            }

            else
            {
                phicen = phi[i][j][0];
                phiB=phi[i][j-1][0];
                phiBB=phi[i][j-2][0];
                phiBBB=phi[i][j-3][0];
                phiT=phi[i][j+1][0];
                phiTT=phi[i][j+2][0];
                phiTTT=phi[i][j+3][0];
            }


            if(bub_conv_scheme == 1)
            {
                grady[i][j] = phi_at_edge(vcen[i][j+1][0], vcen[i][j][0], phiT, phicen, phiTT, phiB);
            }
            else if(bub_conv_scheme == 2)
            {
                grady[i][j] = hj_weno(vcen[i][j+1][0], vcen[i][j][0], phicen, phiB, phiBB, phiBBB, phiT, phiTT, phiTTT, area[i][j][0][0]);
            }
        }


    }

    /*Now calculate the fluxes at the faces. Note v velocity is not yet regarded*/
    #pragma omp parallel for schedule(dynamic)
    for (int i=1; i<xelem-1; i++)
    {
        for (int j=1; j< yelem-1; j++)
        {
            if(bub_conv_scheme == 2)
            {
                rhsx[i][j] = -ucen[i][j][0]*(gradx[i][j]);
                rhsy[i][j] = -vcen[i][j][0]*(grady[i][j]);
            }
            else if(bub_conv_scheme == 1)
            {
                rhsx[i][j] = -(0.5*(ucen[i][j][0] + ucen[i+1][j][0])*gradx[i][j]*area[i][j][0][0] - 0.5*(ucen[i][j][0] + ucen[i-1][j][0])*gradx[i-1][j]*area[i-1][j][0][0])/vol[i][j];
                rhsy[i][j] = -(0.5*(vcen[i][j][0] + vcen[i][j+1][0])*grady[i][j]*area[i][j][1][1] - 0.5*(vcen[i][j][0] + vcen[i][j-1][0])*grady[i][j-1]*area[i][j-1][0][0])/vol[i][j];
            }



        }
        //exit(0);
    }
}

#endif /* RHS_BUB_H */

