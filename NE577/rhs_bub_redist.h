/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   rhs_bub_redist.h
 * Author: nsaini3
 *
 * Created on November 7, 2016, 8:30 PM
 */

#ifndef RHS_BUB_REDIST_H
#define RHS_BUB_REDIST_H

void rhs_redist(vector< vector<double> > &rhsx, vector< vector<double> > &rhsy, vector< vector<vector<double> > > ucen, vector< vector<vector<double> > > vcen, vector< vector<vector<double> > > phi)
{
    /*First of, find the value of phi at cell edges and store it*/
    vector< vector<double> > gradx(xelem, vector<double> (yelem,0.0));
    vector< vector<double> > grady(xelem, vector<double> (yelem,0.0));
    for (int i=0; i < xelem-1; i++)
    {
        for(int j=1; j < yelem-1; j++)
        {
           //cout<<i<<" "<<j<<endl;
            double phicen, phiL, phiLL, phiLLL, phiR, phiRR, phiRRR;
            if (i == 0)
            {
                phicen = phi[i][j][0];
                phiL=phi[xelem-3][j][0];
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
                phiRR=phi[2][j][0];
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


            gradx[i][j] = phi_at_edge(ucen[i+1][j][0], ucen[i][j][0], phiR, phicen, phiRR, phiL);

        }


    }


    for (int i=1; i < xelem-1; i++)
    {
        for(int j=0; j < yelem-1; j++)
        {
           //cout<<i<<" "<<j<<endl;
            double phicen, phiT, phiTT, phiTTT, phiB, phiBB, phiBBB;
            if (j == 0) ///gradient must be 0
            {
                phicen = phi[i][j+1][0];
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
                phiT=phicen;
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


            grady[i][j] = phi_at_edge(vcen[i][j+1][0], vcen[i][j][0], phiT, phicen, phiTT, phiB);

        }
    }


    for (int i=1; i<xelem-1; i++)
    {
        for (int j=1; j< yelem-1; j++)
        {
            rhsx[i][j] = -(0.5*(ucen[i][j][0] + ucen[i+1][j][0])*gradx[i][j]*area[i][j][0][0] - 0.5*(ucen[i][j][0] + ucen[i-1][j][0])*gradx[i-1][j]*area[i-1][j][0][0])/vol[i][j];
            rhsy[i][j] = -(0.5*(vcen[i][j][0] + vcen[i][j+1][0])*grady[i][j]*area[i][j][1][1] - 0.5*(vcen[i][j][0] + vcen[i][j-1][0])*grady[i][j-1]*area[i][j-1][0][0])/vol[i][j];

        }
    }
}


double dplus(double phi[5], int index)
{
    double result = phi[index+1] - phi[index];
    return result;
}


double dminus(double phi[5], int index)
{
    double result = phi[index] - phi[index-1];
    return result;
}

double signof(double a)
{
    double result = a/fabs(a);
    return result;
}

void rhs_redist2(vector< vector<double> > &rhs, vector< vector<vector<double> > > phi2, vector< vector<vector<double> > > phi)
{
    double dxbarplus[xelem][yelem];
    double dxbarminus[xelem][yelem];
    double dxbar[xelem][yelem];

    #pragma omp parallel for schedule(dynamic)
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            double temp_phi[5];
            if(i == 1)
            {
                temp_phi[2] = phi2[i][j][0];
                temp_phi[3] = phi2[i+1][j][0];
                temp_phi[4] = phi2[i+2][j][0];
                temp_phi[1] = phi2[i-1][j][0];
                if(x_bound == 1 || x_bound ==2)
                {
                    temp_phi[0] = temp_phi[1];
                }
                else if(x_bound ==3)
                {
                    temp_phi[0] = phi2[xelem-3][j][0];
                }

            }
            else if (i == xelem-2)
            {
                temp_phi[2] = phi2[i][j][0];
                temp_phi[3] = phi2[i+1][j][0];
                if(x_bound == 1 || x_bound ==2)
                {
                    temp_phi[4] = temp_phi[3];
                }
                else if(x_bound ==3)
                {
                    temp_phi[4] = phi2[2][j][0];
                }

                temp_phi[1] = phi2[i-1][j][0];
                temp_phi[0] = phi2[i-2][j][0];
            }
            else
            {
                temp_phi[2] = phi2[i][j][0];
                temp_phi[3] = phi2[i+1][j][0];
                temp_phi[4] = phi2[i+2][j][0];
                temp_phi[1] = phi2[i-1][j][0];
                temp_phi[0] = phi2[i-2][j][0];
            }

            dxbarplus[i][j] = dplus(temp_phi, 2) -0.5*minmod(dplus(temp_phi,2) - dplus(temp_phi,1) , dplus(temp_phi,3) - dplus(temp_phi,2));
            dxbarminus[i][j] = dminus(temp_phi, 2) +0.5*minmod(dplus(temp_phi,2) - dplus(temp_phi,1) , dplus(temp_phi,1) - dplus(temp_phi,0));

           if(signof(phi[i][j][0])*dplus(temp_phi,2) < 0.0  && signof(phi[i][j][0])*dminus(temp_phi,2) < -signof(phi[i][j][0])*dplus(temp_phi,2))
           {
               dxbar[i][j] = dxbarplus[i][j];
           }
           else if(signof(phi[i][j][0])*dminus(temp_phi,2) > 0.0  && signof(phi[i][j][0])*dplus(temp_phi,2) > -signof(phi[i][j][0])*dminus(temp_phi,2))
           {
               dxbar[i][j] = dxbarminus[i][j];
           }
           else
           {
               dxbar[i][j] = 0.5*(dxbarplus[i][j] + dxbarminus[i][j]);
           }
        }

    }

    double dybarplus[xelem][yelem];
    double dybarminus[xelem][yelem];
    double dybar[xelem][yelem];

    #pragma omp parallel for schedule(dynamic)
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            double temp_phi[5];
            if(j == 1)
            {
                temp_phi[2] = phi2[i][j][0];
                temp_phi[3] = phi2[i][j+1][0];
                temp_phi[4] = phi2[i][j+2][0];
                temp_phi[1] = phi2[i][j-1][0];
                if(y_bound == 1 || y_bound ==2)
                {
                    temp_phi[0] = phi2[i][j-1][0];
                }
                else if(y_bound ==3)
                {
                    temp_phi[0] = phi2[i][yelem-3][0];
                }

            }
            else if (j == yelem-2)
            {
                temp_phi[2] = phi2[i][j][0];
                temp_phi[3] = phi2[i][j+1][0];
                if(y_bound == 1 || y_bound ==2)
                {
                    temp_phi[4] = phi2[i][j+1][0];
                }
                else if(y_bound ==3)
                {
                    temp_phi[4] = phi2[i][2][0];
                }

                temp_phi[1] = phi2[i][j-1][0];
                temp_phi[0] = phi2[i][j-2][0];
            }
            else
            {
                temp_phi[2] = phi2[i][j][0];
                temp_phi[3] = phi2[i][j+1][0];
                temp_phi[4] = phi2[i][j+2][0];
                temp_phi[1] = phi2[i][j-1][0];
                temp_phi[0] = phi2[i][j-2][0];
            }

            dybarplus[i][j] = dplus(temp_phi, 2) -0.5*minmod(dplus(temp_phi,2) - dplus(temp_phi,1) , dplus(temp_phi,3) - dplus(temp_phi,2));
            dybarminus[i][j] = dminus(temp_phi, 2) +0.5*minmod(dplus(temp_phi,2) - dplus(temp_phi,1) , dplus(temp_phi,1) - dplus(temp_phi,0));

           if(signof(phi[i][j][0])*dplus(temp_phi,2) < 0.0  && signof(phi[i][j][0])*dminus(temp_phi,2) < -signof(phi[i][j][0])*dplus(temp_phi,2))
           {
               dybar[i][j] = dybarplus[i][j];
           }
           else if(signof(phi[i][j][0])*dminus(temp_phi,2) > 0.0  && signof(phi[i][j][0])*dplus(temp_phi,2) > -signof(phi[i][j][0])*dminus(temp_phi,2))
           {
               dybar[i][j] = dybarminus[i][j];
           }
           else
           {
               dybar[i][j] = 0.5*(dybarplus[i][j] + dybarminus[i][j]);
           }
        }

    }


    double eps = epsilon*max(xlen/(xelem-2), ylen/(yelem-2));
    #pragma omp parallel for schedule(dynamic)
    for(int i=1; i < xelem-1; i++)
    {
        for(int j=1; j < yelem-1; j++)
        {
            double sign_phi;
            if(phi[i][j][0] >= eps)
            {
                sign_phi = 1.0;
            }
            else if(phi[i][j][0] <= -eps)
            {
                sign_phi =-1.0;
            }
            else
            {
                sign_phi = (phi[i][j][0]/ eps) - (1.0/PI)*sin(PI*phi[i][j][0]/eps);
            }

            rhs[i][j] = sign_phi*(1-sqrt(pow(dxbar[i][j]/area[i][j][1][1],2.0) + pow(dybar[i][j]/area[i][j][0][0],2.0)));
        }
    }
}

#endif /* RHS_BUB_REDIST_H */

