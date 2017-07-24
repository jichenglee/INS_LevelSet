/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   bub_advect.h
 * Author: nsaini3
 *
 * Created on October 18, 2016, 7:42 PM
 */

#ifndef BUB_ADVECT_H
#define BUB_ADVECT_H




void bub_advect(elemsclr &sclr, int iter, double deltat)
{


    ///Interpolate velocity at cell edges to cell centers
    vector< vector<vector<double> > > ucen(xelem, vector<vector<double> > (yelem,vector<double> (zelem,0.0))); //Stored at cell centers
    vector< vector<vector<double> > > vcen(xelem, vector<vector<double> > (yelem,vector<double> (zelem,0.0)));
    for(int i=1; i < xelem-1; i++)
    {
        for(int j=1; j < yelem-1; j++)
        {
            ucen[i][j][0]=0.5*(sclr.u[i][j][0] + sclr.u[i-1][j][0]);
            vcen[i][j][0]=0.5*(sclr.v[i][j][0] + sclr.v[i][j-1][0]);
        }
    }
    cell_center_vel_BC(ucen,vcen);



//    Main bubble iteration loop
    if(iter == 0)
    {
        output_vtk(sclr,iter);
    }




    vector< vector<double> > rhsx(xelem, vector<double> (yelem,0.0));
    vector< vector<double> > rhsy(xelem, vector<double> (yelem,0.0));
    //Calculate the fluxes
    rhs_bub(rhsx, rhsy, ucen, vcen, sclr.phi);

    vector< vector<vector<double> > > phistar(xelem, vector<vector<double> > (yelem,vector<double> (zelem,0.0)));
    #pragma omp parallel for schedule(dynamic)
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            phistar[i][j][0] = sclr.phi[i][j][0] + deltat * (rhsx[i][j] + rhsy[i][j]);
        }
    }

    level_setBC(phistar);

    vector< vector<double> > rhstarx(xelem, vector<double> (yelem,0.0));
    vector< vector<double> > rhstary(xelem, vector<double> (yelem,0.0));
    //Calculate the star fluxes
    rhs_bub(rhstarx, rhstary, ucen, vcen, phistar);
    #pragma omp parallel for schedule(dynamic)
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            //cout<<phi[i][j][0];
            sclr.phi[i][j][0] = sclr.phi[i][j][0] + 0.5*deltat *(rhsx[i][j] + rhstarx[i][j] + rhsy[i][j] + rhstary[i][j]);
            //cout<<" "<<phi[i][j][0]<<endl;
        }
    }

    level_setBC(sclr.phi);


}

#endif /* BUB_ADVECT_H */

