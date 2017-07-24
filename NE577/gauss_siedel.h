/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   gauss_siedel.h
 * Author: nsaini3
 *
 * Created on November 3, 2016, 7:18 PM
 */

#ifndef GAUSS_SIEDEL_H
#define GAUSS_SIEDEL_H

void gs_solver(vector< vector< vector<double> > > &a, vector< vector<double> > &b, vector< vector< vector<double> > > &p)
{
    vector< vector< vector<double> > > tempp(xelem, (vector< vector<double> >(yelem, vector<double>(zelem,0.0))));
    vector< vector<double> > delp(xelem, vector<double> (yelem,0.0));
    vector< vector<double> > d(xelem, vector<double> (yelem,0.0));

    double ires;

    for(int j=0; j<yelem; j++)
    {
        for(int i=0; i<xelem; i++)
        {
            tempp[i][j][0] = p[i][j][0];
            //cout<<b[i][j]<<" ";
        }
        //cout<<endl;
    }

    //exit(0);
    int iter;
    for(iter=0; iter< 10000; iter++)
    {
        //cout<<tempp[1][1][0]<<endl;
        #pragma omp parallel for schedule(dynamic)
        for(int i=1; i<xelem-1; i++)
        {
            for(int j=1; j<yelem-1; j++)
            {
                double res = b[i][j] - (a[i][j][0]*tempp[i-1][j][0] + a[i][j][1]*tempp[i][j-1][0] + a[i][j][3]*tempp[i][j+1][0] + a[i][j][4]*tempp[i+1][j][0]);
                delp[i][j] = res/a[i][j][2];

            }
            /*if(i==1)
            {
                cout<<b[1][1]<<" "<<delp[1][1]<<endl;
            }*/
        }

        for(int i=1; i<xelem-1; i++)
        {
            for(int j=1; j<yelem-1; j++)
            {
                tempp[i][j][0] = delp[i][j];
            }
        }
        //cout<<tempp[1][1][0]<<endl;
        pressureBC(tempp);

        double resnorm=0.0;
        for(int i=1; i<xelem-1; i++)
        {
            for(int j=1; j<yelem-1; j++)
            {
                double res = b[i][j] - (a[i][j][0]*tempp[i-1][j][0] + a[i][j][1]*tempp[i][j-1][0] + a[i][j][3]*tempp[i][j+1][0] + a[i][j][4]*tempp[i+1][j][0] + a[i][j][2]*tempp[i][j][0]);

                resnorm = resnorm + pow(res,2.0)*vol[i][j];
            }
        }

        if(iter == 0)
        {
            ires = resnorm;
            //cout<<"Initial residual "<<ires<<endl;
            //exit(0);
        }
        else
        {
            //cout<<"Pressure Step: "<<iter<<"residual: "<<resnorm/ires<<endl;
            if(resnorm / ires < ptol)
            {
                cout<<"Pressure converged in "<<iter<<" "<<endl;
                break;
            }
        }
    }

    pressureBC(tempp);
    cout<<"Pressure iterations "<<iter<<endl;
    for(int i=0; i<xelem; i++)
    {
        for(int j=0; j<yelem; j++)
        {
            p[i][j][0] = tempp[i][j][0];
        }
    }

}

#endif /* GAUSS_SIEDEL_H */

