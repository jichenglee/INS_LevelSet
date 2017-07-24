/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   output.h
 * Author: User
 *
 * Created on September 27, 2016, 1:46 AM
 */

#ifndef OUTPUT_H
#define OUTPUT_H

void output(elemsclr sclr,int iter)
{

    double sum=0;
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            double exact=2400*yc[i][j]*(0.01-yc[i][j]);
            double err=(sclr.u[i][j][0]-exact)/exact;
            sum=sum+pow(err,2.0)*vol[i][j];
        }
    }
    double error=sqrt(sum);

    ofstream out;
    ofstream compare;
    out.open("output.dat",ios::trunc);
    out<<"variables ="<<" x,"<<" y,"<<" u,"<<" v,"<<" p,"<<" phi"<<endl;//" p"<<endl;
    out<<"zone "<<"i="<<xnode-2<<", "<<"j="<<ynode-2<<", "<<"f=point"<<endl;
    for (int j=1;j<ynode-1;j++)
    {
        for (int i=1;i<xnode-1;i++)
        {
            double unode=0.5*(sclr.u[i-1][j][0]+sclr.u[i-1][j-1][0]);
            double vnode=0.5*(sclr.v[i-1][j-1][0] + sclr.v[i][j-1][0]);
            //double vnode=0.25*(v[i][j]+v[i+1][j]+v[i+1][j+1]+v[i][j+1]);
            double pnode=0.25*(sclr.p[i][j][0]+sclr.p[i-1][j][0]+sclr.p[i-1][j-1][0]+sclr.p[i][j-1][0]);
            double phinode=0.25*(sclr.phi[i][j][0]+sclr.phi[i-1][j][0]+sclr.phi[i-1][j-1][0]+sclr.phi[i][j-1][0]);

            out<<x[i][j]<<" "<<y[i][j]<<" "<<unode<<" "<<vnode<<" "<<pnode<<" "<<phinode<<endl;//" "<<pnode<<" "<<endl;


        }
    }
    out.close();
    //ofstream output2;
    //output2.open("info.txt",ios::trunc);
    if(solnread == 0)
    {
        ofstream soln;
        soln.open("steady_state_soln.txt",ios::trunc);
        for(int j=0; j<yelem; j++)
        {
            for(int i=0; i<xelem ; i++)
            {
                soln<<sclr.u[i][j][0]<<" "<<sclr.v[i][j][0]<<" "<<sclr.p[i][j][0]<<endl;
            }
        }
        soln.close();
    }

    cout<<"Total Iterations for Convergence = "<<iter<<endl;

    cout<<"2-Norm of Error  = "<<error<<endl;
    cout<<"Overall Order = "<<-log10(fabs(error))<<endl;
    //output2.close();
}



void output_vtk(elemsclr &sclr,int iter)//, vector< vector< vector<double> > > &stx, vector< vector< vector<double> > > &sty)
{
    ofstream out;
    out.open(getexepath()+"../../../../output/out_00"+inttostr(iter)+".vts",ios::trunc);

    //Write Headers
    out<<"# vtk DataFile Version 3.0"<<endl;
    out<<"vtk output"<<endl;
    out<<"ASCII"<<endl;
    out<<"DATASET STRUCTURED_GRID"<<endl;
    out<<"DIMENSIONS "<<xnode-2<<" "<<ynode-2<<" "<<1<<endl;
    out<<"POINTS "<<(xnode-2)*(ynode-2)<<" double"<<endl;

    for(int j=1; j<ynode-1; j++)
    {
        for(int i=1; i<xnode-1; i++)
        {
            out<<x[i][j]<<" "<<y[i][j]<<" "<<0<<endl;
        }
    }

    out<<"POINT_DATA "<<(xnode-2)*(ynode-2)<<endl;
    out<<"SCALARS u double"<<endl;
    out<<"LOOKUP_TABLE default"<<endl;

    for(int j=1; j<ynode-1; j++)
    {
        for(int i=1; i<xnode-1; i++)
        {

            double unode=0.5*(sclr.u[i-1][j][0]+sclr.u[i-1][j-1][0]);
            out<<unode<<endl;
        }
    }

    out<<"SCALARS v double"<<endl;
    out<<"LOOKUP_TABLE default"<<endl;

    for(int j=1; j<ynode-1; j++)
    {
        for(int i=1; i<xnode-1; i++)
        {

            double vnode=0.5*(sclr.v[i-1][j][0]+sclr.v[i-1][j-1][0]);
            out<<vnode<<endl;
        }
    }

    out<<"SCALARS p double"<<endl;
    out<<"LOOKUP_TABLE default"<<endl;

    for(int j=1; j<ynode-1; j++)
    {
        for(int i=1; i<xnode-1; i++)
        {
            double pnode=0.25*(sclr.p[i][j][0]+sclr.p[i-1][j][0]+sclr.p[i-1][j-1][0]+sclr.p[i][j-1][0]);
            out<<pnode<<endl;
        }
    }


    out<<"SCALARS phi double"<<endl;
    out<<"LOOKUP_TABLE default"<<endl;

    for(int j=1; j<ynode-1; j++)
    {
        for(int i=1; i<xnode-1; i++)
        {
            double phinode=0.25*(sclr.phi[i][j][0]+sclr.phi[i-1][j][0]+sclr.phi[i-1][j-1][0]+sclr.phi[i][j-1][0]);
            out<<phinode<<endl;
        }
    }

    out<<"SCALARS rho double"<<endl;
    out<<"LOOKUP_TABLE default"<<endl;

    for(int j=1; j<ynode-1; j++)
    {
        for(int i=1; i<xnode-1; i++)
        {
            double rhonode=0.25*(sclr.rho[i][j][0]+sclr.rho[i-1][j][0]+sclr.rho[i-1][j-1][0]+sclr.rho[i][j-1][0]);
            out<<rhonode<<endl;
        }
    }

    out<<"SCALARS mu double"<<endl;
    out<<"LOOKUP_TABLE default"<<endl;

    for(int j=1; j<ynode-1; j++)
    {
        for(int i=1; i<xnode-1; i++)
        {
            double munode=0.25*(sclr.mu[i][j][0]+sclr.mu[i-1][j][0]+sclr.mu[i-1][j-1][0]+sclr.mu[i][j-1][0]);
            out<<munode<<endl;
        }
    }

   /* out<<"SCALARS st_forcex double"<<endl;
    out<<"LOOKUP_TABLE default"<<endl;

    for(int j=1; j<ynode-1; j++)
    {
        for(int i=1; i<xnode-1; i++)
        {
            double stx_node=0.25*(stx[i][j][0]+stx[i-1][j][0]+stx[i-1][j-1][0]+stx[i][j-1][0]);
            out<<stx_node<<endl;
        }
    }

    out<<"SCALARS st_forcey double"<<endl;
    out<<"LOOKUP_TABLE default"<<endl;

    for(int j=1; j<ynode-1; j++)
    {
        for(int i=1; i<xnode-1; i++)
        {
            double sty_node=0.25*(sty[i][j][0]+sty[i-1][j][0]+sty[i-1][j-1][0]+sty[i][j-1][0]);
            out<<sty_node<<endl;
        }
    }*/

    out<<"CELL_DATA "<<(xelem-2)*(yelem-2)<<endl;
    out.close();

}
#endif /* OUTPUT_H */

