/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   read_write.h
 * Author: nsaini3
 *
 * Created on November 11, 2016, 5:13 PM
 */

#ifndef READ_WRITE_H
#define READ_WRITE_H

void read(elemsclr &sclr)
{
    if(startstep != 0)
    {
        ifstream startu;
        ifstream startp;
        ifstream startv;
        ifstream startphi;
        ifstream startrho;
        ifstream startmu;
        startu.open(getexepath()+"../../../../laststep/u.dat");
        startv.open(getexepath()+"../../../../laststep/v.dat");
        startp.open(getexepath()+"../../../../laststep/p.dat");
        startphi.open(getexepath()+"../../../../laststep/phi.dat");
        startrho.open(getexepath()+"../../../../laststep/rho.dat");
        startmu.open(getexepath()+"../../../../laststep/mu.dat");
        for(int j=0; j<yelem; j++)
        {
            for(int i=0; i<xelem; i++)
            {
                startu>>sclr.u[i][j][0];
                startv>>sclr.v[i][j][0];
                startp>>sclr.p[i][j][0];
                startphi>>sclr.phi[i][j][0];
                startrho>>sclr.rho[i][j][0];
                startmu>>sclr.mu[i][j][0];
            }
        }
        startu.close();
        startv.close();
        startp.close();
        startrho.close();
        startphi.close();
        startmu.close();

    }
}

void write(elemsclr sclr)
{
    ofstream startu;
    ofstream startp;
    ofstream startv;
    ofstream startphi;
    ofstream startrho;
    ofstream startmu;
    startu.open(getexepath()+"../../../../laststep/u.dat",ios::trunc);
    startv.open(getexepath()+"../../../../laststep/v.dat",ios::trunc);
    startp.open(getexepath()+"../../../../laststep/p.dat",ios::trunc);
    startphi.open(getexepath()+"../../../../laststep/phi.dat",ios::trunc);
    startrho.open(getexepath()+"../../../../laststep/rho.dat",ios::trunc);
    startmu.open(getexepath()+"../../../../laststep/mu.dat",ios::trunc);
    for(int j=0; j<yelem; j++)
    {
        for(int i=0; i<xelem; i++)
        {
            startu<<sclr.u[i][j][0]<<endl;
            startv<<sclr.v[i][j][0]<<endl;
            startp<<sclr.p[i][j][0]<<endl;
            startphi<<sclr.phi[i][j][0]<<endl;
            startrho<<sclr.rho[i][j][0]<<endl;
            startmu<<sclr.mu[i][j][0]<<endl;
        }
    }
    startu.close();
    startv.close();
    startp.close();
    startrho.close();
    startphi.close();
    startmu.close();
}

#endif /* READ_WRITE_H */

