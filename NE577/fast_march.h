/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   fast_march.h
 * Author: User
 *
 * Created on November 25, 2016, 6:53 PM
 */

#ifndef FAST_MARCH_H
#define FAST_MARCH_H

#include "heap_functions.h"

void interface_nodes(vector< vector< vector<double> > > &phi, vector< vector<int> > &int_cells)
{
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            if(phi[i][j][0]* phi[i+1][j][0] < 0.0)
            {
                    
                /*Check if the cell i,j is already in list*/
                if(int_cells.size() > 1)
                {
                    bool inlist = false;
                    for(int k=int_cells.size()-1; k > 0; k--)
                    {
                        if(i == int_cells[k][0] && j == int_cells[k][1])
                        {
                            inlist = true;
                        }
                    }
                    if(inlist == false)
                    {
                        int_cells.resize(int_cells.size() + 1);
                        int_cells[int_cells.size()-1].resize(2);
                        
                        int_cells[int_cells.size()-1][0] = i;
                        int_cells[int_cells.size()-1][1] = j;
                    }
                }
                else //If list is empty just add the cell to list
                {
                    int_cells.resize(int_cells.size() + 1);
                    int_cells[int_cells.size()-1].resize(2);

                    int_cells[int_cells.size()-1][0] = i;
                    int_cells[int_cells.size()-1][1] = j;
                }
                
                /*Check if the cell i+1,j is already in list*/
                bool inlist = false;
                for(int k=int_cells.size()-1; k > 0; k--)
                {
                    if(i+1 == int_cells[k][0] && j == int_cells[k][1])
                    {
                        inlist = true;
                    }
                }
                if(inlist == false)
                {
                    int_cells.resize(int_cells.size() + 1);
                    int_cells[int_cells.size()-1].resize(2);

                    int_cells[int_cells.size()-1][0] = i+1;
                    int_cells[int_cells.size()-1][1] = j;
                }
            }
            
            
            /*Now for the j direction*/
            if(phi[i][j][0]* phi[i][j+1][0] < 0.0)
            {
                    
                /*Check if the cell i,j is already in list*/
                if(int_cells.size() > 1)
                {
                    bool inlist = false;
                    for(int k=int_cells.size()-1; k > 0; k--)
                    {
                        if(i == int_cells[k][0] && j == int_cells[k][1])
                        {
                            inlist = true;
                        }
                    }
                    if(inlist == false)
                    {
                        int_cells.resize(int_cells.size() + 1);
                        int_cells[int_cells.size()-1].resize(2);
                        
                        int_cells[int_cells.size()-1][0] = i;
                        int_cells[int_cells.size()-1][1] = j;
                    }
                }
                else //If list is empty just add the cell to list
                {
                    int_cells.resize(int_cells.size() + 1);
                    int_cells[int_cells.size()-1].resize(2);

                    int_cells[int_cells.size()-1][0] = i;
                    int_cells[int_cells.size()-1][1] = j;
                }
                
                /*Check if the cell i+1,j is already in list*/
                bool inlist = false;
                for(int k=int_cells.size()-1; k > 0; k--)
                {
                    if(i == int_cells[k][0] && j+1 == int_cells[k][1])
                    {
                        inlist = true;
                    }
                }
                if(inlist == false)
                {
                    int_cells.resize(int_cells.size() + 1);
                    int_cells[int_cells.size()-1].resize(2);

                    int_cells[int_cells.size()-1][0] = i;
                    int_cells[int_cells.size()-1][1] = j+1;
                }
            }

        }
    }
}

void find_ne(vector<int> &ne_i,  vector <int> &ne_j, int i, int j)
{
    /*Left and right neighbours*/
    if(i == xelem-2 && x_bound == 3)
    {
        ne_i[0] = 1;
        ne_j[0] = j;
        
        ne_i[2] = i-1;
        ne_j[2] = j;
    }
    else if(i == 1 && x_bound == 3)
    {
        ne_i[0] = i+1;
        ne_j[0] = j;
        
        ne_i[2] = xelem-2;
        ne_j[2] = j;
    }
    else if(i == xelem-2 && x_bound !=3)
    {
        ne_i[0] = 10000;
        ne_j[0] = 10000;
        
        ne_i[2] = i-1;
        ne_j[2] = j;
    }
    else if(i == 1 && x_bound != 3)
    {
        ne_i[0] = i+1;
        ne_j[0] = j;
        
        ne_i[2] = 10000;
        ne_j[2] = 10000;
    }
    else
    {
        ne_i[0] = i+1;
        ne_j[0] = j;
        
        ne_i[2] = i-1;
        ne_j[2] = j;
    }
    
    
    /*Up and down neighbours*/
    if(j == yelem-2 && y_bound == 3)
    {
        ne_i[1] = i;
        ne_j[1] = 1;
        
        ne_i[3] = i;
        ne_j[3] = j-1;
    }
    else if(j == 1 && y_bound == 3)
    {
        ne_i[1] = i;
        ne_j[1] = j+1;
        
        ne_i[3] = i;
        ne_j[3] = yelem-2;
    }
    else if(j == yelem-2 && y_bound !=3)
    {
        ne_i[1] = 10000;
        ne_j[1] = 10000;
        
        ne_i[3] = i;
        ne_j[3] = j-1;
    }
    else if(j == 1 && y_bound != 3)
    {
        ne_i[1] = i;
        ne_j[1] = j+1;
        
        ne_i[3] = 10000;
        ne_j[3] = 10000;
    }
    else
    {
        ne_i[1] = i;
        ne_j[1] = j+1;
        
        ne_i[3] = i;
        ne_j[3] = j-1;
    }
    
}


void compute_phi(int i, int j, vector< vector<double> > &phi2, vector< vector<int> > &tag, int solve_flag)
{
    vector<int> ne_i(4,0);
    vector<int> ne_j(4,0);
    
    find_ne(ne_i, ne_j, i, j);
    
    double xphi, yphi;
    double xdsq, ydsq;
    double xflag=1.0;
    double yflag=1.0;
    
    double ne_phi[4];
    
    for(int n=0; n<4; n++)
    {
        if(ne_i[n] != 10000)
        {
            if(tag[ne_i[n]][ne_j[n]] == solve_flag)
            {
                ne_phi[n] = phi2[ne_i[n]][ne_j[n]];
            }
            else
            {
                ne_phi[n] = 10000.0;
            }
        }
        else
        {
            ne_phi[n] = 10000.0;
        }
    }
    
    double multiplier;
    if(solve_flag == 1)
    {
        multiplier = 1.0;
    }
    else if(solve_flag == -1)
    {
        multiplier = -1.0;
    }
    if(ne_phi[0] > 9999 && ne_phi[2] > 9999)
    {
        xphi = 0.0;
        xdsq = 1.0;
        xflag = 0.0;
    }
    else
    {
        xphi = multiplier * min(fabs(ne_phi[0]), fabs(ne_phi[2]));
        xdsq = pow(area[1][1][1][1], 2.0);
    }
    
    if(ne_phi[1] > 9999 && ne_phi[3] > 9999)
    {
        yphi = 0.0;
        ydsq = 1.0;
        yflag = 0.0;
    }
    else
    {
        yphi = multiplier * min(fabs(ne_phi[1]), fabs(ne_phi[3]));
        ydsq = pow(area[1][1][0][0], 2.0);
    }
    
    
    double anew = ydsq*xflag + xdsq*yflag;
    double bnew = -2.0*(xphi*ydsq + yphi*xdsq);
    double cnew = xphi*xphi*ydsq + yphi*yphi*xdsq - xdsq*ydsq;
    
    double disc = bnew*bnew - 4.0*anew*cnew;
    
    if(disc < 0.0)
    {
        cout<<"Solution to neighbour not found"<<endl;
        cout<<i<<" "<<j<<" "<<xelem<<" "<<yelem<<endl;
        cout<<xphi<<" "<<xdsq<<" "<<xflag<<endl;
        cout<<yphi<<" "<<ydsq<<" "<<yflag<<endl;
        exit(0);
    }
    
    else if(disc >= 0.0 && solve_flag == 1)
    {
        phi2[i][j] = (-bnew + sqrt(disc))/(2.0*anew);
    }
    
    else if(disc >=0.0 && solve_flag == -1)
    {
        phi2[i][j] = (-bnew - sqrt(disc))/(2.0*anew);
    }
    
}


void fast_march(elemsclr &sclr)
{
    vector< vector<double> > phi2(xelem, vector<double> (yelem,0.0));
    vector< vector<int> > tag(xelem, vector<int> (yelem,0));
    
    vector< vector<int> > int_cells;
    
    /****Calculate interface nodes*****/
    interface_nodes(sclr.phi, int_cells);
    //cout<<"No of cells on interface; "<<int_cells.size()<<endl;
    /*for(int i=0; i<int_cells.size(); i++)
    {
        cout<<int_cells[i][0]<<" "<<int_cells[i][1]<<endl;
    }*/
    
    /*Tag interface nodes as 1*/
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            for(int k=0; k < int_cells.size(); k++)
            {
                if(i == int_cells[k][0] && j == int_cells[k][1] && sclr.phi[i][j][0] > 0.0)
                {
                    tag[i][j] = 1;
                }
                else if(i == int_cells[k][0] && j == int_cells[k][1] && sclr.phi[i][j][0] < 0.0)
                {
                    tag[i][j] = -1;
                }
            }
        }
    }
    
    /*Copy the phi values to phi2*/
    for(int i=0; i<xelem; i++)
    {
        for(int j=0; j<yelem; j++)
        {
            phi2[i][j] = sclr.phi[i][j][0];
        }
    }
    
    int solve_flag = 1;
    
    again:
    
    if(solve_flag == 1)
    {
        cout<<"Solving +ve nodes"<<endl;
    }
    else if(solve_flag == -1)
    {
        cout<<"Solving -ve nodes"<<endl;
    }
    
    /*Tag rest of the nodes as far nodes*/
    int count=0;
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            if(tag[i][j] !=1 && solve_flag == 1 && sclr.phi[i][j][0] > 0.0 ||  tag[i][j] !=-1 && solve_flag == -1 && sclr.phi[i][j][0] < 0.0)
            {
                tag[i][j] = 3;
                count++;
            }
        }
    }
    //cout<<"total far nodes "<<count<<endl;
    
    /*Now compute the phi values of interface nodes*/
    for(int k=0; k < int_cells.size(); k++)
    {
        vector<int> ne_i(4,0);
        vector<int> ne_j(4,0);
        find_ne(ne_i, ne_j, int_cells[k][0], int_cells[k][1]);
        
        for(int n=0; n<4; n++)
        {
            if(ne_i[n] != 10000 && ne_j[n] != 10000)
            {
                if(solve_flag == 1 && sclr.phi[ne_i[n]][ne_j[n]][0] > 0.0 || solve_flag == -1 && sclr.phi[ne_i[n]][ne_j[n]][0] < 0.0)
                {
                    if(tag[ne_i[n]][ne_j[n]] != 1 && tag[ne_i[n]][ne_j[n]] != -1 && tag[ne_i[n]][ne_j[n]] != 2)
                    {
                        compute_phi(ne_i[n], ne_j[n], phi2, tag, solve_flag);
                        tag[ne_i[n]][ne_j[n]] = 2;
                    }
                }
            }
        }
    }
    
    
    vector<double> heap1;
    vector< vector<int> > index;
    
    int tot_tag=0;
    for(int i=1; i < xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            if(tag[i][j] == 2)
            {
                tot_tag++;
            }
        }
    }
    
    heap1.resize(tot_tag);
    index.resize(tot_tag);
    for(int i=0; i < tot_tag; i++)
    {
        index[i].resize(2);
    }
    
    int heap_index=0;
    for(int i=1; i < xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            if(tag[i][j] == 2)
            {
                heap1[heap_index] = abs(phi2[i][j]);
                index[heap_index][0] = i;
                index[heap_index][1] = j;
                heap_index++;
            }
        }
    }
    
    initialize_heap1(heap1, index);
    
    
    
    int band = 1;
    int ncount = 0;
    
    
    
    while(heap1.size() > 0)
    {
        ncount++;
        
        //Tag top of heap as 1 and delete
        tag[index[0][0]][index[0][1]] = solve_flag;
        
        
        
        int new_node[2];
        new_node[0] = index[0][0];
        new_node[1] = index[0][1];
        
        delete_from_heap(heap1,index);
        
        //Find neighbours of new node, tag as band and find new phi values
        vector<int> ne_i(4);
        vector<int> ne_j(4);
        find_ne(ne_i, ne_j, new_node[0], new_node[1]);
        
        for(int n=0; n<4; n++)
        {
            if(ne_i[n] != 10000)
            {
                if(tag[ne_i[n]][ne_j[n]] != 1 && tag[ne_i[n]][ne_j[n]] != -1)
                {
                    compute_phi(ne_i[n], ne_j[n], phi2, tag, solve_flag);
                    
                    if(tag[ne_i[n]][ne_j[n]] == 2)
                    {
                        for(int k=0; k < heap1.size(); k++)
                        {
                            if(index[k][0] == ne_i[n] && index[k][1] == ne_j[n])
                            {
                                heap1[k] = fabs(phi2[ne_i[n]][ne_j[n]]);
                                min_heapify(0, heap1, index);
                                break;
                            }
                        }
                    }
                    else
                    {
                        tag[ne_i[n]][ne_j[n]] = 2;
                        int nodex = ne_i[n];
                        int nodey = ne_j[n];
                        
                        double phivalue = fabs(phi2[ne_i[n]][ne_j[n]]);
                        add_to_heap(phivalue, nodex, nodey, heap1, index);
                    }
                    
                }
            }
        }
        
        //cout<<heap1.size()<<endl;
        //cout<<index.size()<<endl;
    }
    
    //cout<<"Total nodes processed "<<ncount<<endl;
    
    heap1.resize(0);
    index.resize(0);
    
    if(solve_flag == 1)
    {
        solve_flag = -1;
        goto again;
    }
    
    //cout<<heap1.size()<<endl;
    //cout<<index.size()<<endl;
    
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            sclr.phi[i][j][0] = phi2[i][j];
        }
    }
    level_setBC(sclr.phi);
    
    //cout<<"here"<<endl;
    //exit(0);
            
}

#endif /* FAST_MARCH_H */

