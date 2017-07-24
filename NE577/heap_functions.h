/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   heap_functions.h
 * Author: User
 *
 * Created on November 25, 2016, 10:39 PM
 */

#ifndef HEAP_FUNCTIONS_H
#define HEAP_FUNCTIONS_H

void swap(int min_index, int i, vector <double> &heap1, vector< vector<int> > &index)
{
    int temp[2];
    double tempvalue;
    
    temp[0] = index[min_index][0];
    temp[1] = index[min_index][1];
    
    index[min_index][0] = index[i][0];
    index[min_index][1] = index[i][1];
    
    index[i][0] = temp[0];
    index[i][1] = temp[1];
    
    tempvalue = heap1[min_index];
    heap1[min_index] = heap1[i];
    heap1[i] = tempvalue;
}

void min_heapify(int i, vector <double> &heap1, vector< vector<int> > &index)
{
    i++;
    int l,r,n,min_index;
    
    n = heap1.size();
    l = 2*i;
    r = 2*i+1;
    
    if(l <= n && heap1[l-1] < heap1[i-1])
    {
        min_index = l-1;
    }
    else
    {
        min_index = i-1;
    }
    
    if(r <= n && heap1[r-1] < heap1[min_index])
    {
        min_index = r-1;
    }
    
    if(min_index != i-1)
    {
        swap(min_index, i-1, heap1, index);
        min_heapify(min_index, heap1, index);
    }
}

void initialize_heap1(vector <double> &heap1, vector< vector<int> > &index)
{
    int n,i;
    
    n = heap1.size();
    i=n/2;
    
    while(i >= 1)
    {
        min_heapify(i-1, heap1, index);
        i--;
    }
}

void delete_from_heap(vector <double> &heap1, vector< vector<int> > &index)
{
    int n, i, node_del;
    
    n=heap1.size();
    swap(0,n-1,heap1,index);
    
    heap1.resize(n-1);
    index[n-1].resize(0);
    index.resize(n-1);
    
    i=0;
    min_heapify(i, heap1, index);
}


void add_to_heap(double phivalue, int nodex, int nodey, vector <double> &heap1, vector< vector<int> > &index)
{
    int n,i;
    
    n=heap1.size();
    
    heap1.resize(n+1);
    index.resize(n+1);
    index[n].resize(2);
    
    heap1[n] = phivalue;
    index[n][0] = nodex;
    index[n][1] = nodey;
    
    i = (n+1)/2;
    
    while(i >= 1)
    {
        min_heapify(i-1, heap1,index);
        i--;
    }
    
}

#endif /* HEAP_FUNCTIONS_H */

