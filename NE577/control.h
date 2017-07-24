/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   control.h
 * Author: nsaini3
 *
 * Created on September 26, 2016, 4:35 PM
 */

#ifndef CONTROL_H
#define CONTROL_H

void control()
{
    ifstream controlfile;
    controlfile.open("control.txt");
    string line;
    
    while(getline(controlfile,line))
    {
        string junk;
        istringstream ss(line);
        ss>>junk;
        string junk2;
        
        if(junk.compare("xlen") == 0)
        {
            ss>>xlen;
        }
        else if(junk.compare("ylen") == 0)
        {
            ss>>ylen;
        }
        else if(junk.compare("zlen") == 0)
        {
            ss>>zlen;
        }
        else if(junk.compare("xelem") == 0)
        {
            ss>>xelem;
        }
        else if(junk.compare("yelem") == 0)
        {
            ss>>yelem;
        }
        else if(junk.compare("zelem") == 0)
        {
            ss>>zelem;
        }
        else if(junk.compare("Max_Iterations") == 0)
        {
            ss>>itermax;
        }
        else if(junk.compare("Kinematic_viscosity") == 0)
        {
            ss>>nu;
        }
        else if(junk.compare("Tolerance_x") == 0)
        {
            ss>>tol;
        }
        else if(junk.compare("Bubble_radius") == 0)
        {
            ss>>rb_in;
        }
        else if(junk.compare("x_pos_of_bubble") == 0)
        {
            ss>>xb_in;
        }
        else if(junk.compare("y_pos_of_bubble") == 0)
        {
            ss>>yb_in;
        }
        else if(junk.compare("advect_steps") == 0)
        {
            ss>>advect_steps;
        }
        else if(junk.compare("advect_deltat") == 0)
        {
            ss>>advect_deltat;
        }
        else if(junk.compare("Solution_read") == 0)
        {
            ss>>solnread;
        }
        else if(junk.compare("bub_advect_scheme") == 0)
        {
            ss>>junk2;
            if(junk2.compare("book") == 0)
            {
                bub_conv_scheme=1;
            }
            else if(junk2.compare("HJ-WENO") == 0)
            {
                bub_conv_scheme=2;
            }
            
        }
        else if(junk.compare("Liquid_density") == 0)
        {
            ss>>rhof;
        }
        else if(junk.compare("Gas_density") == 0)
        {
            ss>>rhog;
        }
        else if(junk.compare("Liquid_viscosity") == 0)
        {
            ss>>muf;
        }
        else if(junk.compare("Gas_viscosity") == 0)
        {
            ss>>mug;
        }
        else if(junk.compare("Epsilon") == 0)
        {
            ss>>epsilon;
        }
        else if(junk.compare("Surface_tension_coefficient") == 0)
        {
            ss>>sf_coeff;
        }
        else if(junk.compare("GS_relaxation_factor") == 0)
        {
            ss>>relax;
        }
        else if(junk.compare("Tolerance_p") == 0)
        {
            ss>>ptol;
        }
        else if(junk.compare("Re_distance_timestep") == 0)
        {
            ss>>re_time;
        }
        else if(junk.compare("Re_distance_loops") == 0)
        {
            ss>>re_loops;
        }
        else if(junk.compare("print_gap") == 0)
        {
            ss>>print_gap;
        }
        
        else if(junk.compare("Start_step") == 0)
        {
            ss>>startstep;
        }
        
        else if(junk.compare("Surface_tension") == 0)
        {
            ss>>sf_toggle;
        }
        else if(junk.compare("Solve_flow") == 0)
        {
            ss>>flow_solve;
        }
        else if(junk.compare("Variable_density_pressure_solver") == 0)
        {
            ss>>p_solver;
        }
        else if(junk.compare("Advect_bubble") == 0)
        {
            ss>>advect_solve;
        }
        else if(junk.compare("Void_fraction_control") == 0)
        {
            ss>>vf_control;
        }
        else if(junk.compare("Re-distance_method") == 0)
        {
            ss>>redist_method;
        }
        else if(junk.compare("x-boundary") == 0)
        {
            ss>>junk2;
            if(junk2.compare("no-slip") == 0)
            {
                x_bound=1;
            }
            else if(junk2.compare("slip") == 0)
            {
                x_bound=2;
            }
            else if(junk2.compare("periodic") == 0)
            {
                x_bound=3;
            }
            
        }
        
        else if(junk.compare("y-boundary") == 0)
        {
            ss>>junk2;
            if(junk2.compare("no-slip") == 0)
            {
                y_bound=1;
            }
            else if(junk2.compare("slip") == 0)
            {
                y_bound=2;
            }
            else if(junk2.compare("periodic") == 0)
            {
                y_bound=3;
            }
            
        }
        
        else if(junk.compare("Solver_type") == 0)
        {
            ss>>junk2;
            if(junk2.compare("steady-state") == 0)
            {
                sol_type=0;
            }
            else if(junk2.compare("transient") == 0)
            {
                sol_type=1;
            }
            
        }
        
        else if(junk.compare("Time_control") == 0)
        {
            ss>>junk2;
            if(junk2.compare("CFL-based") == 0)
            {
                time_control=1;
            }
            else if(junk2.compare("constant_time") == 0)
            {
                time_control=2;
            }
            
        }
        
        else if(junk.compare("max_CFL") == 0)
        {
            ss>>max_cfl;
        }
        else if(junk.compare("gx") == 0)
        {
            ss>>gx;
        }
        else if(junk.compare("gy") == 0)
        {
            ss>>gy;
        }
        
        else if(junk.compare("Case") == 0)
        {
            ss>>junk2;
            if(junk2.compare("vortex") == 0)
            {
                case_tog=1;
            }
            else if(junk2.compare("bubble_rise") == 0)
            {
                case_tog=2;
            }
            else if(junk2.compare("bubble_break") == 0)
            {
                case_tog=3;
            }
            else if(junk2.compare("zalesak") == 0)
            {
                case_tog=4;
            }
            else if(junk2.compare("homework") == 0)
            {
                case_tog=5;
            }
        }
    }
    controlfile.close();
}

#endif /* CONTROL_H */

