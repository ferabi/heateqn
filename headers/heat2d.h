#pragma once
#include "vector.h"
#include "matrix.h"
#include "vector_operations.h"
#include "conjugategradient.h"
#include <math.h>
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

class Heat2D
{
    private :
        //diffusion coefficient = alpha
        double alpha; 
        //number of points per dimension
        unsigned int m;
        //timestep
        double dt;
        //iteration matrix
        Matrix<double> M;
    
    public:
        //constructor takes three arguments
        Heat2D(const double in_alpha,const unsigned int in_m, const double in_dt);

        //exact solution at time t
        Vector<double> exact(double in_time) const;

        //numerical solution at time t_end
        Vector<double> solve(double time_end) const;
};



