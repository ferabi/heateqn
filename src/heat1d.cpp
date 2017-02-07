#include "heat1d.h"

inline double laplace_op(const unsigned int i, const unsigned int j)
{
    if(i == j)
    { 
        return -2.0;
    }
    if (i-j == -1 || i-j == 1)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

inline Vector<double> initial_solvec(const unsigned int points)
{
    Vector<double> result(points);

    double dx = 1.0 /(points + 1.0);
    double x;
    for (auto i = 0; i < points; i++)
    {
        x = dx * (i+1);
        result[i] = sin(M_PI * x);
    }

    return result;
}




Heat1D::Heat1D(const double in_alpha,const unsigned int in_m, const double in_dt):
    alpha(in_alpha), m(in_m), dt(in_dt), M(in_m)
{
    //initialisiing the iteration matrix
    double dx = 1.0 / ( this->m + 1.0 );
    double dxsquared = dx * dx;
    double expr;

    for (auto i = 0; i < m; i++)
    {
        for(auto j = i; j < m; j++)
        {
            double ones = (i==j)?1.0:0.0;
            expr = ones - (alpha * (dt/dxsquared) * laplace_op(i,j));
            
            if (expr != 0)
            {
                M[{i,j}] = expr;
                M[{j,i}] = M[{i,j}];
            }
        }
    }
		
	//output of the matrix built
    //std::cout << M << std::endl;
}

Vector<double> Heat1D::exact(double in_time) const
{
    //making a vector to keep the exact solution in
    Vector<double> exact_sol(m);

    Vector<double> in_sol = initial_solvec(m);
    
    //filling the vector
    for (auto i = 0; i < m; i++)
    {
        auto exponent = - 1 * M_PI * M_PI* alpha * in_time;
        
		exact_sol[i] = exp(exponent) * in_sol[i];
    }
    
    return exact_sol;
}

Vector<double> Heat1D::solve(double time_end) const
{
    //check if the input time_end is a multiple of dt
    if ( remainder(time_end,dt) > 1e-10 && time_end > 0.0)
    {
        std::cout << "Input time is not a multiple of dt" << std::endl;
        
		return Vector<double>(m);
    }
 
    //initial solution vector w_l when l=0 i.e. w_0
    Vector<double> w_l = initial_solvec(m);

    //make the solution vector w_(l+1)
    Vector<double> w_l1 = w_l;

    for(auto i = 0; i < time_end/dt ; i++)
    {
        auto steps = cg(M, w_l, w_l1, 1e-17, 100000);

        if (steps == -1)
            break;
        
		w_l = w_l1;
    }
    
	return w_l1;
}



    





