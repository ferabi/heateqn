#include <ctime>
#include "heat2d.h"
#include <cmath>
template <unsigned int k>
inline double laplace_op_2d(const unsigned int i, const unsigned int j, const unsigned int m)
{
    if (i == j)
    {
        return -2.0;
    }
    if ( j-i == pow(m,k) )
    {
        if(k==0 && j%m ==0)
            return 0.0;
        else
            return 1.0;
    }
    else
    {
        return 0.0;
    }

}

inline Vector<double> initial_solvec_2d(const unsigned int points)
{
    Vector<double> result(points*points);
	double dx = 1.0 /(points + 1.0);
    double x;
    
	double dy = 1.0 /(points + 1.0);
    double y;
    
	for (auto i = 0; i < points; i++)
    {
        x = dx * (i+1);
        
        for (auto j = 0; j < points; j++)
        {
            y = dy * (j+1);
           	result[ points * i + j] = sin(M_PI * x) * sin(M_PI * y);
        }
    }
    
    return result;
}



//the matrix will now have m^2 rows and columns
Heat2D::Heat2D(const double in_alpha,const unsigned int in_m, const double in_dt):
    alpha(in_alpha), m(in_m), dt(in_dt), M(in_m * in_m)
{
    const clock_t begin_time = std::clock();
    
	//initialising the iteration matrix
	double dx = 1.0 / ( this->m + 1.0 );
    double dxsquared = dx * dx;
    double expr;
    
    for (auto i = 0; i < M.get_rows(); i++)
    {
        for(auto j = i; j < M.get_columns(); j++)
        {
            //its a symmetric matrix therefore the loop is
            //only half of the entire possibilities
						
		    //the identity matrix values
            double ones = (i==j)?1.0:0.0;
			//individual expression values
            expr = ones - (alpha * (dt/dxsquared) * 
						    (laplace_op_2d<0>(i,j, this->m) + 
							laplace_op_2d<1>(i,j, this->m)
							)
						  );
            
			//if the expr isn't zero, store it
			if (expr != 0)
            {
                M[{i,j}] = expr;
                M[{j,i}] = M[{i,j}];
            }
        }
    }

    std::cout << "time elapsed for making the matrix is : " 
			  << float(std::clock()- begin_time) 
			  << std::endl; 
		
		//output of the matrix built
		//std::cout << M << std::endl;
}

Vector<double> Heat2D::exact(double in_time) const
{
    //making a vector to keep the exact solution 
    Vector<double> exact_sol(m*m);
    Vector<double> in_sol = initial_solvec_2d(m);
    
    //filling the vector
    for (auto i = 0; i < m*m; i++)
    {
        auto exponent = - 2 * M_PI * M_PI* alpha * in_time;
        exact_sol[i] = exp(exponent) * in_sol[i];
    }
    
		return exact_sol;
}

Vector<double> Heat2D::solve(double time_end) const
{
    const clock_t begin_time = std::clock();
    
    //check if the input time_end is a multiple of dt
    if ( remainder(time_end,dt) > 1e-15 && time_end > 0.0)
    {
        std::cout << "Input time is not a multiple of dt" << std::endl;
        return Vector<double>(m*m);
    }
    
    //initial solution vector w_l when l=0 i.e. w_0
    Vector<double> w_l = initial_solvec_2d(m);

    //make the solution vector w_(l+1)
    Vector<double> w_l1 = w_l;

    for(auto i = 0; i < time_end/dt ; i++)
    {
        std::cout << "solving" << std::endl;
        auto steps = cg(M, w_l, w_l1, 1e-16, 200);
        std::cout << "number of steps till convergence: " << steps << std::endl;
				if (steps == -1)
        {
            return w_l1;
        }
        w_l = w_l1;
    }
    
    std::cout << "time to get the answer" << float(std::clock() - begin_time)<< std::endl;
	return w_l1;
}

