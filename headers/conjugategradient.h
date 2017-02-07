#pragma once
#include "matrix.h"
#include "vector.h"
#include "vector_operations.h"

/*
 * Solving the equation Ax=b
 * with a tolerance tol and maximum
 * iterations of maxiter. 
 * The initial guess is Vector<T>& x,
 * and it is modified to give the end result as well.
 */

template<typename T>
int cg( const Matrix<T>& A, const Vector<T>& b,
        Vector<T>& x, T tol, unsigned int maxiter)
{
	//initial guess residual_0
    Vector<T> residual_k = b - matvec(A,x);
    
	// initial guess p_0
    Vector<T> p_k = residual_k;
    
	// A * p_k
    Vector<T> Ap_k = matvec(A, p_k);
    
	int steps = 1;
    
	for (auto k = 0; k < (maxiter-1); k++)
    { 
        auto dot_Ap_k_pk = dot(Ap_k, p_k);
		auto dot_rk_rk = dot(residual_k,residual_k);
    	auto alpha_k = dot_rk_rk / dot_Ap_k_pk;
				
		// x_k+1 = x_k + alpha_k * p_k
        x = x + alpha_k * p_k ;   
	 			
		//new residual			
		residual_k = residual_k - alpha_k * Ap_k;  
		auto dot_rknew_rknew = dot(residual_k, residual_k);

        if ( dot_rknew_rknew < tol*tol)
        {
			//stop and return the num of steps taken
			return steps;        
		}
        
        auto beta_k = dot_rknew_rknew / dot_rk_rk;
				
		// p_k+1 = r_k+1 + beta_k* p_k
		p_k = residual_k + beta_k * p_k;   
		//Uodating the dot product of matrix and p_k			
        Ap_k = matvec(A, p_k);
        
		// number of steps taken
        steps +=1;    
	}
    
    return -1;
}
