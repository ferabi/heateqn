#include "vector.h"
#include "vector_operations.cpp"
#include <iostream>
#include "matrix.h"
#include "conjugategradient.cpp"
#include "heat2d.h"
#include <fstream>
int main()
{
    /*
    //initialising by len
    Vector<int> a((int)23);
    std::cout << a.get_length() << std::endl;
    //intialising by a vector
    Vector<int> b(a);
    std::cout << b.get_length() << std::endl;
    //intialising by a move constructor
    //Vector<int> c(a+b);

    //initialising by list
    Vector<int> c({1,2,3});
    std::cout << c.get_length() << std::endl;

    //copy assignment
    //c=b;
    //move assignment
    //c= std::move(b);
    //std::cout << b.get_length() << std::endl;

    //overloading of compound operators
    Vector<int> d({1,2,3});
    Vector<int> e({2,4,6,8});
    Vector<int> f({2,4,6,8});
    Vector<int> g = d * 2;
    g *= 2;
    std::cout << c[2] << std::endl;
    std::cout << g[2] << std::endl;
    std::cout << dot(c,d)<< std::endl;

    return 0;
    */
    //matrix constructor
    Vector<float> b({1,2});
    Matrix<float> A(2,2);
    Vector<float> x;
    A[{0,0}] = 4;
    A[{0,1}] = 1;
    A[{1,0}] = 1;
    A[{1,1}] = 3;
    x = matvec(A,b);
    std::cout << x << std::endl;
    //A[{1,1}] = 4;
    //A[{1,2}] = 2;
    //A[{2,0}] = 1;
    //A[{2,1}] = 2;
    //A[{2,2}] = 3;
   
    float tol = 1e-10;
    unsigned int max = 100;
    int j = cg(A, b, x, tol, max);
    //Vector<float> f =  x * tol;
    
    std::cout << x[0]<< std::endl;
    std::cout << x[1]<< std::endl;
    //std::cout << x[2]<< std::endl;
    std::cout << j << std::endl;
    //Matrix<float> i(3);
    //i(3);
    std::ofstream data;
    Heat2D heat_sol(0.3125, 99,0.001);
    Vector<double> v =heat_sol.solve(1);
    data.open("datasolve.txt");
    std::cout  << "from numerical solving:" << std::endl;
    data << v << std::endl;
    data.close();
    Vector<double> g = heat_sol.exact(1);
    std::cout << "from exact solution" << std::endl;
    data.open("dataexact.txt");
    data << g << std::endl;
    data.close();
    //std::cout << remainder(1,0.1) << std::endl;
    return 0;

}
