#pragma once
#include "iostream"
#include "map"
#include "array" 
#include "vector.h"

typedef std::array<int,2> Array;

//forward declaration to have a declaration for friend function
template<typename T>
class Matrix; 

template<typename T>
Vector<T> matvec(const Matrix<T>& mat, const Vector<T>& vec);

template<typename T>
std::ostream& operator<<(std::ostream& os, Matrix<T>& mat);

template<typename T>
class Matrix
{
    /* The entries are stored in the following map:
     * {Row --> {Column--->Value}}
     */

    typedef std::map<unsigned int,T> Column;
    typedef std::map<unsigned int, Column> Entry;

    private:

        unsigned int rows;
        unsigned int columns;
        Entry data;
        
    public:
        //get the number of rows;
        unsigned int get_rows() const;

        //get the number of columns
        unsigned int get_columns() const;
        
        //default constructor
        Matrix();

        //a constructor that initialises a square matrix
        Matrix(const unsigned int& length);

        //a constructor that initializes matrix by length
        Matrix(const unsigned int& numrows, const unsigned int& numcols);
        
        //only for reading stuff from the array
        inline T operator[](const Array& matcoord) const ;

        //only for inputting stuff, therefore the reference
        inline T& operator[](const Array& matcoord);
        
        //ofstream operator
        friend std::ostream& operator<< <T> (std::ostream& os, Matrix<T>& mat);
        
        //non member member function for matrix and vector multiplication
        friend Vector<T> matvec<T>(const Matrix<T>& mat, const Vector<T>& vec);

				//Destructor
        ~Matrix();
};

template<typename T>
unsigned int Matrix<T>::get_rows() const { return rows; }

template<typename T>
unsigned int Matrix<T>::get_columns() const {return columns; }

template<typename T>
Matrix<T>::Matrix() : rows(0), columns(0), data() {}

template<typename T>
Matrix<T>::Matrix(const unsigned int& length) : rows(length), columns(length), data() {}

template<typename T>
Matrix<T>:: Matrix(const unsigned int& numrows, const unsigned int& numcols):
    rows(numrows), columns(numcols), data() {}

template<typename T>
T Matrix<T>::operator[](const Array& matcoord) const
{
    //use by starting with (0,0) instead of (1,1)
    if(matcoord[0] > this->get_rows()-1 || matcoord[1] > this->get_columns()-1 )
        throw;
    return data[matcoord[0]][matcoord[1]];
}

template<typename T>
T& Matrix<T>::operator[](const Array& matcoord)
{
    if(matcoord[0] > this->get_rows()-1 || matcoord[1] > this->get_columns()-1)
        throw;
    return data[matcoord[0]][matcoord[1]];
}

template<typename T>
std::ostream& operator<<(std::ostream& os, Matrix<T>& mat)
{
    Array matcoord;
    for (auto i = 0; i < mat.get_rows(); i++)
    {
        os << "[";
        for (auto j =0; j < mat.get_columns(); j++)
        {
            matcoord = {i,j};
            os << mat[matcoord] << ",\t ";
        }
        os << "]" << std::endl;
    }
    return os;
}

template<typename T>
Vector<T> matvec(const Matrix<T>& mat, const Vector<T>& vec)
{
    //returning by value instead of a new vector
    Vector<T> result(mat.get_rows());
    
		T sum;

    if ( mat.get_columns() != vec.get_length() )
    {
        std::cerr << "The matrix and vector are not the same\
                     size for multiplication" 
                  << std::endl;
				
				return result;
    }

    else
    {
         for (auto i = mat.data.begin(); i != mat.data.end(); ++i)
         {
             sum =0;

             for(auto j = i->second.begin(); j != i->second.end(); ++j)
             {
                 sum += j->second * vec[j->first];
             }
             result[i->first] = sum;
         }
    }
    return result;
}

template<typename T>
Matrix<T>::~Matrix()
{
    rows=0;
    columns=0;
    data.clear();
}



