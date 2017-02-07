#pragma once
#include <iostream>
#include <memory>

//forward declaration for use in friend functions
template<typename T>
class Vector;

//multiplying with a scalar on the right
template<typename T>
Vector<T> operator*(Vector<T>& vec, const T& scalar);

//multiplying with a scalar on the left
template<typename T>
Vector<T> operator*(const T& scalar, Vector<T>& vec);

//overloading the insertion operator
template<typename T>
std::ostream& operator<<(std::ostream& os, Vector<T>& vec);

template<typename T>
class Vector
{
	private:
		unsigned int length;
		T* element;
	
	public:
		
		//get the length of the vector
		unsigned int get_length() const;
		
		//default constructors
		Vector();
        
        //A constructor that is intialised by the length
        Vector(const unsigned int& in_length);
        
        //A copy constructor that copies another vector
        Vector(const Vector<T>& other);
       
        //constructor that takes initializer list
        Vector(const std::initializer_list<T>& list);

        //copy assignment from one vector to the other
        Vector<T>& operator=(const Vector<T>& other);

        //move assignment from one vector to the other
        Vector<T>& operator=( Vector<T>&& other);

        //overloading of compound assignment operators

        Vector<T>& operator+=(const Vector<T>& other);

        Vector<T>& operator-=(const Vector<T>& other);

        //overloading of binary arithmatic operators

        Vector<T> operator+(const Vector<T>& other) const;

        Vector<T> operator-(const Vector<T>& other) const;
       
        //multiplying with a scalar, overloading compound operator
        Vector<T>& operator*=(const T& scalar);

        //multiplying with a scalar on the right
        friend Vector operator*<T>(Vector& vec, const T& scalar);

        //multiplying with a scalar on the left
        friend Vector operator*<T>(const T& scalar, Vector& vec);

        //ofstream operator for the vector class
        friend std::ostream& operator<< <T> (std::ostream& os, Vector<T>& vec);

        //overloading [] operator for getting the elements of
        //the vector directly
        inline T operator[](const unsigned int& i) const;
        
        //For writing with [] operator
        inline T& operator[](const unsigned int& i);
        
        //Destructor
        ~Vector();

};


template<typename T>
unsigned int Vector<T>::get_length() const { return length; }

template<typename T>
Vector<T>::Vector() : length(0), element(nullptr) {}

template<typename T>
Vector<T>::Vector(const unsigned int& in_length) : length(in_length), element( new T[length]() ) {}

template<typename T>
Vector<T>::Vector(const Vector<T>& other) : Vector(other.length) 
{
    for (auto i = 0; i < length; i++)
    {
        element[i] = other.element[i];
    }
}

template<typename T>
Vector<T>::Vector(const std::initializer_list<T>& list) : Vector((unsigned int)list.size())
{
    std::uninitialized_copy(list.begin(),list.end(),element);
}


template<typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& other)
{
    if(this != &other)
    {
        delete[] element;
        element = new T[other.length];
        length = other.length;
        for(auto i=0; i < other.length; i++)
        {
            element[i] = other.element[i];
        }
    }
    return *this;
}

template<typename T>
Vector<T>& Vector<T>::operator=(Vector<T>&& other)
{
    if(this != &other)
    {
        delete[] element;
        length = other.length;
        element = other.element;
        other.length = 0;
        other.element = nullptr ;
    }
    return *this;
}

template<typename T>
Vector<T>& Vector<T>::operator+=(const Vector<T>& other)
{ 
    if(this->get_length() != other.get_length())
    {
        std::cerr << "Vectors being added are not \
                     of the same length" 
                  << std::endl;
    }
    else
    {
        for(auto i = 0; i < this->get_length(); i++)
        {
            this->element[i] += other.element[i];
        }
    }
    return *this;
}

template<typename T>
Vector<T>& Vector<T>::operator-=(const Vector<T>& other)
{
    if(this->get_length() != other.get_length())
    {
        std::cerr << "Vectors being subtracted are not\
                     of the same length" 
                  << std::endl;
    }
    else
    {
        for(auto i = 0; i < this->get_length(); i++)
        {
            this->element[i] -= other.element[i];
        }
    }
    return *this;
}


template<typename T>
Vector<T> Vector<T>::operator+(const Vector<T>& other) const
{
    return Vector<T>(*this) += other;
}
            
template<typename T>
Vector<T> Vector<T>::operator-(const Vector<T>& other) const
{
    return Vector<T>(*this) -= other;
}

template<typename T>
Vector<T>& Vector<T>::operator*=(const T& scalar)
{
    for (auto i =0; i < this->get_length(); i++)
    {
        this->element[i] = this->element[i] * scalar;
    }
    return *this;
}

template<typename T>
Vector<T> operator*(Vector<T>& vec, const T& scalar)
{
    return Vector<T>(vec) *= scalar; 
}

template<typename T>
Vector<T> operator*(const T& scalar, Vector<T>& vec)
{
    return vec * scalar;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, Vector<T>& vec)
{
    for(auto i=0; i < vec.get_length(); i++)
    {
        os << vec[i] << std::endl;
    }
    
		return os;
}
template<typename T>
T Vector<T>::operator[](const unsigned int& i) const
{
    try
    {
        if (i > this->get_length()-1)
            throw "Accessing element out of bounds of vector";
    }
    catch(const char* exception)
    {
        std::cerr << "Error: " << exception << std::endl;
    }

    return element[i];
}
template<typename T>
T& Vector<T>::operator[](const unsigned int& i)
{
    try
    {
        if (i > this->get_length()-1)
            throw "Accessing element out of bounds of vector";
    }
    catch(const char* exception)
    {
        std::cerr << "Error: " << exception << std::endl;
    }

    return element[i];
}
template<typename T>
Vector<T>::~Vector()
{
    delete[] element;
    length = 0;
}
