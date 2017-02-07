#pragma once
#include "iostream"
#include "vector.h"

template<typename T>
T dot(const Vector<T>& left, const Vector<T>& right)
{
    T result = 0;
    if( left.get_length() != right.get_length())
    {
        std::cerr << "Dot product of two vectors with \
                     different lengths is being attempted"
                  << std::endl;
				return result;
    }
    else
    {
        for (auto i=0; i < left.get_length(); i++)
        {
            result += left[i] * right[i];
        }
    }
    return result;
}

