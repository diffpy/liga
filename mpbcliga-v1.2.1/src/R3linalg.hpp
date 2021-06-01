/***********************************************************************
* Short Title: linear algebra functions on R3
*
* Comments: declaration of necessary linear algebra functions for
*     blitz::TinyVector  and  blitz::TinyMatrix
*
* <license text>
***********************************************************************/

#ifndef R3LINALG_HPP_INCLUDED
#define R3LINALG_HPP_INCLUDED

#include <blitz/array.h>

#include "Counter.hpp"

namespace R3 {

////////////////////////////////////////////////////////////////////////
// Declarations
////////////////////////////////////////////////////////////////////////

// constants

const int Ndim = 3;
using blitz::all;

// types

typedef blitz::TinyMatrix<double,Ndim,Ndim> Matrix;
typedef blitz::TinyVector<double,Ndim> Vector;

// functions

double determinant(const Matrix& A);
Matrix inverse(const Matrix& A);
Matrix transpose(const Matrix& A);
const Matrix& product(const Matrix&, const Matrix&);

template <class V> double norm(const V&);
template <class V> double distance(const V& u, const V& v);
template <class V> double dot(const V& u, const V& v);
template <class V> Vector cross(const V& u, const V& v);
const Vector& product(const Vector&, const Matrix&);

template <class M>
    bool MatricesAlmostEqual(const M& A, const M& B, double precision=0.0);

template <class V>
    bool VectorsAlmostEqual(const V& A, const V& B, double precision=0.0);


////////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////////


template <class V>
inline double norm(const V& u)
{
    static Counter* R3_norm_calls = Counter::getCounter("R3_norm_calls");
    R3_norm_calls->count();
    return sqrt(R3::dot(u, u));
}


template <class V>
inline double distance(const V& u, const V& v)
{
    static Counter* R3_distance_calls =
        Counter::getCounter("R3_distance_calls");
    R3_distance_calls->count();
    static R3::Vector duv;
    duv[0] = u[0] - v[0];
    duv[1] = u[1] - v[1];
    duv[2] = u[2] - v[2];
    return R3::norm(duv);
}


template <class V>
inline double dot(const V& u, const V& v)
{
    return (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]);
}


template <class V>
inline Vector cross(const V& u, const V& v)
{
    Vector res;
    res[0] = u[1]*v[2] - u[2]*v[1];
    res[1] = u[2]*v[0] - u[0]*v[2];
    res[2] = u[0]*v[1] - u[1]*v[0];
    return res;
}


inline const Vector& product(const Vector& u, const Matrix& M)
{
    static Vector res;
    res[0] = u[0]*M(0,0) + u[1]*M(1,0)+ u[2]*M(2,0);
    res[1] = u[0]*M(0,1) + u[1]*M(1,1)+ u[2]*M(2,1);
    res[2] = u[0]*M(0,2) + u[1]*M(1,2)+ u[2]*M(2,2);
    return res;
}


template <class M>
bool MatricesAlmostEqual(const M& A, const M& B, double precision)
{
    for (int i = 0; i < Ndim; ++i)
    {
        for (int j = 0; j < Ndim; ++j)
        {
            if (fabs(A(i,j) - B(i,j)) > precision)  return false;
        }
    }
    return true;
}


template <class V>
bool VectorsAlmostEqual(const V& u, const V& v, double precision)
{
    for (int i = 0; i < Ndim; ++i)
    {
        if (fabs(u[i] - v[i]) > precision)  return false;
    }
    return true;
}


}       // End of namespace R3

#endif  // R3LINALG_HPP_INCLUDED
