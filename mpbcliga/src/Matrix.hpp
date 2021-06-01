/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2006 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* template class Matrix
*
* Comments: optimized from original vector of vectors
*
***********************************************************************/

#ifndef MATRIX_HPP_INCLUDED
#define MATRIX_HPP_INCLUDED

#include <iostream>
#include <vector>

template <class T> class Matrix
{
    private:

        // Data Methods
        T* mdata;
        size_t mrows;
        size_t mcols;
        size_t msize;

    public:

        // Constructors
        Matrix() :
            mdata(NULL), mrows(0), mcols(0), msize(0)
        { }

        Matrix(const Matrix& src) :
            mdata(new T[src.mrows*src.mcols]),
            mrows(src.mrows), mcols(src.mcols), msize(src.msize)
        {
            std::copy(src.mdata, src.mdata + src.msize, mdata);
        }

        Matrix(size_t m, size_t n) :
            mrows(m), mcols(n), msize(mrows*mcols)
        {
            mdata = new T[mrows*mcols];
            std::fill(mdata, mdata + msize, T(0));
        }

        // Destructor
        virtual ~Matrix()
        {
            delete[] mdata;
        }

        // Methods
        Matrix& operator=(const Matrix& src)
        {
            if (this == &src)   return *this;
            if (msize != src.msize)
            {
                delete[] mdata;
                mdata = new T[src.msize];
            }
            std::copy(src.mdata, src.mdata + src.msize, mdata);
            mrows = src.mrows;
            mcols = src.mcols;
            msize = src.msize;
            return *this;
        }

        void fill(const T& value)
        {
            std::fill_n(mdata, msize, value);
        }

        void resize(size_t m, size_t n, const T* value=NULL)
        {
            if (m == mrows && n == mcols)   return;
            T* resized = new T[m*n];
            if (value)
            {
                std::fill(resized, resized + m*n, *value);
            }
            for (   size_t i = 0, offset = 0;
                    i != std::min(m, mrows); ++i, offset += n )
            {
                for (size_t j = 0; j != std::min(n, mcols); ++j)
                {
                    resized[offset+j] = (*this)(i, j);
                }
            }
            delete[] mdata;
            mdata = resized;
            mrows = m;
            mcols = n;
            msize = m*n;
        }

        void resize(size_t m, size_t n, const T& value)
        {
            resize(m, n, &value);
        }

        void clear()
        {
            delete[] mdata;
            mdata = NULL;
            mrows = mcols = msize = 0;
        }

        virtual std::vector<T> rowVector(size_t i) const
        {
            return std::vector<T>(mdata + mcols*i, mdata + mcols*(i+1));
        }

        virtual std::vector<T> columnVector(size_t j) const
        {
            std::vector<T> column(mrows);
            typename std::vector<T>::iterator vii;
            vii = column.begin();
            T* pt = mdata + j;
            for (; vii != column.end(); ++vii, pt += mcols)   *vii = *pt;
            return column;
        }

        inline size_t rows() const
        {
            return mrows;
        }

        inline size_t columns() const
        {
            return mcols;
        }

        inline T& operator()(size_t i, size_t j) const
        {
            return *(mdata + i*mcols + j);
        }

        Matrix<T> transposed()
        {
            Matrix<T> mxt(mcols, mrows);
            T* pt = mxt.mdata;
            for (size_t j = 0; j != mcols; ++j)
            {
                for (size_t i = 0; i != mrows; ++i, ++pt)
                {
                    *pt = (*this)(i, j);
                }
            }
            return mxt;
        }

        Matrix<T> operator*(Matrix<T>& B)
        {
            Matrix<T>& A = *this;
            Matrix<T> C(A.mrows, B.mcols);
            if (A.mcols != B.mrows)
            {
                const char* emsg = "Inconsistent matrix multiplication";
                throw emsg;
            }
            for (size_t i = 0; i != A.mrows; ++i)
            {
                for (size_t j = 0; j != B.mcols; ++j)
                {
                    T& Cij = C(i,j);
                    Cij = 0;
                    for (size_t k = 0; k != A.mcols; ++k)
                    {
                        Cij += A(i,k)*B(k,j);
                    }
                }
            }
            return C;
        }

        template <class T1>
            friend std::ostream& operator<<(std::ostream &out, Matrix<T1> &mx);

        void print()
        {
            using namespace std;
            for (size_t i = 0; i != mrows; ++i)
            {
                cout << i << ": ";
                vector<T> rowi = rowVector(i);
                for (size_t j = 0; j != mcols; ++j)
                {
                    cout << rowi[j];
                    cout << (j != mcols - 1 ? " " : "\n");
                }
            }
        }
};

template <class T> class SymmetricMatrix : public Matrix<T>
{
    public:

        // Constructors

        SymmetricMatrix() : Matrix<T>() { }
        SymmetricMatrix(const Matrix<T>& src) : Matrix<T>(src) { }
        SymmetricMatrix(size_t m) : Matrix<T>(m, m) { }

        // Overloaded methods

        SymmetricMatrix& operator=(const SymmetricMatrix& src)
        {
            this->Matrix<T>::operator=(src);
            return *this;
        }

        inline T& operator()(size_t i, size_t j) const
        {
            T& rv = (i < j) ?
                Matrix<T>::operator()(i, j) :
                Matrix<T>::operator()(j, i);
            return rv;
        }

        virtual std::vector<T> rowVector(size_t i) const
        {
            std::vector<T> rv(Matrix<T>::columns());
            for (size_t j = 0; j != Matrix<T>::columns(); ++j)
            {
                rv[j] = operator()(i, j);
            }
            return rv;
        }

        virtual std::vector<T> columnVector(size_t j) const
        {
            std::vector<T> rv = rowVector(j);
            return rv;
        }

        SymmetricMatrix<T> transposed()
        {
            return SymmetricMatrix<T>(*this);
        }

};

template <class T>
std::ostream& operator<<(std::ostream &out, Matrix<T> &mx)
{
    using namespace std;
    out << "(";
    for (size_t i = 0; i != mx.rows(); ++i)
    {
        out << "(";
        vector<T> rowi = mx.rowVector(i);
        for (size_t j = 0; j != mx.columns(); ++j)
        {
            out << rowi[j];
            out << ((j != mx.columns() - 1) ? ", " : ")");
        }
        if (i != mx.rows() - 1)     out << ", ";
    }
    out << ")\n";
    return out;
}

#endif  // MATRIX_HPP_INCLUDED
