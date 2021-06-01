/***********************************************************************
* Short Title: LIGA specific exceptions
*
* Comments:
*
* <license text>
***********************************************************************/

#ifndef EXCEPTIONS_HPP_INCLUDED
#define EXCEPTIONS_HPP_INCLUDED

#include <string>
#include <stdexcept>


class IOError : public std::runtime_error
{
    public:

        IOError(const std::string msg="") : std::runtime_error(msg)
        { }
};


class TriangulationError : public std::runtime_error
{
    public:

        TriangulationError(const std::string msg="") : std::runtime_error(msg)
        { }
};


class InvalidDistanceTable : public std::runtime_error
{
    public:

        InvalidDistanceTable(const std::string msg="") : std::runtime_error(msg)
        { }
};


class InvalidMolecule : public std::runtime_error
{
    public:

        InvalidMolecule(const std::string msg="") : std::runtime_error(msg)
        { }
};

#endif  // EXCEPTIONS_HPP_INCLUDED
