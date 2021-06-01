#ifndef LIGAUTILS_HPP_INCLUDED
#define LIGAUTILS_HPP_INCLUDED

#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <valarray>
#include <limits>
#include <algorithm>

////////////////////////////////////////////////////////////////////////
// Declarations
////////////////////////////////////////////////////////////////////////

// constants

const double DOUBLE_MAX = std::numeric_limits<double>().max();
const double DOUBLE_EPS = std::numeric_limits<double>().epsilon();

namespace NS_LIGA {
const double eps_cost = sqrt(DOUBLE_EPS);
const double eps_distance = sqrt(DOUBLE_EPS);
}   // namespace NS_LIGA

// functions

// round-off handling
bool eps_eq(const double& x, const double& y);
bool eps_gt(const double& x, const double& y);
bool eps_lt(const double& x, const double& y);

// valarray operations
double vdnorm(const std::valarray<double>&);
double vddot(const std::valarray<double>&,
        const std::valarray<double>&);
std::valarray<double> vdcross(const std::valarray<double>&,
                              const std::valarray<double>&);
// fitness
double convertCostToFitness(const double& cst);
double costToFitness(const double& cst);
std::valarray<double> costToFitness(const std::valarray<double>& vc);
template <typename Container>
    Container costToFitness(const Container& vc);

// file utilities
// similar to mkstemp(3)
std::ofstream& mktempofstream(std::ofstream& out, std::string& writefile);
bool read_header(std::istream& fid, std::string& header);
bool read_header(std::istream& fid);
template<typename T> bool read_data(std::istream& fid, std::vector<T>& v);


////////////////////////////////////////////////////////////////////////
// Definition for inline and template functions
////////////////////////////////////////////////////////////////////////

// round-off handling

inline bool eps_eq(const double& x, const double& y)
{
    return fabs(x-y) < NS_LIGA::eps_cost;
}

inline bool eps_gt(const double& x, const double& y)
{
    return x > y + NS_LIGA::eps_cost;
}

inline bool eps_lt(const double& x, const double& y)
{
    return x < y - NS_LIGA::eps_cost;
}


inline
bool isNearZeroRoundOff(const double& x, double tol=NS_LIGA::eps_cost)
{
    bool rv = (x < tol);
    return rv;
}

// fitness

inline double convertCostToFitness(const double& cst)
{
    double ftn;
    ftn = 1.0 / (cst + DOUBLE_EPS);
    return ftn;
}

inline double costToFitness(const double& cst)
{
    return convertCostToFitness(cst);
}

inline std::valarray<double> costToFitness(const std::valarray<double>& vc)
{
    std::valarray<double> vf(vc);
    double* first = &(vf[0]);
    double* last = &(vf[vf.size()]);
    std::transform(first, last, first, convertCostToFitness);
    return vf;
}

template <typename Container>
Container costToFitness(const Container& vc)
{
    Container vf(vc);
    std::transform(vf.begin(), vf.end(), vf.begin(), convertCostToFitness);
    return vf;
}

// file utilities

template<typename T>
bool read_data(std::istream& fid, std::vector<T>& v)
{
    // prepare v
    T x;
    while (fid >> x)
    {
        v.push_back(x);
    }
    return !(fid.rdstate() & std::ios::badbit);
}

#endif  // LIGAUTILS_HPP_INCLUDED
