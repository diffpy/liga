/***********************************************************************
* Short Title: random number generation used for LIGA algorithm
*
* Comments:
*
* <license text>
***********************************************************************/

#ifndef RANDOM_HPP_INCLUDED
#define RANDOM_HPP_INCLUDED

#include <numeric>
#include <vector>
#include <stdexcept>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "R3linalg.hpp"

namespace NS_LIGA {

////////////////////////////////////////////////////////////////////////
// Types and global variables
////////////////////////////////////////////////////////////////////////

typedef std::vector<size_t> PickType;

// random number generator state
extern gsl_rng* rng;

////////////////////////////////////////////////////////////////////////
// Function declarations
////////////////////////////////////////////////////////////////////////

inline void randomSeed(unsigned long int seed);
inline double randomFloat();
inline size_t randomInt(size_t N);
inline int plusminus();
inline double randomBeta(const double& a, const double& b);
R3::Vector randomDir();
const PickType& randomPickFew(size_t k, size_t N);
const PickType& randomPickWithRepeat(size_t k, size_t N);
template <class T>
    const PickType& randomWeighedPick(size_t k, size_t N, T* wts);
template <class T>
    size_t randomWeighedInt(size_t N, T* wts);
template <class Container> void randomShuffle(Container& c);

////////////////////////////////////////////////////////////////////////
// Class declarations
////////////////////////////////////////////////////////////////////////

class RandomWeighedGenerator
{
    public:

        // class methods
        static RandomWeighedGenerator* instance();

        // constructors
        RandomWeighedGenerator();
        template <class Iter> RandomWeighedGenerator(Iter first, Iter last);

        // methods
        template <class Iter> void setWeights(Iter first, Iter last);
        inline size_t numChoices() const;
        const PickType& weighedPick(size_t k) const;
        size_t weighedInt() const;

    private:

        // data
        std::vector<double> _weight;
        std::vector<double> _cumul_weight;
        double _total_weight;
        mutable PickType _pick_index;
        mutable PickType _pick;
};

////////////////////////////////////////////////////////////////////////
// Definitions of inline functions
////////////////////////////////////////////////////////////////////////

inline void randomSeed(unsigned long int seed)
{
    gsl_rng_set(rng, seed);
}

inline double randomFloat()
{
    return gsl_rng_uniform(rng);
}

inline size_t randomInt(size_t N)
{
    return gsl_rng_uniform_int(rng, N);
}

inline int plusminus()
{
    return (randomInt(2) != 0) ? +1 : -1;
}

inline double randomBeta(const double& a, const double& b)
{
    return gsl_ran_beta(rng, a, b);
}

////////////////////////////////////////////////////////////////////////
// Definitions of template functions
////////////////////////////////////////////////////////////////////////

template <class T>
const PickType& randomWeighedPick(size_t k, size_t N, T* wts)
{
    RandomWeighedGenerator* rwg = RandomWeighedGenerator::instance();
    rwg->setWeights(wts, wts + N);
    return rwg->weighedPick(k);
}

template <class T>
size_t randomWeighedInt(size_t N, T* wts)
{
    RandomWeighedGenerator* rwg = RandomWeighedGenerator::instance();
    rwg->setWeights(wts, wts + N);
    return rwg->weighedInt();
}

template <class Container>
void randomShuffle(Container& c)
{
    size_t N = c.size();
    for (size_t k = N; k > 1; --k)
    {
        size_t i = randomInt(k);
        std::swap(c[i], c[k-1]);
    }
}

////////////////////////////////////////////////////////////////////////
// class RandomWeighedGenerator  -  definitions of template methods
////////////////////////////////////////////////////////////////////////

// constructor

template <class Iter>
RandomWeighedGenerator::RandomWeighedGenerator(Iter first, Iter last)
{
    setWeights(first, last);
}

// public template methods

template <class Iter>
void RandomWeighedGenerator::setWeights(Iter first, Iter last)
{
    this->_weight.resize(last - first);
    std::vector<double>::iterator wii = this->_weight.begin();
    for (Iter ii = first; ii != last; ++ii, ++wii)
    {
        if (*ii < 0.0)
        {
            const char* emsg = "setWeights(): negative choice probability";
            throw std::out_of_range(emsg);
        }
        *wii = *ii;
    }
    this->_cumul_weight.resize(this->_weight.size());
    partial_sum(this->_weight.begin(), this->_weight.end(),
            this->_cumul_weight.begin());
    this->_total_weight = this->_cumul_weight.back();
}

inline size_t RandomWeighedGenerator::numChoices() const
{
    return this->_weight.size();
}


}   // namespace NS_LIGA

#endif  // RANDOM_HPP_INCLUDED
