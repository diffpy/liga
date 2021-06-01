/***********************************************************************
* Short Title: random number generation used for LIGA algorithm
*
* Comments:
*
* <license text>
***********************************************************************/

#include <cassert>
#include <map>
#include <memory>
#include <algorithm>
#include "Random.hpp"

using namespace std;

namespace NS_LIGA {


////////////////////////////////////////////////////////////////////////
// Constants
////////////////////////////////////////////////////////////////////////

// RNG type and state
gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);


////////////////////////////////////////////////////////////////////////
// Function definitions
////////////////////////////////////////////////////////////////////////

R3::Vector randomDir()
{
    double phi = 2*M_PI*randomFloat();
    double z = 2*randomFloat() - 1.0;
    double rxy = sqrt(1.0 - z*z);
    R3::Vector rdir(rxy*cos(phi), rxy*sin(phi), z);
    return rdir;
}

const PickType& randomPickFew(size_t k, size_t N)
{
    // cannot pick more items than available
    if (k > N)
    {
        throw out_of_range("randomPickFew(): too many items to pick");
    }
    static PickType pick;
    pick.resize(k);
    // check trivial case
    if (k == 0)
    {
        return pick;
    }
    size_t Nrem = N;
    map<size_t,size_t> tr;
    PickType::iterator vi = pick.begin();
    for (size_t i = 0; i < k; ++i, --Nrem)
    {
        size_t j = randomInt(Nrem);
        for (size_t c = 0;  tr.count(j) == 1;  j = tr[j], ++c)
        {
            // number of translations must be smaller than N
            assert(c < N);
        }
        *(vi++) = j;
        tr[j] = Nrem - 1;
    }
    return pick;
}

const PickType& randomPickWithRepeat(size_t k, size_t N)
{
    static PickType pick;
    pick.resize(k);
    for (PickType::iterator vi = pick.begin(); vi != pick.end(); ++vi)
    {
        *vi = randomInt(N);
    }
    return pick;
}

////////////////////////////////////////////////////////////////////////
// class RandomWeighedGenerator  -  definitions
////////////////////////////////////////////////////////////////////////

// class methods

RandomWeighedGenerator* RandomWeighedGenerator::instance()
{
    static auto_ptr<RandomWeighedGenerator> the_rwg(new RandomWeighedGenerator);
    return the_rwg.get();
}

// constructors

RandomWeighedGenerator::RandomWeighedGenerator()
{
    double w = 1.0;
    setWeights(&w, &w + 1);
}

// public methods

const PickType& RandomWeighedGenerator::weighedPick(size_t k) const
{
    size_t N = numChoices();
    // check arguments
    if (k > N)
    {
        const char* emsg = "weighedPick(): too many items to pick";
        throw out_of_range(emsg);
    }
    this->_pick.resize(k);
    // check trivial case
    if (k == 0)
    {
        return this->_pick;
    }
    // here we can do some real work
    // create working copies of _cumul_weight and _total_weight
    double cwts[N];
    copy(this->_cumul_weight.begin(), this->_cumul_weight.end(), cwts);
    double totwts = this->_total_weight;
    // save original indices
    this->_pick_index.resize(N);
    for (size_t i = 0; i != N; ++i)
    {
        this->_pick_index[i] = i;
    }
    // main loop
    PickType::iterator vi = this->_pick.begin();
    for (size_t i = 0, Nrem = N; i != k; ++i, --Nrem, ++vi)
    {
        assert(Nrem > 0);
        // probabilities are uniform when totwts == 0.0
        if (totwts == 0.0)
        {
            size_t isel = randomInt(Nrem);
            *vi = this->_pick_index[isel];
            // overwrite this element with the last number
            this->_pick_index[isel] = this->_pick_index[Nrem-1];
        }
        // otherwise we need to do binary search on cwts
        else
        {
            double* p;
            p = lower_bound(cwts, cwts + Nrem, totwts*randomFloat());
            size_t isel = p - cwts;
            *vi = this->_pick_index[isel];
            // overwrite this element with the last number
            this->_pick_index[isel] = this->_pick_index[Nrem-1];
            // and update cwts
            totwts = (p == cwts) ? 0.0 : *(p-1);
            for (size_t j = isel; j != Nrem - 1; ++j, ++p)
            {
                totwts += this->_weight[this->_pick_index[j]];
                *p = totwts;
            }
        }
    }
    return this->_pick;
}

size_t RandomWeighedGenerator::weighedInt() const
{
    // check _total_weight
    if (this->_total_weight == 0.0)
    {
        return randomInt(numChoices());
    }
    // here this->_total_weight > 0.0
    double ranval = this->_total_weight*randomFloat();
    size_t idx;
    idx = lower_bound(this->_cumul_weight.begin(), this->_cumul_weight.end(),
            ranval) - this->_cumul_weight.begin();
    return idx;
}


////////////////////////////////////////////////////////////////////////
// Function definitions
////////////////////////////////////////////////////////////////////////


}   // namespace NS_LIGA

// End of file
