/***********************************************************************
* Short Title: declaration of TrialDistributor and related classes
*
* Comments:
*
* <license text>
***********************************************************************/

#include "TrialDistributor.hpp"
#include "RunPar_t.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////
// class TrialDistributor
////////////////////////////////////////////////////////////////////////

// class data

const size_t TrialDistributor::histsize = 10;

// class methods

TrialDistributor* TrialDistributor::create(RunPar_t* rp)
{
    const string& tstp = rp->trialsharing;
    if (!isType(tstp))
    {
        ostringstream emsg;
        emsg << "TrialDistributor '" << tstp << "' not defined.";
        throw runtime_error(emsg.str());
    }
    TrialDistributor* td = create(distributorsRegistry()[tstp]);
    // copy data from rp
    td->resize(rp->mol->getMaxAtomCount() + 1);
    td->tolcost = rp->tolcost;
    td->base_level = rp->base_level;
    return td;
}

TrialDistributor* TrialDistributor::create(DistributorType tp)
{
    switch (tp)
    {
        case EQUAL:     return new TrialDistributorEqual();
        case SIZE:      return new TrialDistributorSize();
        case SUCCESS:   return new TrialDistributorSuccess();
        default:
            ostringstream emsg;
            emsg << "Unhandled value of DistributorType " << tp;
            throw invalid_argument(emsg.str());
    }
    return NULL;
}

list<string> TrialDistributor::getTypes()
{
    list<string> rv;
    map<string,DistributorType>::iterator ii = distributorsRegistry().begin();
    for (; ii != distributorsRegistry().end(); ++ii)
    {
        rv.push_back(ii->first);
    }
    return rv;
}

bool TrialDistributor::isType(const string& tp)
{
    return distributorsRegistry().count(tp);
}

// public methods

bool TrialDistributor::Register()
{
    distributorsRegistry()[typeStr()] = type();
    return true;
}

void TrialDistributor::setLevelBadness(size_t lv, double bd)
{
    BadnessHistory& hist = lvbadlog[lv];
    hist.push_back(bd);
    if (hist.size() > histsize)     hist.pop_front();
}

void TrialDistributor::setLevelFillRate(size_t lv, double fr)
{
    fillrate[lv] = fr;
}

void TrialDistributor::resize(size_t sz)
{
    // clear history at every level
    lvbadlog.clear();
    lvbadlog.resize(sz);
    fillrate.resize(sz, 0.0);
    tshares.resize(sz, 0.0);
    top_level = sz - 1;
}

// protected methods

// private class methods

map<string,TrialDistributor::DistributorType>&
TrialDistributor::distributorsRegistry()
{
    static map<string,TrialDistributor::DistributorType> the_distributors;
    return the_distributors;
}


////////////////////////////////////////////////////////////////////////
// class TrialDistributorEqual
////////////////////////////////////////////////////////////////////////

void TrialDistributorEqual::share(int seasontrials)
{
    tshares = 0.0;
    if (base_level >= top_level)    return;
    tshares = 1.0 * seasontrials / (top_level - base_level);
    tshares[top_level] = 0.0;
}


////////////////////////////////////////////////////////////////////////
// class TrialDistributorSize
////////////////////////////////////////////////////////////////////////

void TrialDistributorSize::share(int seasontrials)
{
    tshares = 0.0;
    if (base_level >= top_level)    return;
    for (int lv = base_level; lv < top_level; ++lv)     tshares[lv] = lv;
    tshares[top_level] = 0.0;
    tshares *= seasontrials/tshares.sum();
}


////////////////////////////////////////////////////////////////////////
// class TrialDistributorSuccess
////////////////////////////////////////////////////////////////////////

// class data

void TrialDistributorSuccess::share(int seasontrials)
{
    tshares = 0.0;
    if (base_level >= top_level)    return;
    // calculate improvement ratio at each level
    valarray<double> scwt(0.0, size());
    for (int lv = base_level; lv <= top_level; ++lv)
    {
        BadnessHistory& hist = lvbadlog[lv];
        double tothwt = 0.0;
        int hwt = histsize - 1;
        for (int i = hist.size() - 1;  i > 0 && hwt > 0.0;  --i, --hwt)
        {
            double improvement = (hist[i] < hist[i-1]) ?
                (hist[i-1] - hist[i])/tolcost : 0;
            scwt[lv] += hwt * improvement;
            tothwt += hwt;
        }
        if (tothwt > 0.0)   scwt[lv] /= tothwt;
    }
    // split half of success weight to the levels below
    // for the last level move all success weitht to the former level
    for (int lo = 0, hi = 1; hi <= top_level; ++lo, ++hi)
    {
        double scshift = (hi < top_level) ? scwt[hi]/2 : scwt[hi];
        scwt[lo] += scshift;
        scwt[hi] -= scshift;
    }
    // calculate tshares
    // get size share weights
    valarray<double> szwt(0.0, size());
    for (int lv = base_level; lv < top_level; ++lv)     szwt[lv] = lv;
    szwt /= szwt.sum();
    // average success and size shares
    tshares = scwt + szwt;
    tshares[top_level] = 0.0;
    tshares *= seasontrials/tshares.sum();
}


////////////////////////////////////////////////////////////////////////
// registration of trial distributors
////////////////////////////////////////////////////////////////////////

namespace {

bool regTrialDistributorEqual = TrialDistributorEqual().Register();
bool regTrialDistributorSize = TrialDistributorSize().Register();
bool regTrialDistributorSuccess = TrialDistributorSuccess().Register();

}   // namespace


// End of TrialDistributor.cpp
