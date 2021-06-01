/***********************************************************************
 * Short Title: AtomCost - functor for external atom cost calculation
 *
 * Comments: implementation of class AtomCost
 *
 * <license text>
 ***********************************************************************/

#include <cassert>
#include <cmath>
#include "AtomSequence.hpp"
#include "AtomCost.hpp"
#include "Molecule.hpp"
#include "LigaUtils.hpp"
#include "Counter.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////

double penalty(const double& dd, const double& desd)
{
    static Counter* penalty_calls = Counter::getCounter("penalty_calls");
    penalty_calls->count();
    double ddw = dd / desd;
    return ddw * ddw;
}

double penalty_gradient(const double& dd, const double& desd)
{
    return 2 * dd / desd / desd;
}


////////////////////////////////////////////////////////////////////////
// class AtomCost
////////////////////////////////////////////////////////////////////////

// constructor

AtomCost::AtomCost(const Molecule* m) : arg_atom(NULL)
{
    this->setScale(1.0);
    resetFor(m);
}

// public methods

void AtomCost::resetFor(const Molecule* m)
{
    arg_cluster = m;
    const DistanceTable& dtgt = arg_cluster->getDistanceTable();
    noCutoff();
    use_distances = !arg_cluster->getDistReuse();
    if (use_distances && dtgt.size() > useflag.size())
    {
        useflag.resize(dtgt.size(), false);
    }
}

double AtomCost::eval(const Atom_t& a, int flags)
{
    return eval(&a, flags);
}

double AtomCost::eval(const Atom_t* pa, int flags)
{
    // assign arguments
    this->arg_atom = pa;
    this->_selfcost_flag = flags & SELFCOST;
    this->_gradient_flag = flags & GRADIENT;
    // begin calculation
    resizeArrays();
    resetUseFlags();
    resetGradient();
    total_cost = 0.0;
    // selfcost is always zero in non-periodic materials
    if (this->_selfcost_flag)  return this->totalCost();
    vector<double>::iterator tgdii = target_distances.begin();
    vector<double>::iterator ptcii = partial_costs.begin();
    const DistanceTable& dtgt = arg_cluster->getDistanceTable();
    for (AtomSequenceIndex seq(arg_cluster); !seq.finished(); seq.next())
    {
        // assertion checks
        assert(tgdii < target_distances.end());
        assert(ptcii < partial_costs.end());
        // calculation
        double d = R3::distance(arg_atom->r, seq.ptr()->r);
        size_t nearidx = nearDistanceIndex(d);
        const double& dnear = dtgt[nearidx];
        const double& desd = dtgt.getesd(dnear);
        *(tgdii++) = dnear;
        double dd = dnear - d;
        double pcost = this->penaltyScaled(dd, desd);
        *(ptcii++) = pcost;
        total_cost += pcost;
        if (use_distances)
        {
            useflag[nearidx] = true;
            useflag_indices.push_back(nearidx);
            useatom_indices.push_back(seq.idx());
        }
        bool cutitoff = !this->_gradient_flag && this->apply_cutoff &&
            this->total_cost + arg_atom->Badness() > this->cutoff_cost;
        if (cutitoff)   break;
        if (this->_gradient_flag && d > NS_LIGA::eps_distance)
        {
            static R3::Vector g_dd_xyz;
            g_dd_xyz = (-1.0/d) * (arg_atom->r - seq.ptr()->r);
            double g_pcost_dd = penalty_gradient(dd, desd) * this->getScale();
            this->_gradient += g_pcost_dd * g_dd_xyz;
        }
    }
    if (apply_cutoff && arg_atom->Badness() + total_cost < lowest_cost)
    {
        lowest_cost = arg_atom->Badness() + total_cost;
        cutoff_cost = min(cutoff_cost, lowest_cost + cutoff_range);
    }
    this->_gradient_cached = this->_gradient_flag;
    return total_cost;
}


const R3::Vector& AtomCost::gradient()
{
    assert(this->_gradient_cached);
    return this->_gradient;
}


double AtomCost::lowest() const
{
    return lowest_cost;
}

double AtomCost::cutoff() const
{
    return cutoff_cost;
}

void AtomCost::setCutoff(double cf)
{
    apply_cutoff = true;
    cutoff_cost = cf;
    lowest_cost = DOUBLE_MAX;
}

double AtomCost::cutoffRange() const
{
    return cutoff_range;
}

void AtomCost::setCutoffRange(double cutrng)
{
    apply_cutoff = true;
    cutoff_range = cutrng;
    lowest_cost = DOUBLE_MAX;
}

void AtomCost::noCutoff()
{
    apply_cutoff = false;
    lowest_cost = DOUBLE_MAX;
    cutoff_cost = DOUBLE_MAX;
    cutoff_range = DOUBLE_MAX;
}

double AtomCost::totalCost() const
{
    return this->total_cost;
}

const vector<double>& AtomCost::partialCosts() const
{
    return this->partial_costs;
}

const vector<double>& AtomCost::targetDistances() const
{
    return target_distances;
}

const vector<int>& AtomCost::usedTargetDistanceIndices() const
{
    return useflag_indices;
}

const vector<int>& AtomCost::usedTargetAtomIndices() const
{
    return useatom_indices;
}


void AtomCost::setScale(double sc)
{
    mscale = sc;
}


const double& AtomCost::getScale() const
{
    return mscale;
}


double AtomCost::penaltyScaled(const double& dd, const double& esd) const
{
    double rv = penalty(dd, esd) * this->getScale();
    return rv;
}

// protected methods

const vector<Atom_t*>& AtomCost::getClusterAtoms() const
{
    return arg_cluster->atoms;
}

void AtomCost::resizeArrays()
{
    bool isresized = arg_cluster->countAtoms() == int(partial_costs.size());
    if (isresized)  return;
    partial_costs.resize(arg_cluster->countAtoms());
    target_distances.resize(arg_cluster->countAtoms());
}

void AtomCost::resetUseFlags()
{
    if (useflag_indices.empty())    return;
    for (vector<int>::iterator ii = useflag_indices.begin();
            ii != useflag_indices.end(); ++ii)
    {
        assert(*ii < int(useflag.size()));
        useflag[*ii] = false;
    }
    useflag_indices.clear();
    useatom_indices.clear();
}

void AtomCost::resetGradient()
{
    this->_gradient = 0.0;
    // gradient is always zero when calculating self cost
    this->_gradient_cached = this->_selfcost_flag;
}

size_t AtomCost::nearDistanceIndex(const double& d) const
{
    const DistanceTable& dtgt = arg_cluster->getDistanceTable();
    int idx = dtgt.find_nearest(d) - dtgt.begin();
    if (use_distances && useflag[idx])
    {
        int hi, lo, nidx = -1;
        int sz = dtgt.size();
        for (hi = idx + 1; hi < sz && useflag[hi]; ++hi)
        { }
        if (hi < sz)
        {
            nidx = hi;
        }
        for (lo = idx - 1; lo >= 0 && useflag[lo]; --lo)
        { }
        if (lo >= 0 && (nidx < 0 || d - dtgt[lo] < dtgt[nidx] - d))
        {
            nidx = lo;
        }
        idx = nidx;
    }
    return idx;
}

double AtomCost::nearDistance(const double& d) const
{
    const DistanceTable& dtgt = arg_cluster->getDistanceTable();
    double dfind = this->use_distances ?
        dtgt[this->nearDistanceIndex(d)] :
        *(dtgt.find_nearest(d));
    return dfind;
}

// End of AtomCost.cpp
