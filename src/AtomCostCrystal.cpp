/***********************************************************************
 * Short Title: AtomCostCrystal - atom cost calculation with PBC
 *
 * Comments: implementation of class AtomCostCrystal.
 *      AtomCostCrystal inherits from AtomCost.
 *
 * <license text>
 ***********************************************************************/

#include <cassert>
#include <cmath>
#include "AtomCostCrystal.hpp"
#include "AtomSequence.hpp"
#include "PointsInSphere.hpp"
#include "LigaUtils.hpp"
#include "Lattice.hpp"
#include "Crystal.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////
// class AtomCostCrystal
////////////////////////////////////////////////////////////////////////

// constructor

AtomCostCrystal::AtomCostCrystal(const Crystal* cluster) : AtomCost(cluster)
{
    use_distances = false;
}

// public methods - overloaded

void AtomCostCrystal::resetFor(const Molecule* clust)
{
    this->AtomCost::resetFor(clust);
    this->arg_cluster = static_cast<const Crystal*>(clust);
    assert(this->arg_cluster == this->AtomCost::arg_cluster);
    assert(!this->use_distances);
    this->_rmax = arg_cluster->getRmax();
    pair<double,double> rext =
        arg_cluster->getRExtent(0.0, this->arg_cluster->getRmax());
    _sph.reset(new PointsInSphere(rext.first, rext.second,
                arg_cluster->getLattice()) );
}


double AtomCostCrystal::eval(const Atom_t& a, int flags)
{
    return eval(&a, flags);
}


double AtomCostCrystal::eval(const Atom_t* pa, int flags)
{
    // assign arguments
    this->arg_atom = pa;
    this->crst_atom = this->arg_atom;
    this->_selfcost_flag = flags & SELFCOST;
    this->_gradient_flag = flags & GRADIENT;
    // begin calculation
    resizeArrays();
    resetGradient();
    vector<double>::iterator ptcii = this->partial_costs.begin();
    vector<int>::iterator pcntii = this->pair_counts.begin();
    // reset result data
    this->total_cost = 0.0;
    this->total_pair_count = 0;
    pair<double,int> costcount(0.0, 0);
    // short circuit for selfcost
    if (this->_selfcost_flag)
    {
        R3::Vector zeros(0.0, 0.0, 0.0);
        costcount = this->pairCostCount(zeros);
        assert(0.0 == R3::norm(this->gradient()));
        this->total_cost = costcount.first;
        this->total_pair_count = costcount.second;
        return this->totalCost();
    }
    // Cartesian separation vector mapped to unit cell
    R3::Vector rcv;
    for (AtomSequenceIndex seq(arg_cluster); !seq.finished(); seq.next())
    {
        // assertion checks
        assert(ptcii < this->partial_costs.end());
        assert(pcntii < this->pair_counts.end());
        // calculation
        this->crst_atom = seq.ptr();
        rcv = this->arg_atom->r - this->crst_atom->r;
        costcount = (this->arg_atom == this->crst_atom) ? make_pair(0.0, 0) :
                this->pairCostCount(rcv);
        *(ptcii++) = costcount.first;
        *(pcntii++) = costcount.second;
        this->total_cost += costcount.first;
        this->total_pair_count += costcount.second;
        bool cutitoff = !this->_gradient_flag && this->apply_cutoff &&
            this->total_cost + this->arg_atom->Badness() > this->cutoff_cost;
        if (cutitoff)   break;
    }
    bool islowestcost = this->apply_cutoff &&
        this->arg_atom->Badness() + this->total_cost < this->lowest_cost;
    if (islowestcost)
    {
        this->lowest_cost = this->arg_atom->Badness() + this->total_cost;
        this->cutoff_cost =
            min(this->cutoff_cost, this->lowest_cost + this->cutoff_range);
    }
    this->_gradient_cached = this->_gradient_flag;
    return this->totalCost();
}


int AtomCostCrystal::totalPairCount() const
{
    return this->total_pair_count;
}


const vector<int>& AtomCostCrystal::pairCounts() const
{
    return this->pair_counts;
}


// public methods - specific

// Evaluate cost of Cartesian separation between two sites
// The separation vector is transformed to unit cell vector.
// This method also sets

pair<double,int>
AtomCostCrystal::pairCostCount(const R3::Vector& cv)
{
    const Lattice& lat = arg_cluster->getLattice();
    static R3::Vector ucv;
    ucv = lat.ucvCartesian(cv);
    R3::Vector rc_dd;
    double paircost = 0.0;
    int paircount = 0;
    int loopscale = this->_selfcost_flag ? 1 : 2;
    double penaltyscale = loopscale * this->getScale();
    for (_sph->rewind(); !_sph->finished(); _sph->next())
    {
        rc_dd = ucv + lat.cartesian(_sph->mno());
        double d = R3::norm(rc_dd);
        if (d > this->_rmax)    continue;
        if (this->_selfcost_flag && d == 0.0)  continue;
        const pair<double,double>& ddesd = this->pairDistanceDifference(d);
        paircost += penalty(ddesd.first, ddesd.second) * penaltyscale;
        paircount += loopscale;
        if (this->_gradient_flag && d > NS_LIGA::eps_distance)
        {
            static R3::Vector g_dd_xyz;
            g_dd_xyz = (-1.0/d) * rc_dd;
            double g_pcost_dd =
                penalty_gradient(ddesd.first, ddesd.second) * penaltyscale;
            this->_gradient += g_pcost_dd * g_dd_xyz;
        }
    }
    return make_pair(paircost, paircount);
}

// protected methods

const pair<double,double>&
AtomCostCrystal::pairDistanceDifference(const double& d) const
{
    static pair<double,double> rv;
    const DistanceTable& dtgt = arg_cluster->getDistanceTable();
    double dn = this->nearDistance(d);
    rv.first = dn - d;
    rv.second = dtgt.getesd(dn);
    return rv;
}


void AtomCostCrystal::resizeArrays()
{
    size_t sz = arg_cluster->countAtoms();
    bool isresized = (sz == this->partial_costs.size());
    if (isresized)  return;
    this->partial_costs.resize(sz);
    this->pair_counts.resize(sz);
}

// End of file
