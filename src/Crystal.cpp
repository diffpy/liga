/***********************************************************************
* Short Title: object definitions for Biosphere Genetic Algorithm
*
* Comments:
*
* $Id$
*
* <license text>
***********************************************************************/

#include <sstream>
#include <stdexcept>
#include <boost/foreach.hpp>
#include "Crystal.hpp"
#include "Lattice.hpp"
#include "Atom_t.hpp"
#include "AtomCostCrystal.hpp"
#include "AtomOverlapCrystal.hpp"
#include "LigaUtils.hpp"
#include "AtomSequence.hpp"

using namespace std;
using namespace NS_LIGA;

////////////////////////////////////////////////////////////////////////
// class Crystal
////////////////////////////////////////////////////////////////////////


// constructors

Crystal::Crystal() : Molecule()
{
    init();
}

Crystal::Crystal(const Crystal& crs) : Molecule()
{
    init();
    *this = crs;
}

Crystal::~Crystal()
{ }


// public methods


Crystal& Crystal::operator=(const Molecule& mol)
{
    const Crystal& crs = dynamic_cast<const Crystal&>(mol);
    return operator=(crs);
}


Crystal& Crystal::operator=(const Crystal& crs)
{
    if (this == &crs)   return *this;
    // execute base assignment
    Molecule::operator=(crs);
    // copy Crystal specific data
    this->_lattice = crs._lattice;
    this->_lattice_max_ucd = crs._lattice_max_ucd;
    this->_rmin = crs._rmin;
    this->_rmax = crs._rmax;
    this->_cost = crs._cost;
    this->_full_distance_table = crs._full_distance_table;
    this->pmx_pair_counts = crs.pmx_pair_counts;
    this->_count_pairs = crs._count_pairs;
    this->_cost_data_cached = crs._cost_data_cached;
    return *this;
}


Molecule* Crystal::clone() const
{
    Molecule* pclone = new Crystal(*this);
    return pclone;
}


void Crystal::setDistanceTable(const DistanceTable& dtbl)
{
    vector<double> dtbl_unique = dtbl.unique();
    this->_full_distance_table.reset(new DistanceTable(dtbl_unique));
    this->cropDistanceTable();
    this->uncacheCostData();
    this->CheckIntegrity();
}

void Crystal::setDistReuse(bool flag)
{
    if (!flag)
    {
        const char* emsg = "Crystal requires distreuse is true.";
        throw range_error(emsg);
    }
    this->Molecule::setDistReuse(flag);
}

void Crystal::setLattice(const Lattice& lat)
{
    this->_lattice.reset(new Lattice(lat));
    this->_lattice_max_ucd = this->_lattice->ucMaxDiagonalLength();
    this->uncacheCostData();
}

const Lattice& Crystal::getLattice() const
{
    return *_lattice;
}

void Crystal::setRmax(double rmax)
{
    this->_rmax = rmax;
    this->cropDistanceTable();
    this->uncacheCostData();
}

double Crystal::getRmax() const
{
    double rv = (this->_rmax > 0) ? this->_rmax :
        (this->_full_distance_table->maxDistance() + eps_distance);
    return rv;
}


// r-range extended by circum diameter of the primitive cell
// The equivalent unit cell vectors are mapped to fractional coordinates
// from 0 <= xi < 1.  Due to round-off errors, some vectors can end very close
// to 1.  Therefore the r-extent has to equal the largest diagonal in the cell.

pair<double,double> Crystal::getRExtent(double rmin, double rmax) const
{
    return make_pair(rmin - _lattice_max_ucd, rmax + _lattice_max_ucd);
}


double Crystal::cost() const
{
    if (!this->_cost_data_cached)   recalculate();
    return this->Molecule::cost();
}


int Crystal::countPairs() const
{
    if (!this->_cost_data_cached)   recalculate();
    return this->_count_pairs;
}


void Crystal::recalculate() const
{
    // reset molecule
    this->ResetBadness();
    this->pmx_partial_costs.fill(0.0);
    this->pmx_pair_counts.fill(0);
    this->_count_pairs = 0;
    // reset all atoms
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
        seq.ptr()->ResetBadness();
    }
    AtomCostCrystal* atomcost;
    atomcost = static_cast<AtomCostCrystal*>(getAtomCostCalculator());
    // fill in self contributions
    if (this->countAtoms())
    {
        atomcost->eval(this->getAtom(0), AtomCost::SELFCOST);
    }
    double diagpaircost = atomcost->totalCost();
    int diagpaircount = atomcost->totalPairCount();
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
        int idx = seq.ptr()->pmxidx;
        this->pmx_partial_costs(idx, idx) = diagpaircost;
        this->IncBadness(diagpaircost);
        seq.ptr()->IncBadness(diagpaircost);
        this->pmx_pair_counts(idx, idx) = diagpaircount;
        this->_count_pairs += diagpaircount;
    }
    // fill off-diagonal elements
    R3::Vector rc_dd;
    pair<double,int> costcount;
    for (AtomSequence seq0(this); !seq0.finished(); seq0.next())
    {
        AtomSequence seq1 = seq0;
        seq1.next();
        for (; !seq1.finished(); seq1.next())
        {
            int idx0 = seq0.ptr()->pmxidx;
            int idx1 = seq1.ptr()->pmxidx;
            assert(idx0 != idx1);
            rc_dd = seq1.ptr()->r - seq0.ptr()->r;
            costcount = atomcost->pairCostCount(rc_dd);
            // pmx is SymmetricMatrix, it is sufficient to make one assignment
            this->pmx_partial_costs(idx0, idx1) = costcount.first;
            double halfcost = costcount.first / 2;
            this->IncBadness(costcount.first);
            seq0.ptr()->IncBadness(halfcost);
            seq1.ptr()->IncBadness(halfcost);
            // pmx is SymmetricMatrix, it is sufficient to make one assignment
            this->pmx_pair_counts(idx0, idx1) = costcount.second;
            this->_count_pairs += costcount.second;
        }
    }
    this->_cost_data_cached = true;
    this->recalculateOverlap();
}


AtomCost* Crystal::getAtomCostCalculator() const
{
    static AtomCostCrystal the_acc(this);
    the_acc.resetFor(this);
    return &the_acc;
}


AtomCost* Crystal::getAtomOverlapCalculator() const
{
    static AtomOverlapCrystal the_overlap_calculator(this);
    the_overlap_calculator.resetFor(this);
    return &the_overlap_calculator;
}


void Crystal::Shift(const R3::Vector& drc)
{
    this->Molecule::Shift(drc);
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
	Atom_t* pa = seq.ptr();
  	pa->r = this->ucvCartesianAdjusted(pa->r);
    }
}


void Crystal::Clear()
{
    this->Molecule::Clear();
    uncacheCostData();
}


const pair<int*,int*>& Crystal::Evolve(const int* est_triang)
{
    const pair<int*,int*>& acc_tot = this->Molecule::Evolve(est_triang);
    this->shiftToOrigin();
    return acc_tot;
}


void Crystal::Degenerate(int Npop)
{
    this->Molecule::Degenerate(Npop);
    this->shiftToOrigin();
}

// protected methods

void Crystal::AddInternal(Atom_t* pa)
{
    Molecule::AddInternal(pa);
    // shift to unit cell and take care of zero round-off
    pa->r = this->ucvCartesianAdjusted(pa->r);
}


void Crystal::addNewAtomPairs(Atom_t* pa)
{
    if (!this->_cost_data_cached)   return;
    AtomCostCrystal* atomcost;
    atomcost = static_cast<AtomCostCrystal*>(getAtomCostCalculator());
    atomcost->eval(pa);
    // store partial costs
    const vector<double>& ptcs = atomcost->partialCosts();
    const vector<int>& pcnt = atomcost->pairCounts();
    assert(atoms.size() == ptcs.size());
    assert(atoms.size() == pcnt.size());
    for (AtomSequenceIndex seq(this); !seq.finished(); seq.next())
    {
	double paircost = ptcs[seq.idx()];
	int idx0 = pa->pmxidx;
	int idx1 = seq.ptr()->pmxidx;
        assert(idx0 != idx1);
	this->pmx_partial_costs(idx0, idx1) = paircost;
	double paircosthalf = paircost / 2.0;
	seq.ptr()->IncBadness(paircosthalf);
	pa->IncBadness(paircosthalf);
        int paircount = pcnt[seq.idx()];
        this->pmx_pair_counts(idx0, idx1) = paircount;
    }
    this->IncBadness(atomcost->totalCost());
    this->_count_pairs += atomcost->totalPairCount();
    // add self contribution:
    // calculates self cost and self pair count
    double diagpaircost;
    int diagpaircount;
    if (this->atoms.empty())
    {
        R3::Vector zeros(0.0, 0.0, 0.0);
        Atom_t adummy("", 0.0, 0.0, 0.0);
        atomcost->eval(adummy, AtomCost::SELFCOST);
        diagpaircost = atomcost->totalCost();
        diagpaircount = atomcost->totalPairCount();
    }
    else
    {
        int idx1 = atoms.front()->pmxidx;
        diagpaircost = this->pmx_partial_costs(idx1, idx1);
        diagpaircount = this->pmx_pair_counts(idx1, idx1);
    }
    int idx0 = pa->pmxidx;
    this->pmx_partial_costs(idx0, idx0) = diagpaircost;
    pa->IncBadness(diagpaircost);
    this->IncBadness(diagpaircost);
    this->pmx_pair_counts(idx0, idx0) = diagpaircount;
    this->_count_pairs += diagpaircount;
    // take care of small round offs
    if (this->Badness() < NS_LIGA::eps_cost)    this->ResetBadness();
    // add overlap contributions
    this->applyOverlapContributions(pa, ADD);
}


void Crystal::removeAtomPairs(Atom_t* pa)
{
    // remove associated pair costs and pair counts
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
	int idx0 = pa->pmxidx;
	int idx1 = seq.ptr()->pmxidx;
	// remove pair costs
	double paircost = pmx_partial_costs(idx0, idx1);
	double paircosthalf = paircost/2.0;
	pa->DecBadness(paircosthalf);
	seq.ptr()->DecBadness(paircosthalf);
	this->DecBadness(paircost);
        // remove pair counts
        this->_count_pairs -= this->pmx_pair_counts(idx0, idx1);
    }
    if (this->Badness() < NS_LIGA::eps_cost)    this->ResetBadness();
    assert(this->_count_pairs >= 0);
    // remove overlap contributions
    this->applyOverlapContributions(pa, REMOVE);
}


const Crystal::TriangulationAnchor&
Crystal::getLineAnchor(const RandomWeighedGenerator& rwg)
{
    assert(countAtoms() >= 1);
    static TriangulationAnchor anch;
    anch.count = 2;
    anch.B0 = anyOffsetAtomSite(rwg);
    anch.B1 = anyOffsetAtomSite(rwg);
    return anch;
}


const Crystal::TriangulationAnchor&
Crystal::getPlaneAnchor(const RandomWeighedGenerator& rwg)
{
    assert(countAtoms() >= 1);
    static TriangulationAnchor anch;
    anch.count = 3;
    anch.B0 = anyOffsetAtomSite(rwg);
    anch.B1 = anyOffsetAtomSite(rwg);
    anch.B2 = anyOffsetAtomSite(rwg);
    return anch;
}


const Crystal::TriangulationAnchor&
Crystal::getPyramidAnchor(const RandomWeighedGenerator& rwg)
{
    return Crystal::getPlaneAnchor(rwg);
}


void Crystal::resizePairMatrices(int sz)
{
    Molecule::resizePairMatrices(sz);
    size_t sznew = this->pmx_partial_costs.rows();
    this->pmx_pair_counts.resize(sznew, sznew, 0);
}


boost::python::object Crystal::newDiffPyStructure() const
{
    boost::python::object stru;
    stru = this->Molecule::newDiffPyStructure();
    const Lattice& L = this->getLattice();
    boost::python::object lattice;
    lattice = stru.attr("lattice");
    boost::python::call_method<void>(lattice.ptr(), "setLatPar",
            L.a(), L.b(), L.c(), L.alpha(), L.beta(), L.gamma());
    return stru;
}

void Crystal::setFromDiffPyStructure(boost::python::object stru)
{
    namespace python = boost::python;
    boost::python::object lattice;
    lattice = stru.attr("lattice");
    double a, b, c, alpha, beta, gamma;
    a = python::extract<double>(lattice.attr("a"));
    b = python::extract<double>(lattice.attr("b"));
    c = python::extract<double>(lattice.attr("c"));
    alpha = python::extract<double>(lattice.attr("alpha"));
    beta = python::extract<double>(lattice.attr("beta"));
    gamma = python::extract<double>(lattice.attr("gamma"));
    Lattice L(a, b, c, alpha, beta, gamma);
    this->setLattice(L);
    this->Molecule::setFromDiffPyStructure(stru);
}


// private class methods


boost::shared_ptr<Lattice> Crystal::getDefaultLattice()
{
    static boost::shared_ptr<Lattice> default_lattice;
    if (!default_lattice.get())
    {
        default_lattice.reset(new Lattice());
    }
    return default_lattice;
}


// private methods


void Crystal::init()
{
    this->_rmin = 0.0;
    this->_rmax = 0.0;
    this->_lattice = Crystal::getDefaultLattice();
    this->_lattice_max_ucd = this->_lattice->ucMaxDiagonalLength();
    this->_full_distance_table = this->Molecule::_distance_table;
    // the default _distance_table should exist and be blank
    assert(this->_full_distance_table.get());
    this->uncacheCostData();
    this->setDistReuse(true);
}


void Crystal::uncacheCostData() const
{
    this->_cost_data_cached = false;
    if (countAtoms() == 0)
    {
        this->ResetBadness();
        this->_count_pairs = 0;
        this->pmx_pair_counts.fill(0);
        this->_cost_data_cached = true;
    }
}


void Crystal::cropDistanceTable()
{
    DistanceTable::const_iterator lo, hi;
    lo = this->_full_distance_table->begin();
    hi = upper_bound(this->_full_distance_table->begin(),
            this->_full_distance_table->end(), this->getRmax());
    vector<double> vcropped(lo, hi);
    // Molecule::setDistanceTable must be passed DistanceTable
    // otherwise there would be infinite loop.
    DistanceTable dcropped(vcropped);
    this->Molecule::setDistanceTable(dcropped);
}


// position of random atom site offset by random lattice vector
const R3::Vector&
Crystal::anyOffsetAtomSite(const RandomWeighedGenerator& rwg) const
{
    assert(countAtoms() >= 1);
    static R3::Vector rv;
    int idx = rwg.weighedInt();
    R3::Vector mno(randomInt(2), randomInt(2), randomInt(2));
    rv = this->atoms[idx]->r + getLattice().cartesian(mno);
    return rv;
}


// shift cartesian vector to unit cell and adjust fractional
// coordinates that are very close to zero
R3::Vector Crystal::ucvCartesianAdjusted(const R3::Vector& cv) const
{
    const Lattice& lat = this->getLattice();
    R3::Vector ucl;
    ucl = lat.ucvFractional(lat.fractional(cv));
    R3::Vector ucs0, ucs1;
    for (int i = 0; i != R3::Ndim; ++i)
    {
        ucs0 = ucs1 = ucl;
        ucs0[i] = 0.0;
        ucs1[i] = 1.0;
        bool nearzero =
            lat.distance(ucl, ucs0) < NS_LIGA::eps_distance ||
            lat.distance(ucl, ucs1) < NS_LIGA::eps_distance;
        if (nearzero)   ucl[i] = 0.0;
    }
    R3::Vector res = lat.cartesian(ucl);
    return res;
}


// Move first atom to the origin of the lattice coordinate system.
void Crystal::shiftToOrigin()
{
    if (this->countAtoms() == 0)    return;
    // make sure the first atom is at the origin
    const R3::Vector rc0 = this->getAtom(0).r;
    if (R3::norm(rc0) > NS_LIGA::eps_distance)  this->Shift(-rc0);
}


// End of file
