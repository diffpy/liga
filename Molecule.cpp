/***********************************************************************
* Short Title: class Molecule - definitions
*
* Comments:
*
* $Id$
*
* <license text>
***********************************************************************/

#include <sstream>
#include <cassert>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include "Molecule.hpp"
#include "LigaUtils.hpp"
#include "AtomFilter_t.hpp"
#include "AtomSequence.hpp"
#include "AtomCost.hpp"
#include "Counter.hpp"
#include "Random.hpp"
#include "R3linalg.hpp"
#include "Exceptions.hpp"

using namespace std;
using namespace NS_LIGA;

////////////////////////////////////////////////////////////////////////
// Molecule definitions
////////////////////////////////////////////////////////////////////////

// static members
double Molecule::tol_nbad = 0.05*0.05;
double Molecule::tol_r = 1.0e-8;
double Molecule::promotefrac = 0.1;
bool Molecule::promotejump = true;
bool Molecule::promoterelax = false;
bool Molecule::demoterelax = false;
vector<AtomFilter_t*> Molecule::atom_filters;
double Molecule::lookout_prob = 0.0;
string Molecule::output_format = "rawxyz";

// class methods

long Molecule::getUniqueId()
{
    static long uid = 0;
    return uid++;
}

// constructors

Molecule::Molecule() : id(Molecule::getUniqueId())
{
    init();
    const static DistanceTable nodistances;
    setDistanceTable(nodistances);
}

Molecule::Molecule(const DistanceTable& dtab) : id(Molecule::getUniqueId())
{
    init();
    setDistanceTable(dtab);
}

Molecule::Molecule(const DistanceTable& dtab,
	const vector<double>& vx, const vector<double>& vy,
	const vector<double>& vz
	) : id(Molecule::getUniqueId())
{
    init();
    setDistanceTable(dtab);
    if (vx.size() != vy.size() || vx.size() != vz.size())
    {
        ostringstream emsg;
	emsg << "E: invalid coordinate vectors";
	throw InvalidMolecule(emsg.str());
    }
    for (size_t i = 0; i != vx.size(); ++i)
    {
        R3::Vector rc(vx[i], vy[i], vz[i]);
	AddAt(rc);
    }
}

Molecule::Molecule(const Molecule& M) : id(Molecule::getUniqueId())
{
    init();
    *this  = M;
}

Molecule& Molecule::operator=(const Molecule& M)
{
    if (this == &M) return *this;
    // Clear() must be the first statement
    Clear();
    // assign distance table - share reference when distances are reused
    if (M.getDistReuse())
    {
        this->_distance_table = M._distance_table;
    }
    else
    {
        this->setDistanceTable(*M._distance_table);
    }
    // duplicate source atoms
    atoms.resize(M.atoms.size(), NULL);
    for (AtomSequenceIndex seq(&M); !seq.finished(); seq.next())
    {
	atoms[seq.idx()] = new Atom_t(seq.ref());
    }
    pmx_used_distances = M.pmx_used_distances;
    pmx_partial_costs = M.pmx_partial_costs;
    free_pmx_slots = M.free_pmx_slots;
    // finished duplication
    this->_max_atom_count = M._max_atom_count;
    this->_badness = M._badness;
    this->_distreuse = M._distreuse;
    // IO helpers
    opened_file = M.opened_file;
    trace = M.trace;
    return *this;
}


Molecule* Molecule::clone() const
{
    Molecule* pclone = new Molecule(*this);
    return pclone;
}


void Molecule::init()
{
    this->_badness = 0;
    this->_max_atom_count = -1;
    this->_distreuse = false;
}

Molecule::~Molecule()
{
    // we must call Clear() to delete Atom_t objects
    Clear();
}

void Molecule::setDistanceTable(const DistanceTable& dtbl)
{
    this->_distance_table.reset(new DistanceTable(dtbl));
}

void Molecule::setDistanceTable(const vector<double>& dvec)
{
    DistanceTable dtbl(dvec);
    setDistanceTable(dtbl);
}

const DistanceTable& Molecule::getDistanceTable() const
{
    assert(this->_distance_table.get() != NULL);
    return *_distance_table;
}

void Molecule::setDistReuse(bool flag)
{
    this->_distreuse = flag;
}

bool Molecule::getDistReuse() const
{
    return this->_distreuse;
}

//////////////////////////////////////////////////////////////////////////
//// Molecule badness/fitness evaluation
//////////////////////////////////////////////////////////////////////////

double Molecule::maxTableDistance() const
{
    bool hasdistance =
        this->_distance_table.get() && !this->_distance_table->empty();
    double mx = hasdistance ? this->_distance_table->back() : 0.0;
    return mx;
}

void Molecule::recalculate() const
{
    if (countAtoms() > getMaxAtomCount())
    {
	ostringstream emsg;
	emsg << "E: molecule too large in recalculate()";
	throw InvalidMolecule(emsg.str());
    }
    // reset molecule
    this->ResetBadness();
    // reset all atoms
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
	seq.ptr()->ResetBadness();
    }
    // create array of all badness contributions with atom index
    for (AtomSequenceIndex seq0(this); !seq0.finished(); seq0.next())
    {
	AtomSequenceIndex seq1 = seq0;
	for (seq1.next(); !seq1.finished(); seq1.next())
	{
	    int i0 = seq0.ptr()->pmxidx;
	    int i1 = seq1.ptr()->pmxidx;
	    double& paircost = pmx_partial_costs(i0,i1);
            this->IncBadness(paircost);
            double paircosthalf = paircost / 2.0;
 	    seq0.ptr()->IncBadness(paircosthalf);
 	    seq1.ptr()->IncBadness(paircosthalf);
	}
    }
}

// Local helpers for Molecule::reassignPairs()

namespace {

class PairMatrixElement
{
    public:

        // class methods
        static bool compareDistance(
                const PairMatrixElement& p0, const PairMatrixElement& p1)
        {
            return p0.d01 < p1.d01;
        }

        // public data
        double d01;
        int i0;
        int i1;
};

}   // namespace

void Molecule::reassignPairs()
{
    if (getDistReuse())  return;
    vector<PairMatrixElement> pmx_elements;
    pmx_elements.reserve(countPairs());
    vector<double> used_distances;
    used_distances.reserve(countPairs());
    PairMatrixElement pme;
    for (AtomSequenceIndex seq0(this); !seq0.finished(); seq0.next())
    {
	AtomSequenceIndex seq1 = seq0;
        for (seq1.next(); !seq1.finished(); seq1.next())
        {
            pme.i0 = seq0.ptr()->pmxidx;
	    pme.i1 = seq1.ptr()->pmxidx;
            pme.d01 = R3::distance(seq0.ptr()->r, seq1.ptr()->r);
            pmx_elements.push_back(pme);
            assert(pmx_used_distances(pme.i0, pme.i1) >= 0.0);
            used_distances.push_back(pmx_used_distances(pme.i0, pme.i1));
        }
    }
    double orgbadness = Badness();
    sort(pmx_elements.begin(), pmx_elements.end(),
            PairMatrixElement::compareDistance);
    sort(used_distances.begin(), used_distances.end());
    vector<PairMatrixElement>::iterator pmii = pmx_elements.begin();
    vector<double>::iterator udii = used_distances.begin();
    assert(pmx_elements.size() == used_distances.size());
    for (; pmii != pmx_elements.end(); ++pmii, ++udii)
    {
        pmx_used_distances(pmii->i0, pmii->i1) = *udii;
    }
    recalculate();
    // increase orgbadness to avoid failures from round-off errors
    orgbadness = (1 + 1e-6)*orgbadness + 1e-6;
    assert(Badness() < orgbadness);
}

double Molecule::cost() const
{
    double rv;
    rv = countPairs() > 0 ? Badness()/countPairs() : 0.0;
    return rv;
}

const double& Molecule::Badness() const
{
    return this->_badness;
}

void Molecule::IncBadness(const double& db) const
{
    this->_badness += db;
}

void Molecule::DecBadness(const double& db) const
{
    this->_badness -= db;
}

void Molecule::ResetBadness(double b) const
{
    this->_badness = b;
}

bool Molecule::full() const
{
    bool rv = countAtoms() >= getMaxAtomCount();
    return rv;
}


int Molecule::countAtoms() const   
{
    return this->atoms.size();
}

int Molecule::countPairs() const
{
    int N = countAtoms();
    int rv = N * (N - 1) / 2;
    return rv;
}


bool comp_Atom_Badness(const Atom_t& lhs, const Atom_t& rhs)
{
    return lhs.Badness() < rhs.Badness();
}

bool comp_pAtom_Badness(const Atom_t* lhs, const Atom_t* rhs)
{
    return lhs->Badness() < rhs->Badness();
}


bool comp_pAtom_FreeBadness(const Atom_t* lhs, const Atom_t* rhs)
{
    return lhs->FreeBadness() < rhs->FreeBadness();
}


//////////////////////////////////////////////////////////////////////////
// Molecule operators
//////////////////////////////////////////////////////////////////////////

bool operator==(const Molecule& m0, const Molecule& m1)
{
    if (&m0 == &m1)
	return true;
    if (m0.getMaxAtomCount() != m1.getMaxAtomCount())
	return false;
    AtomSequence seq0(&m0), seq1(&m1);
    for (; !seq0.finished() && !seq1.finished() && *seq0.ptr() == *seq1.ptr();
	    seq0.next(), seq1.next() )
    { }
    return seq0.finished() && seq1.finished();
}

void Molecule::setMaxAtomCount(int cnt)
{
    if (cnt > this->_distance_table->estNumAtoms() && !getDistReuse())
    {
	ostringstream emsg;
	emsg << "E: not enough distances for max atom count = " <<
            cnt << '.' << "  Forgot to set distreuse?";
	throw InvalidMolecule(emsg.str());
    }
    else if (cnt < 1)
    {
	ostringstream emsg;
	emsg << "E: invalid value of max atom count = " << cnt;
	throw InvalidMolecule(emsg.str());
    }
    else if (cnt < countAtoms())
    {
	ostringstream emsg;
	emsg << "E: molecule too large in setMaxAtomCount()";
	throw InvalidMolecule(emsg.str());
    }
    _max_atom_count = cnt;
}

int Molecule::getMaxAtomCount() const
{
    int cnt;
    cnt = this->_max_atom_count < 0 ?
        this->getDistanceTable().estNumAtoms() :
        _max_atom_count;
    return cnt;
}

void Molecule::Shift(double dx, double dy, double dz)
{
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
	Atom_t* pa = seq.ptr();
	pa->r[0] += dx;
	pa->r[1] += dy;
	pa->r[2] += dz;
    }
}

void Molecule::Center()
{
    double avg_rx = 0.0, avg_ry = 0.0, avg_rz = 0.0;
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
	Atom_t* pa = seq.ptr();
	avg_rx += pa->r[0];
	avg_ry += pa->r[1];
	avg_rz += pa->r[2];
    }
    avg_rx /= countAtoms();
    avg_ry /= countAtoms();
    avg_rz /= countAtoms();
    Shift(-avg_rx, -avg_ry, -avg_rz);
}

void Molecule::Pop(const int aidx)
{
    if (aidx < 0 || aidx >= countAtoms())
    {
	ostringstream emsg;
	emsg << "Molecule::Pop(const int) invalid index " << aidx << '.';
	throw range_error(emsg.str());
    }
    Atom_t* pa = atoms[aidx];
    // Pop should never get called on fixed atom
    assert(!pa->fixed);
    removeAtomPairs(pa);
    free_pmx_slots.insert(pa->pmxidx);
    delete pa;
    atoms.erase(atoms.begin() + aidx);
}

void Molecule::Pop(const list<int>& cidx)
{
    // create a sorted set of indices of atoms to be popped
    set<int> popped(cidx.begin(), cidx.end());
    set<int>::reverse_iterator rii;
    for (rii = popped.rbegin(); rii != popped.rend(); ++rii)  Pop(*rii);
}

void Molecule::Clear()
{
    returnUsedDistances();
    // remove all atoms
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
	delete seq.ptr();
    }
    atoms.clear();
    free_pmx_slots.clear();
    ResetBadness();
}

void Molecule::Add(const Molecule& M)
{
    for (AtomSequence seq(&M); !seq.finished(); seq.next())
    {
	Add(seq.ref());
    }
}


void Molecule::AddAt(double rx0, double ry0, double rz0)
{
    Atom_t a(rx0, ry0, rz0);
    Add(a);
}


void Molecule::AddAt(const R3::Vector& rc)
{
    Atom_t a(rc);
    Add(a);
}


void Molecule::Add(const Atom_t& atom)
{
    if (countAtoms() == getMaxAtomCount())
    {
        ostringstream emsg;
	emsg << "E: molecule too large in Add()";
	throw InvalidMolecule(emsg.str());
    }
    // create new atom
    Atom_t* pnew_atom;
    pnew_atom = new Atom_t(atom);
    pnew_atom->ResetBadness();
    pnew_atom->pmxidx = getPairMatrixIndex();
    // create new pairs while summing up the costs
    addNewAtomPairs(pnew_atom);
    atoms.push_back(pnew_atom);
    if (full())     reassignPairs();
}


void Molecule::Fix(const int cidx)
{
    if (cidx < 0 || cidx >= countAtoms())
    {
	ostringstream emsg;
	emsg << "Molecule::Fix(const int) invalid index " << cidx << '.';
	throw range_error(emsg.str());
    }
    atoms[cidx]->fixed = true;
}

inline bool pAtom_is_fixed(const Atom_t* pa) { return pa->fixed; }

int Molecule::NFixed() const
{
    return count_if(atoms.begin(), atoms.end(), pAtom_is_fixed);
}

valarray<int> Molecule::good_neighbors_count(const vector<Atom_t>& vta)
{
    valarray<int> cnt(0, vta.size());
    // high limit for badness of a pair with good neighbor
    double hi_pbad = tol_nbad / 10.0;
    typedef vector<Atom_t>::const_iterator VAcit;
    VAcit a1i, a2i;
    int* pcnt1; int* pcnt2;
    for (  a1i = vta.begin(), pcnt1 = &(cnt[0]);
	    a1i != vta.end(); ++a1i, ++pcnt1 )
    {
	for (   a2i = a1i+1, pcnt2 = pcnt1+1; a2i != vta.end();
		++a2i, ++pcnt2 )
	{
	    double d = R3::distance(a1i->r, a2i->r);
	    double dd = *(this->_distance_table->find_nearest(d)) - d;
	    if (penalty(dd) < hi_pbad)
	    {
		++(*pcnt1);
		++(*pcnt2);
	    }
	}
    }
    return cnt;
}

void Molecule::filter_good_atoms(vector<Atom_t>& vta,
	double evolve_range, double hi_abad)
{
    if (countAtoms() == getMaxAtomCount())
    {
	ostringstream emsg;
	emsg << "E: Molecule too large in filter_good_atoms()";
	throw InvalidMolecule(emsg.str());
    }
    typedef vector<Atom_t>::iterator VAit;
    VAit gai;
    // first check if atoms pass through atom_filters
    gai = vta.begin();
    if ( !atom_filters.empty() )
    {
	for (VAit tai = vta.begin(); tai != vta.end(); ++tai)
	{
	    if (check_atom_filters(&(*tai)))
	    {
		*(gai++) = *tai;
	    }
	}
	vta.erase(gai, vta.end());
    }
    // obtain badness of every test atom, do lazy evaluation when badness
    // is higher than (minimum + evolve_range)
    AtomCost* atomcost = getAtomCostCalculator();
    atomcost->setCutoff(hi_abad);
    atomcost->setCutoffRange(evolve_range);
    for (VAit tai = vta.begin(); tai != vta.end(); ++tai)
    {
	tai->IncBadness(atomcost->eval(*tai));
    }
    // atom cost cutoff is available here,
    // let us keep only good atoms
    gai = vta.begin();
    for (VAit tai = vta.begin(); tai != vta.end(); ++tai)
    {
	if (tai->Badness() > atomcost->cutoff())    continue;
	*(gai++) = *tai;
    }
    vta.erase(gai, vta.end());
}

bool Molecule::check_atom_filters(Atom_t* pa)
{
    bool isgood = true;
    vector<AtomFilter_t*>::iterator pafi = atom_filters.begin();
    for (; isgood && pafi != atom_filters.end(); ++pafi)
    {
        AtomFilter_t* filter = *pafi;
	isgood = filter->Check(pa, this);
    }
    return isgood;
}

struct rxa_par
{
    typedef vector<Atom_t*> VPA;
    typedef valarray<double> VAD;
    VPA* atoms;
    VAD* ad0;
    VAD* wt;
};

int rxa_fdf(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
{
    static Atom_t ta(0.0, 0.0, 0.0);
    copy(x->data, x->data + x->size, ta.r.data());
    AtomCost* atomcost = static_cast<AtomCost*>(params);
    atomcost->eval(ta);
    for (size_t m = 0; m != atomcost->lsqComponentsSize(); ++m)
    {
	gsl_vector_set(f, m, atomcost->lsqComponents()[m]);
	for (size_t n = 0; n != atomcost->lsqParametersSize(); ++n)
	{
	    gsl_matrix_set(J, m, n, atomcost->lsqJacobianGet(m, n));
	}
    }
    return GSL_SUCCESS;
}

int rxa_f(const gsl_vector* x, void* params, gsl_vector* f)
{
    static Atom_t ta(0.0, 0.0, 0.0);
    copy(x->data, x->data + x->size, ta.r.data());
    AtomCost* atomcost = static_cast<AtomCost*>(params);
    atomcost->eval(ta);
    for (size_t m = 0; m != atomcost->lsqComponentsSize(); ++m)
    {
	gsl_vector_set(f, m, atomcost->lsqComponents()[m]);
    }
    return GSL_SUCCESS;
}

int rxa_df(const gsl_vector* x, void* params, gsl_matrix* J)
{
    static Atom_t ta(0.0, 0.0, 0.0);
    copy(x->data, x->data + x->size, ta.r.data());
    AtomCost* atomcost = static_cast<AtomCost*>(params);
    atomcost->eval(ta);
    for (size_t m = 0; m != atomcost->lsqComponentsSize(); ++m)
    {
	for (size_t n = 0; n != atomcost->lsqParametersSize(); ++n)
	{
	    gsl_matrix_set(J, m, n, atomcost->lsqJacobianGet(m, n));
	}
    }
    return GSL_SUCCESS;
}

void Molecule::RelaxAtom(const int cidx)
{
    if (cidx < 0 || cidx >= countAtoms())
    {
	ostringstream emsg;
	emsg << "Molecule::RelaxAtom(const int) invalid index " << cidx << '.';
	throw range_error(emsg.str());
    }
    // RelaxAtom should never get called on a fixed atom
    assert(!atoms[cidx]->fixed);
    Atom_t ta = getAtom(cidx);
    Pop(cidx);
    RelaxExternalAtom(ta);
    Add(ta);
}

void Molecule::RelaxExternalAtom(Atom_t& ta)
{
    const int maximum_relaxations = 20;
    const int maximum_iterations = 500;
    // pj: this seems to be crashing when countAtoms() < 3
    if (countAtoms() < 3)   return;
    // do relaxation on a copy of ta
    Atom_t rta(ta);
    // loop while badness is improved
    double lo_abad = DOUBLE_MAX;
    AtomCost* atomcost = getAtomCostCalculator();
    for (int nrelax = 0; nrelax < maximum_relaxations; ++nrelax)
    {
	double tbad = atomcost->eval(rta);
	// get out if tbad did not improve
	if (!eps_lt(tbad, lo_abad))	break;
	lo_abad = tbad;
	ta = rta;
	// get out if lo_abad is very low
	if (lo_abad < NS_LIGA::eps_cost)    break;
	// carry out relaxation otherwise
	// define function to be minimized
	gsl_multifit_function_fdf f;
	f.f = &rxa_f;
	f.df = &rxa_df;
	f.fdf = &rxa_fdf;
	f.n = atomcost->lsqComponentsSize();
	f.p = atomcost->lsqParametersSize();
	f.params = static_cast<void*>(atomcost);
	// bind rta coordinates to vector x
	gsl_vector_view x = gsl_vector_view_array(rta.r.data(), 3);
	// allocate solver
	gsl_multifit_fdfsolver* lms = gsl_multifit_fdfsolver_alloc(
		gsl_multifit_fdfsolver_lmsder, f.n, f.p);
	gsl_multifit_fdfsolver_set(lms, &f, &x.vector);
	// allocate vector for gradient test
	gsl_vector* G = gsl_vector_alloc(3);
	// minimize atom badness
	int iter = 0, status;
	do
	{
	    ++iter;
	    status = gsl_multifit_fdfsolver_iterate(lms);
	    if (status)
	    {
		if ( status != GSL_ETOLF && status != GSL_ETOLX &&
			status != GSL_CONTINUE )
		{
		    cerr << "LM solver status = " << status << " " <<
			gsl_strerror(status) << endl;
		}
		break;
	    }
	    // scale f_i with atom fitness
	    // pj: test gradient or do a simplex search?
	    gsl_multifit_gradient(lms->J, lms->f, G);
	    status = gsl_multifit_test_gradient(G, NS_LIGA::eps_cost/tol_r);
//	    status = gsl_multifit_test_delta(lms->dx, lms->x, tol_r, tol_r);
	}
	while (status == GSL_CONTINUE && iter < maximum_iterations);
	for (int i = 0; i < 3; ++i)
	    rta.r[i] = gsl_vector_get(lms->x, i);
	gsl_vector_free(G);
	gsl_multifit_fdfsolver_free(lms);
    }
}

AtomCost* Molecule::getAtomCostCalculator() const
{
    static AtomCost the_acc(this);
    the_acc.resetFor(this);
    return &the_acc;
}

void Molecule::addNewAtomPairs(Atom_t* pa)
{
    // calculate the cost
    AtomCost* atomcost = getAtomCostCalculator();
    atomcost->eval(pa);
    // store partial costs
    const vector<double>& ptcs = atomcost->partialCosts();
    assert(atoms.size() == ptcs.size());
    for (AtomSequenceIndex seq(this); !seq.finished(); seq.next())
    {
	double pairbadness = ptcs[seq.idx()];
	int idx0 = pa->pmxidx;
	int idx1 = seq.ptr()->pmxidx;
	pmx_partial_costs(idx0,idx1) = pairbadness;
	double badnesshalf = pairbadness / 2.0;
	seq.ptr()->IncBadness(badnesshalf);
	pa->IncBadness(badnesshalf);
    }
    this->IncBadness(atomcost->totalCost());
    if (this->Badness() < NS_LIGA::eps_cost)     this->ResetBadness();
    // remember used distances in pmx_used_distances
    if (!getDistReuse())
    {
        const vector<int>& didcs = atomcost->usedTargetDistanceIndices();
        const vector<int>& aidcs = atomcost->usedTargetAtomIndices();
        assert(didcs.size() == aidcs.size());
        vector<int>::const_iterator dii = didcs.begin();
        vector<int>::const_iterator aii = aidcs.begin();
        for (; dii != didcs.end(); ++dii, ++aii)
        {
            int idx0 = pa->pmxidx;
            int idx1 = atoms[*aii]->pmxidx;
            pmx_used_distances(idx0,idx1) = this->_distance_table->at(*dii);
        }
        // remove used distances from target table
        set<int> uidcs(didcs.begin(), didcs.end());
        set<int>::reverse_iterator ii;
        for (ii = uidcs.rbegin(); ii != uidcs.rend(); ++ii)
        {
            DistanceTable::iterator pos = this->_distance_table->begin() + *ii;
            this->_distance_table->erase(pos);
        }
    }
}

void Molecule::removeAtomPairs(Atom_t* pa)
{
    // remove associated pair costs
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
	if (pa == seq.ptr())	continue;
	int idx0 = pa->pmxidx;
	int idx1 = seq.ptr()->pmxidx;
	// remove pair costs
	double pairbadness = pmx_partial_costs(idx0,idx1);
	double badnesshalf = pairbadness/2.0;
	pa->DecBadness(badnesshalf);
	seq.ptr()->DecBadness(badnesshalf);
	this->DecBadness(pairbadness);
    }
    // return any associated used distances
    if (!getDistReuse())
    {
        for (AtomSequence seq(this); !seq.finished(); seq.next())
        {
            if (pa == seq.ptr())	continue;
            // return any used distances
            int idx0 = pa->pmxidx;
            int idx1 = seq.ptr()->pmxidx;
            double& udst = pmx_used_distances(idx0, idx1);
            if (udst > 0.0)     this->_distance_table->return_back(udst);
            udst = 0.0;
        }
    }
    if (this->Badness() < NS_LIGA::eps_cost)    this->ResetBadness();
}

int Molecule::push_good_distances(
        AtomArray& vta,
        const RandomWeighedGenerator& rwg,
        int ntrials)
{
    using NS_LIGA::eps_distance;
    if (!ntrials)   return 0;
    // add new atom in direction defined by 2 atoms
    if (countAtoms() == getMaxAtomCount())
    {
        ostringstream emsg;
	emsg << "E: molecule too large for finding a new position";
	throw InvalidMolecule(emsg.str());
    }
    else if (countAtoms() < 1)
    {
        ostringstream emsg;
	emsg << "E: empty molecule, no way to push_good_distances()";
	throw InvalidMolecule(emsg.str());
    }
    int push_count = 0;
    for (int nt = 0; nt < ntrials; ++nt)
    {
        const TriangulationAnchor& anch = this->getLineAnchor(rwg);
        R3::Vector direction(0.0, 0.0, 0.0);
        if (anch.count > 1)     direction = anch.B1 - anch.B0;
	// normalize direction if defined
	bool lattice_direction;
	double nm_direction = R3::norm(direction);
	if (nm_direction > eps_distance)
	{
	    direction /= nm_direction;
	    lattice_direction = true;
	}
	// otherwise orient along the z-axis
	else
	{
            direction = 0.0, 0.0, 1.0;
	    lattice_direction = false;
	}
	// pick free distance
        const DistanceTable& dtbl = *this->_distance_table;
	int didx = randomInt(dtbl.size());
	double radius = dtbl[didx];
	// add front atom
        R3::Vector nr;
        nr = anch.B0 + direction*radius;
	Atom_t ad1front(nr);
	ad1front.ttp = LINEAR;
	vta.push_back(ad1front);
	++push_count;
	// check opposite direction when it makes sense
	// this accounts for extra trial
	if (lattice_direction)
	{
	    ++nt;
            nr = anch.B0 - direction*radius;
	    Atom_t ad1back(nr);
	    ad1back.ttp = LINEAR;
	    vta.push_back(ad1back);
	    ++push_count;
	}
    }
    return push_count;
}

int Molecule::push_good_triangles(
	AtomArray& vta,
        const RandomWeighedGenerator& rwg,
        int ntrials)
{
    using NS_LIGA::eps_distance;
    if (!ntrials)   return 0;
    // generate randomly oriented triangles
    if (countAtoms() == getMaxAtomCount())
    {
        ostringstream emsg;
	emsg << "E: molecule too large for finding a new position";
	throw InvalidMolecule(emsg.str());
    }
    int push_count = 0;
    for (int nt = 0; nt < ntrials; ++nt)
    {
        const TriangulationAnchor& anch = getPlaneAnchor(rwg);
        const DistanceTable& dtbl = *this->_distance_table;
	// pick 2 vertex distances
        const PickType& didx = getDistReuse() ?
            randomPickWithRepeat(2, dtbl.size()) :
            randomPickFew(2, dtbl.size());
	double r02 = dtbl[didx[0]];
	double r12 = dtbl[didx[1]];
	double r01 = R3::distance(anch.B0, anch.B1);
	// is triangle base reasonably large?
	if (r01 < eps_distance)    continue;
	// get and store both possible values of xlong
	double xl0 = (r02*r02 + r01*r01 - r12*r12) / (2.0*r01);
	double xlong[2] = { xl0, r01-xl0 };
	// get and store both possible values of xperp
	double xp2 = r02*r02 - xlong[0]*xlong[0];
	double xp = sqrt(fabs(xp2));
	if (xp < eps_distance)  xp = 0.0;
	else if (xp2 < 0.0)  continue;
	double xperp[2] = { -xp, xp };
	// find direction along triangle base:
        R3::Vector longdir;
        longdir = (anch.B1 - anch.B0)/r01;
	// generate direction perpendicular to longdir
        R3::Vector perpdir(0.0, 0.0, 0.0);
	if (anch.count > 2)
	{
            perpdir = anch.B2 - anch.B0;
	    perpdir -= longdir*R3::dot(longdir, perpdir);
	}
	// normalize perpdir if defined
	bool lattice_plane;
	double nm_perpdir = R3::norm(perpdir);
	if (nm_perpdir > eps_distance)
	{
	    perpdir /= nm_perpdir;
	    lattice_plane = true;
	}
	// otherwise generate direction in cartesian plane
	// perpendicular to the smallest component of longdir
	else
	{
            R3::Vector uv(0.0, 0.0, 0.0);
            R3::Vector ald = fabs(longdir);
	    int ijk = minIndex(ald);
	    uv[ijk] = 1.0;
	    perpdir = R3::cross(longdir, uv);
	    perpdir /= R3::norm(perpdir);
	    lattice_plane = false;
	}
	// allocate vallarays for positions of a0 and vertex P
        R3::Vector P;
	// if vertex search has already failed above, nt would increase by 1
	// here we want nt to count number of added vertices
	--nt;
	// loops over all 4 vertices in case of lattice_plane
	for (double* pxl = xlong; pxl != xlong+2; ++pxl)
	{
	    for (double* pxp = xperp; pxp != xperp+2; ++pxp)
	    {
		++nt;
		P = anch.B0 + (*pxl)*longdir + (*pxp)*perpdir;
		Atom_t ad2(P);
		ad2.ttp = PLANAR;
		vta.push_back(ad2);
		++push_count;
		if (!lattice_plane)  break;
	    }
	    if (!lattice_plane)  break;
	}
    }
    return push_count;
}

int Molecule::push_good_pyramids(
        AtomArray& vta,
        const RandomWeighedGenerator& rwg,
        int ntrials)
{
    using NS_LIGA::eps_distance;
    if (!ntrials)   return 0;
    if (countAtoms() == getMaxAtomCount())
    {
        ostringstream emsg;
	emsg << "E: molecule too large for finding a new position";
	throw InvalidMolecule(emsg.str());
    }
    int push_count = 0;
    for (int nt = 0; nt < ntrials;)
    {
	// pick 3 base atoms
        const TriangulationAnchor& anch = getPyramidAnchor(rwg);
        const DistanceTable& dtbl = *this->_distance_table;
	// pick 3 vertex distances
	PickType didx = getDistReuse() ?
            randomPickWithRepeat(3, dtbl.size()) :
            randomPickFew(3, dtbl.size());
	// loop over all permutations of selected distances
	sort(didx.begin(), didx.end());
	do
	{
	    ++nt;
	    double r03 = dtbl[didx[0]];
	    double r13 = dtbl[didx[1]];
	    double r23 = dtbl[didx[2]];
	    // uvi is a unit vector in B0B1 direction
            R3::Vector uvi;
            uvi = anch.B1 - anch.B0;
	    double r01 = R3::norm(uvi);
	    if (r01 < eps_distance)    continue;
	    uvi /= r01;
	    // v02 is B0B2 vector
            R3::Vector v02;
            v02 = anch.B2 - anch.B0;
	    // uvj lies in B0B1B2 plane and is perpendicular to uvi
            R3::Vector uvj = v02;
	    uvj -= uvi*R3::dot(uvi, uvj);
	    double nm_uvj = R3::norm(uvj);
	    if (nm_uvj < eps_distance)  continue;
	    uvj /= nm_uvj;
	    // uvk is a unit vector perpendicular to B0B1B2 plane
            R3::Vector uvk = R3::cross(uvi, uvj);
	    double xP1 = -0.5/(r01)*(r01*r01 + r03*r03 - r13*r13);
	    // Pn are coordinates in pyramid coordinate system
            R3::Vector P1(xP1, 0.0, 0.0);
	    // vT is translation from pyramid to normal system
            R3::Vector vT;
            vT = anch.B0 - xP1*uvi;
	    // obtain coordinates of P3
            R3::Vector P3 = P1;
	    P3[0] += R3::dot(uvi, v02);
	    P3[1] += R3::dot(uvj, v02);
	    P3[2] = 0.0;
	    double xP3 = P3[0];
	    double yP3 = P3[1];
	    // find pyramid vertices
            R3::Vector P4;
	    double h2 = r03*r03 - xP1*xP1;
	    // does P4 belong to B0B1 line?
	    if (fabs(h2) < eps_distance)
	    {
		// is vertex on B0B1
		if (fabs(R3::norm(P3) - r03) > eps_distance)   continue;
		P4 = vT;
		Atom_t ad3(P4);
		ad3.ttp = SPATIAL;
		vta.push_back(ad3);
		++push_count;
		continue;
	    }
	    else if (h2 < 0)
	    {
		continue;
	    }
	    double yP4 = 0.5/(yP3)*(h2 + xP3*xP3 + yP3*yP3 - r23*r23);
	    double z2P4 = h2 - yP4*yP4;
	    // does P4 belong to B0B1B2 plane?
	    if (fabs(z2P4) < eps_distance)
	    {
		P4 = yP4*uvj + vT;
		Atom_t ad3(P4);
		ad3.ttp = SPATIAL;
		vta.push_back(ad3);
		++push_count;
		continue;
	    }
	    else if (z2P4 < 0)
	    {
		continue;
	    }
	    // here we can construct 2 pyramids
	    double zP4 = sqrt(z2P4);
	    // top one
	    P4 = yP4*uvj + zP4*uvk + vT;
	    Atom_t ad3top(P4);
	    ad3top.ttp = SPATIAL;
	    vta.push_back(ad3top);
	    ++push_count;
	    // and bottom one, which makes an extra trial
	    ++nt;
	    P4 = yP4*uvj - zP4*uvk + vT;
	    Atom_t ad3bottom(P4);
	    ad3bottom.ttp = SPATIAL;
	    vta.push_back(ad3bottom);
	    push_count++;
	} while ( next_permutation(didx.begin(), didx.end()) );
    }
    return push_count;
}


const Molecule::TriangulationAnchor&
Molecule::getLineAnchor(const RandomWeighedGenerator& rwg)
{
    assert(countAtoms() >= 1);
    static TriangulationAnchor anch;
    anch.count = min(countAtoms(),2);
    const PickType& aidx = rwg.weighedPick(anch.count);
    anch.B0 = this->atoms[aidx[0]]->r;
    if (anch.count > 1)
    {
        anch.B1 = this->atoms[aidx[1]]->r;
    }
    return anch;
}


const Molecule::TriangulationAnchor&
Molecule::getPlaneAnchor(const RandomWeighedGenerator& rwg)
{
    assert(countAtoms() >= 2);
    static TriangulationAnchor anch;
    // pick 2 atoms for base and 3rd for plane orientation
    anch.count = min(countAtoms(), 3);
    const PickType& aidx = rwg.weighedPick(anch.count);
    anch.B0 = this->atoms[aidx[0]]->r;
    anch.B1 = this->atoms[aidx[1]]->r;
    if (anch.count > 2)     anch.B2 = this->atoms[aidx[2]]->r;
    return anch;
}


const Molecule::TriangulationAnchor&
Molecule::getPyramidAnchor(const RandomWeighedGenerator& rwg)
{
    assert(countAtoms() >= 3);
    static TriangulationAnchor anch;
    // pick 3 base atoms
    const PickType& aidx = rwg.weighedPick(3);
    anch.count = 3;
    anch.B0 = this->atoms[aidx[0]]->r;
    anch.B1 = this->atoms[aidx[1]]->r;
    anch.B2 = this->atoms[aidx[2]]->r;
    return anch;
}


int Molecule::push_second_atoms(vector<Atom_t>& vta, int ntrials)
{
    if (countAtoms() != 1)
    {
        ostringstream emsg;
	emsg << "E: push_second_atoms() must be called with 1-atom molecule";
	throw InvalidMolecule(emsg.str());
    }
    // second atoms will be pushed along z-direction from a0
    Atom_t& a0 = *atoms[0];
    const DistanceTable& dtbl = *this->_distance_table;
    // new position
    R3::Vector nr = a0.r;
    int push_count = 0;
    if (ntrials > 2*dtbl.countUnique())
    {
	// we can push all the unique distances in both directions
	typedef vector<double>::iterator VDit;
	vector<double> dtu = dtbl.unique();
	for (VDit dui = dtu.begin(); dui != dtu.end(); ++dui)
	{
	    // top atom
	    nr[2] = a0.r[2] + *dui;
	    Atom_t adtop(nr); adtop.ttp = LINEAR;
	    vta.push_back(adtop);
	    // bottom atom
	    nr[2] = a0.r[2] - *dui;
	    Atom_t adbot(nr); adbot.ttp = LINEAR;
	    vta.push_back(adbot);
	    push_count += 2;
	}
    }
    else
    {
	// distances will be picked randomly
	for (push_count = 0; push_count < ntrials; ++push_count)
	{
	    int didx = randomInt(dtbl.size());
	    double dz = dtbl[didx];
	    // randomize direction
            dz *= plusminus();
	    nr[2] = a0.r[2] + dz;
	    Atom_t ad(nr); ad.ttp = LINEAR;
	    vta.push_back(ad);
	}
    }
    return push_count;
}

int Molecule::push_third_atoms(vector<Atom_t>& vta, int ntrials)
{
    using NS_LIGA::eps_distance;
    if (countAtoms() != 2)
    {
        ostringstream emsg;
	emsg << "E: push_third_atoms() must be called with 2-atom molecule";
	throw InvalidMolecule(emsg.str());
    }
    const DistanceTable& dtbl = *this->_distance_table;
    // build list of distances from the base atoms
    list<double> d0, d1;
    if (ntrials > 2 * dtbl.countUnique() * dtbl.countUnique())
    {
	// we can push all the unique triangles
	typedef vector<double>::iterator VDit;
	vector<double> dtu = dtbl.unique();
	for (VDit ui0 = dtu.begin(); ui0 != dtu.end(); ++ui0)
	{
	    for (VDit ui1 = dtu.begin(); ui1 != dtu.end(); ++ui1)
	    {
		d0.push_back(*ui0);
		d1.push_back(*ui1);
	    }
	}
    }
    else
    {
	// distances will be picked randomly
	for (int i = 0; i < ntrials; ++i)
	{
            const PickType& didx = getDistReuse() ?
                randomPickWithRepeat(2, dtbl.size()) :
                randomPickFew(2, dtbl.size());
	    d0.push_back( dtbl[didx[0]] );
	    d1.push_back( dtbl[didx[1]] );
	}
    }
    // define some constants and base atoms
    int push_count = 0;
    Atom_t a0 = *atoms[0];
    Atom_t a1 = *atoms[1];
    double r01 = R3::distance(a0.r, a1.r);
    // longitudinal direction along triangle base
    R3::Vector longdir;
    longdir = (a1.r - a0.r)/r01;
    // perpendicular direction to longdir is a cross product of x with longdir
    R3::Vector ux(0.0, 0.0, 0.0);
    ux[0] = 1.0;
    R3::Vector perpdir = R3::cross(ux, longdir);
    double nm_perpdir = R3::norm(perpdir);
    if (nm_perpdir == 0.0)
    {
	// longdir is ux, let us set perpdir = uy
	perpdir = 0.0;
	perpdir[1] = 1.0;
    }
    else
    {
	perpdir /= nm_perpdir;
    }
    R3::Vector Pa0 = a0.r;
    // now loop over list of distances
    typedef list<double>::iterator LDit;
    for (   LDit d0i = d0.begin(), d1i = d1.begin();
	    d0i != d0.end() && d1i != d1.end(); ++d0i, ++d1i )
    {
	double& r02 = *d0i;
	double& r12 = *d1i;
	// is triangle base reasonably large?
	if (r01 < eps_distance)     continue;
	// calculate xlong
	double xlong = (r02*r02 + r01*r01 - r12*r12) / (2.0*r01);
	// calculate xperp
	double xp2 = r02*r02 - xlong*xlong;
	double xperp = sqrt(fabs(xp2));
	if (xperp < eps_distance)   xperp = 0.0;
	else if (xp2 < 0.0)         continue;
	else if (randomInt(2) == 0) xperp = -xperp;
	// add atom
        R3::Vector Pn;
	Pn = Pa0 + xlong*longdir + xperp*perpdir;
	Atom_t ad(Pn); ad.ttp = PLANAR;
	vta.push_back(ad);
	++push_count;
    }
    return push_count;
}

const pair<int*,int*>& Molecule::Evolve(const int* est_triang)
{
    // aliases for input arguments
    const int& nlinear = est_triang[LINEAR];
    const int& nplanar = est_triang[PLANAR];
    const int& nspatial = est_triang[SPATIAL];
    // static result arrays
    static int acc[NTGTYPES];
    static int tot[NTGTYPES];
    static pair<int*,int*> acc_tot(acc, tot);
    // reset result arrays
    fill(acc, acc + NTGTYPES, 0);
    fill(tot, tot + NTGTYPES, 0);
    if (countAtoms() == getMaxAtomCount())
    {
        ostringstream emsg;
	emsg << "E: full-sized molecule cannot Evolve()";
	throw InvalidMolecule(emsg.str());
    }
    // containter for test atoms
    vector<Atom_t> vta;
    bool lookout = lookout_prob && 0 < countAtoms() && countAtoms() <= 2 &&
	randomFloat() < lookout_prob;
    const int lookout_trials = 1500;
    // evolution is trivial for empty or 1-atom molecule
    switch (countAtoms())
    {
	case 0:
	    AddAt(0.0, 0.0, 0.0);
            acc[LINEAR] = 1;
            tot[LINEAR] = 1;
	    return acc_tot;
	case 1:
	    if (lookout)
	    {
		push_second_atoms(vta, lookout_trials);
		break;
	    }
	case 2:
	    if (lookout)
	    {
		push_third_atoms(vta, lookout_trials);
		break;
	    }
	default:
            // create random generator weighed with atom fitnesses
            // calculate array of atom fitnesses
            vector<double> vacost(countAtoms());
            vector<double> vafit(countAtoms());
            // first fill the cost array
            for (AtomSequenceIndex seq(this); !seq.finished(); seq.next())
            {
                vacost[seq.idx()] = seq.ptr()->Badness();
            }
            // then convert to fitness
            vafit = costToFitness(vacost);
            // finally create RandomWeighedGenerator
            RandomWeighedGenerator rwg(vafit.begin(), vafit.end());
	    push_good_distances(vta, rwg, nlinear);
	    push_good_triangles(vta, rwg, nplanar);
	    push_good_pyramids(vta, rwg, nspatial);
    }
    typedef vector<Atom_t>::iterator VAit;
    // count total triangulation attempts
    for (VAit ai = vta.begin(); ai != vta.end(); ++ai)
    {
        ++tot[ai->ttp];
    }
    // set badness range from desired badness
    double evolve_range = countAtoms()*tol_nbad*promotefrac;
    double hi_abad = DOUBLE_MAX;
    // try to add as many atoms as possible
    while (true)
    {
	filter_good_atoms(vta, evolve_range, hi_abad);
        // finished when no test atoms left
	if (vta.empty())   break;
	// calculate fitness of test atoms
	valarray<double> vtafit(vta.size());
	if (lookout)
	{
	    // calculate fitness as the number of good neighbors
	    valarray<int> cnt = good_neighbors_count(vta);
	    int max_cnt = cnt.max();
	    double* pfit = &vtafit[0];
	    int* pcnt = &cnt[0];
	    for ( ; pfit != &vtafit[vtafit.size()]; ++pfit, ++pcnt)
	    {
		*pfit = (*pcnt < max_cnt/2) ? 0.0 : *pcnt;
	    }
	}
	else
	{
	    // calculate fitness as reciprocal value of badness
	    // fill the vtafit array with badness
	    double* pfit = &vtafit[0];
	    for (VAit ai = vta.begin(); ai != vta.end(); ++ai, ++pfit)
	    {
		*pfit = ai->Badness();
	    }
	    // then get the reciprocal value
            double* ftnfirst = &(vtafit[0]);
            double* ftnlast = &(vtafit[vtafit.size()]);
	    transform(ftnfirst, ftnlast, ftnfirst, convertCostToFitness);
	}
	// vtafit is ready here
	int idx = randomWeighedInt(vtafit.size(), &vtafit[0]);
	Add(vta[idx]);
        acc[vta[idx].ttp]++;
	hi_abad = vta[idx].Badness() + evolve_range;
	vta.erase(vta.begin()+idx);
	if (promoterelax)
	{
	    int worst_idx = max_element(atoms.begin(), atoms.end(),
		    comp_pAtom_FreeBadness) - atoms.begin();
            Atom_t* worst = atoms[worst_idx];
	    if (eps_gt(worst->Badness(), 0.0) && !worst->fixed)
	    {
		RelaxAtom(worst_idx);
	    }
	}
	if (full() || !promotejump)     break;
	for (VAit pai = vta.begin(); pai != vta.end(); ++pai)
        {
	    pai->ResetBadness();
        }
    }
    return acc_tot;
}

void Molecule::Degenerate(int Npop)
{
    Npop = min(countAtoms(), Npop);
    if (Npop == 0)  return;
    // build array of atom badnesses
    double freebad[countAtoms()];
    int freeidx[countAtoms()];
    int Nfree = 0;
    for (int i = 0; i != countAtoms(); ++i)
    {
	Atom_t* pai = atoms[i];
	if ( pai->fixed )  continue;
	freebad[Nfree] = pai->Badness();
	freeidx[Nfree] = i;
	Nfree++;
    }
    if (Nfree == 0)  return;
    Npop = min(Nfree, Npop);
    // build list of indices of atoms to pop
    const PickType& idxidx = randomWeighedPick(Npop, Nfree, freebad);
    list<int> ipop;
    PickType::const_iterator ii;
    for (ii = idxidx.begin(); ii != idxidx.end(); ++ii)
    {
	ipop.push_back(freeidx[*ii]);
    }
    Pop(ipop);
    if (demoterelax && countAtoms() > 1)
    {
        int worst_idx = max_element(atoms.begin(), atoms.end(),
                comp_pAtom_FreeBadness) - atoms.begin();
        Atom_t* worst = atoms[worst_idx];
        if (eps_gt(worst->Badness(), 0.0) && !worst->fixed)
        {
            RelaxAtom(worst_idx);
        }
    }
}

int Molecule::getPairMatrixIndex()
{
    int idx;
    if (!free_pmx_slots.empty())
    {
	set<int>::iterator firstfree = free_pmx_slots.begin();
	idx = *firstfree;
	free_pmx_slots.erase(firstfree);
    }
    else
    {
        // no free slots here
        idx = atoms.size();
        resizePairMatrices(idx + 1);
    }
    assert(idx < (int) this->pmx_partial_costs.rows());
    return idx;
}

void Molecule::resizePairMatrices(int sz)
{
    int szcur = this->pmx_partial_costs.rows();
    if (sz <= szcur)    return;
    int sznew1 = sz;
    int sznew2 = min(2*szcur, getMaxAtomCount());
    int sznew = max(sznew1, sznew2);
    assert(sznew <= getMaxAtomCount());
    this->pmx_partial_costs.resize(sznew, sznew, 0.0);
    if (!getDistReuse())
    {
        this->pmx_used_distances.resize(sznew, sznew, 0.0);
    }
}


void Molecule::returnUsedDistances()
{
    if (getDistReuse() || !this->_distance_table.get())
    {
        assert(this->pmx_used_distances.rows() == 0);
        return;
    }
    // return used distances
    for (AtomSequence seq0(this); !seq0.finished(); seq0.next())
    {
	AtomSequence seq1 = seq0;
	for (seq1.next(); !seq1.finished(); seq1.next())
	{
	    int i0 = seq0.ptr()->pmxidx;
	    int i1 = seq1.ptr()->pmxidx;
	    double& udst = pmx_used_distances(i0, i1);
	    if (udst > 0.0)	this->_distance_table->push_back(udst);
	    udst = 0.0;
	}
    }
    sort(this->_distance_table->begin(), this->_distance_table->end());
}


////////////////////////////////////////////////////////////////////////
// Molecule IO functions
////////////////////////////////////////////////////////////////////////

Molecule::ParseHeader::ParseHeader(const string& s) : header(s)
{
    // parse format
    string fmt;
    // initialize members:
    state =
	read_token("LIGA molecule format", fmt) &&
	read_token("NAtoms", NAtoms);
    if (!state)
    {
	if (read_token("Number of particles", NAtoms))
	{
	    format = ATOMEYE;
	    fmt = "atomeye";
	}
	else
	    return;
    }
    else if (fmt == "xyz")
	format = XYZ;
    else if (fmt == "atomeye")
	format = ATOMEYE;
    else
    {
	state = false;
	return;
    }
}

template<typename T> bool Molecule::ParseHeader::read_token(
	const char* token, T& value
	)
{
    const char* fieldsep = ":= ";
    int ltoken = strlen(token);
    string::size_type sp;
    const string::size_type npos = string::npos;
    if (
	    npos == (sp = header.find(token)) ||
	    npos == (sp = header.find_first_not_of(fieldsep, sp+ltoken))
       )
    {
	return false;
    }
    istringstream istrs(header.substr(sp));
    bool result = (istrs >> value);
    return result;
}

istream& Molecule::ReadXYZ(istream& fid)
{
    // read values to integer vector vxyz
    string header;
    vector<double> vxyz;
    bool result = read_header(fid, header) && read_data(fid, vxyz);
    if (!result) return fid;
    int vxyz_NAtoms = vxyz.size()/3;
    // check how many numbers were read
    ParseHeader ph(header);
    if (ph && vxyz_NAtoms != ph.NAtoms)
    {
        ostringstream emsg;
	emsg << "E: " << opened_file << ": expected " << ph.NAtoms <<
	    " atoms, read " << vxyz_NAtoms;
	throw IOError(emsg.str());
    }
    if ( vxyz.size() % 3 )
    {
        ostringstream emsg;
	emsg << "E: " << opened_file << ": incomplete data";
	throw IOError(emsg.str());
    }
    Clear();
    for (size_t i = 0; i + 2 < vxyz.size(); i += 3)
    {
	AddAt(vxyz[i], vxyz[i + 1], vxyz[i + 2]);
    }
    return fid;
}

bool Molecule::ReadXYZ(const char* file)
{
    // open file for reading
    ifstream fid(file);
    if (!fid)
    {
        ostringstream emsg;
	emsg << "E: unable to read '" << file << "'";
	throw IOError(emsg.str());
    }
    opened_file = file;
    bool result = ReadXYZ(fid);
    opened_file.clear();
    fid.close();
    return result;
}

void Molecule::WriteFile(const char* filename)
{
    // invalid structure format can throw exception,
    // check if write operator works first
    ostringstream output;
    output << *this;
    // test if filename is writeable
    ofstream fid(filename, ios_base::out|ios_base::ate);
    if (!fid)
    {
	ostringstream emsg;
	emsg << "WriteFile(): unable to write to '" << filename << "'";
	throw IOError(emsg.str());
    }
    fid.close();
    // write via temporary file
    int filelen = strlen(filename);
    char writefile[filelen+6+1];
    strcpy(writefile, filename);
    memset(writefile+filelen, 'X', 6);
    writefile[filelen+6] = '\0';
    mktempofstream(fid, writefile);
    fid << output.str();
    fid.close();
    rename(writefile, filename);
}

void Molecule::setOutputFormat(const std::string& format)
{
    Molecule::output_format = format;
}

istream& operator>>(istream& fid, Molecule& M)
{
    string header;
    istream::pos_type p = fid.tellg();
    if( !read_header(fid, header) )
    {
	fid.setstate(ios_base::failbit);
	return fid;
    }
    fid.seekg(p);
    Molecule::ParseHeader ph(header);
    if (!ph)
    {
	fid.setstate(ios_base::failbit);
	return fid;
    }
    bool result = true;
    switch (ph.format)
    {
	case Molecule::XYZ:
	    result = M.ReadXYZ(fid);
	    break;
	case Molecule::ATOMEYE:
	    throw runtime_error("reading of atomeye files not implemented");
	    break;
    }
    if (!result)
    {
	fid.setstate(ios_base::failbit);
    }
    return fid;
}

ostream& operator<<(ostream& fid, Molecule& M)
{
    namespace python = boost::python;
    try {
        initializePython();
        string element;
        element = (Molecule::output_format == "rawxyz") ? "" : "C";
        python::object stru = M.newDiffPyStructure();
        typedef vector<Atom_t*>::iterator VPAit;
        for (VPAit pai = M.atoms.begin(); pai != M.atoms.end(); ++pai)
        {
            Atom_t& ai = **pai;
            stru.attr("addNewAtom")(element);
            python::object alast = stru.attr("getLastAtom")();
            python::object xyz_cartn;
            xyz_cartn = python::make_tuple(ai.r[0], ai.r[1], ai.r[2]);
            alast.attr("xyz_cartn") = xyz_cartn;
        }
        string s;
        s = python::call_method<string>(stru.ptr(),
                "writeStr", Molecule::output_format);
        fid << s;
    }
    catch (python::error_already_set) {
        if (PyErr_Occurred())
        {
            PyErr_Print();
        }
        const char* emsg = "Cannot output structure.";
        throw IOError(emsg);
    }
    return fid;
}

boost::python::object Molecule::newDiffPyStructure()
{
    namespace python = boost::python;
    python::object mstru = python::import("diffpy.Structure");
    python::object stru = mstru.attr("Structure")();
    return stru;
}

void Molecule::PrintBadness() const
{
    if (!countAtoms())  return;
    cout << "ABadness() =";
    double mab = (*max_element(atoms.begin(), atoms.end(),
		    comp_pAtom_Badness)) -> Badness();
    bool marked = false;
    vector<Atom_t*>::const_iterator pai;
    for (pai = atoms.begin(); pai != atoms.end(); ++pai)
    {
	cout << ' ';
	if (!marked && (*pai)->Badness() == mab)
	{
	    cout << '+';
	    marked = true;
	}
	cout << (*pai)->Badness();
    }
    cout << endl;
}

void Molecule::PrintFitness()
{
    valarray<double> vacost(countAtoms());
    valarray<double> vafit(countAtoms());
    // first fill the array with badness
    double* pd = &vacost[0];
    vector<Atom_t*>::iterator pai = atoms.begin();
    for (; pai != atoms.end(); ++pai, ++pd)     *pd = (*pai)->Badness();
    // then get the reciprocal value
    vafit = costToFitness(vacost);
    cout << "AFitness() =";
    double mab = vafit.max();
    bool marked = false;
    for (pd = &vafit[0]; pd != &vafit[vafit.size()]; ++pd)
    {
	cout << ' ';
	if (!marked && *pd == mab)
	{
	    cout << '+';
	    marked = true;
	}
	cout << *pd;
    }
    cout << endl;
}

// End of file
