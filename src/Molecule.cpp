/*****************************************************************************
* Short Title: class Molecule - definitions
*
* Comments:
*
* <license text>
*****************************************************************************/

#include <sstream>
#include <cassert>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit_nlin.h>
#include <boost/foreach.hpp>

#include "Molecule.hpp"
#include "LigaUtils.hpp"
#include "AtomFilter_t.hpp"
#include "AtomSequence.hpp"
#include "AtomCost.hpp"
#include "AtomOverlap.hpp"
#include "Counter.hpp"
#include "Random.hpp"
#include "R3linalg.hpp"
#include "Exceptions.hpp"

using namespace std;
using namespace NS_LIGA;

//////////////////////////////////////////////////////////////////////////////
// class Molecule
//////////////////////////////////////////////////////////////////////////////

// Static Data ---------------------------------------------------------------

double Molecule::tol_nbad = 0.05*0.05;
double Molecule::tol_r = 1.0e-8;
double Molecule::promotefrac = 0.1;
bool Molecule::promotejump = true;
bool Molecule::promoterelax = false;
bool Molecule::demoterelax = false;
vector<AtomFilter_t*> Molecule::atom_filters;
string Molecule::output_format = "rawxyz";

// Class Methods -------------------------------------------------------------

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
    this->CheckIntegrity();
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
    this->_atom_radii_table = M._atom_radii_table;
    // duplicate source atoms
    atoms_storage = M.atoms_storage;
    // build dictionary to find equivalent Atom_t pointers
    map<const Atom_t*,Atom_t*> mapatomptr;
    list<Atom_t>::const_iterator asrc = M.atoms_storage.begin();
    list<Atom_t>::iterator adst = atoms_storage.begin();
    for (; adst != atoms_storage.end(); ++asrc, ++adst)
    {
        mapatomptr[&(*asrc)] = &(*adst);
    }
    // translate pointers to the atoms present in the molecule
    atoms.resize(M.atoms.size());
    for (size_t i = 0; i != M.atoms.size(); ++i)
    {
        atoms[i] = mapatomptr[M.atoms[i]];
        assert(atoms[i] != NULL);
    }
    // translate pointers to the unassigned atoms in the bucket
    atoms_bucket.resize(M.atoms_bucket.size());
    for (size_t i = 0; i != M.atoms_bucket.size(); ++i)
    {
        assert(mapatomptr.count(M.atoms_bucket[i]));
        atoms_bucket[i] = mapatomptr[M.atoms_bucket[i]];
        assert(atoms_bucket[i] != NULL);
    }
    pmx_used_distances = M.pmx_used_distances;
    pmx_partial_costs = M.pmx_partial_costs;
    free_pmx_slots = M.free_pmx_slots;
    // finished duplication
    this->_badness = M._badness;
    this->_overlap = M._overlap;
    this->_distreuse = M._distreuse;
    this->_samepairradius = M._samepairradius;
    // IO helpers
    trace = M.trace;
    return *this;
}


Molecule* Molecule::copy() const
{
    Molecule* pclone = new Molecule(*this);
    return pclone;
}


void Molecule::init()
{
    static boost::shared_ptr<DistanceTable>
        empty_distance_table(new DistanceTable());
    this->_distance_table = empty_distance_table;
    this->_badness = 0.0;
    this->_overlap = 0.0;
    this->_distreuse = false;
    this->_samepairradius = -1.0;
}


void Molecule::setDistanceTable(const DistanceTable& dtbl)
{
    this->_distance_table.reset(new DistanceTable(dtbl));
    if (atoms_storage.empty() && !dtbl.empty())
    {
        ChemicalFormula chfm;
        ChemicalFormula::value_type sc("", dtbl.estNumAtoms());
        chfm.push_back(sc);
        this->setChemicalFormula(chfm);
    }
}


void Molecule::setDistanceTable(const vector<double>& dvec)
{
    DistanceTable dtbl(dvec);
    setDistanceTable(dtbl);
}


const DistanceTable& Molecule::getDistanceTable() const
{
    assert(this->_distance_table.get() != NULL);
    return *(this->_distance_table);
}


DistanceTable Molecule::getDistanceTable()
{
    assert(this->_distance_table.get() != NULL);
    return *(this->_distance_table);
}


void Molecule::setDistReuse(bool flag)
{
    this->_distreuse = flag;
}


bool Molecule::getDistReuse() const
{
    return this->_distreuse;
}


double Molecule::getMaxAtomRadius() const
{
    double maxradius = 0.0;
    list<Atom_t>::const_iterator ai;
    for (ai = atoms_storage.begin(); ai != atoms_storage.end(); ++ai)
    {
        if (ai->radius > maxradius)  maxradius = ai->radius;
    }
    return maxradius;
}


// Molecule Badness/Fitness Evaluation ---------------------------------------

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
            const double& paircost = pmx_partial_costs(i0,i1);
            this->IncBadness(paircost);
            double paircosthalf = paircost / 2.0;
            seq0.ptr()->IncBadness(paircosthalf);
            seq1.ptr()->IncBadness(paircosthalf);
        }
    }
    this->recalculateOverlap();
}


AtomCost* Molecule::getAtomCostCalculator() const
{
    static AtomCost the_acc(this);
    the_acc.resetFor(this);
    return &the_acc;
}


AtomCost* Molecule::getAtomOverlapCalculator() const
{
    static AtomOverlap the_overlap_calculator(this);
    the_overlap_calculator.resetFor(this);
    return &the_overlap_calculator;
}


void Molecule::setAtomCostScale(double sc)
{
    this->getAtomCostCalculator()->setScale(sc);
    this->recalculate();
}


void Molecule::setAtomOverlapScale(double sc)
{
    this->getAtomOverlapCalculator()->setScale(sc);
    this->recalculate();
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
    const AtomCost* atomcost = this->getAtomCostCalculator();
    const DistanceTable& dtgt = this->getDistanceTable();
    for (; pmii != pmx_elements.end(); ++pmii, ++udii)
    {
        pmx_used_distances(pmii->i0, pmii->i1) = *udii;
        double dd = *udii - pmii->d01;
        double desd = dtgt.getesd(*udii);
        pmx_partial_costs(pmii->i0, pmii->i1) =
            atomcost->penaltyScaled(dd, desd);
    }
    recalculate();
    // increase orgbadness to avoid failures from round-off errors
    orgbadness = (1 + 1e-6)*orgbadness + 1e-6;
    assert(Badness() < orgbadness);
}


double Molecule::cost() const
{
    double rv = this->costDistance() + this->costOverlap();
    return rv;
}


double Molecule::costDistance() const
{
    double rv;
    rv = countPairs() > 0 ? this->Badness() / this->countPairs() : 0.0;
    return rv;
}


double Molecule::costOverlap() const
{
    double rv;
    rv = countAtoms() > 0 ? this->Overlap() / this->countAtoms() : 0.0;
    return rv;
}


const double& Molecule::Badness() const
{
    return this->_badness;
}


void Molecule::IncBadness(const double& db) const
{
    this->_badness += db;
    // Take care of round-offs, but only if they are very small.
    if (db < 0.0 && isNearZeroRoundOff(this->_badness))
    {
        this->ResetBadness(0.0);
    }
}


void Molecule::DecBadness(const double& db) const
{
    this->IncBadness(-db);
}


void Molecule::ResetBadness(double b) const
{
    this->_badness = b;
}


const double& Molecule::Overlap() const
{
    return this->_overlap;
}


void Molecule::IncOverlap(const double& doverlap) const
{
    this->_overlap += doverlap;
    // Take care of round-offs, but only if they are very small.
    if (doverlap < 0.0 && isNearZeroRoundOff(this->_overlap))
    {
        this->ResetOverlap(0.0);
    }
}


void Molecule::DecOverlap(const double& doverlap) const
{
    this->IncOverlap(-doverlap);
}


void Molecule::ResetOverlap(double overlap) const
{
    this->_overlap = overlap;
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


double Molecule::pairsPerAtom() const
{
    double ppa = countAtoms() ? 1.0 * countPairs() / countAtoms() : 0.0;
    return ppa;
}


double Molecule::pairsPerAtomInc() const
{
    int n = countAtoms();
    double m = n ? (countPairs() / (1.0 * n * n) + 0.5 + 0.5 / n) : 1.0;
    double ppa = (n + 1) * (m - 0.5) - 0.5;
    return ppa;
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


bool comp_pAtom_FreeOverlap(const Atom_t* lhs, const Atom_t* rhs)
{
    return lhs->FreeOverlap() < rhs->FreeOverlap();
}


// Molecule Operators --------------------------------------------------------

void Molecule::setChemicalFormula(const string& s)
{
    ChemicalFormula formula(s);
    this->setChemicalFormula(formula);
}


void Molecule::setChemicalFormula(const ChemicalFormula& formula)
{
    ChemicalFormula::const_iterator ec;
    vector<string> fmexpanded = formula.expand();
    atoms_storage.resize(formula.countElements(), Atom_t("", 0.0, 0.0, 0.0));
    set<const Atom_t*> set_storage;
    set<const Atom_t*> set_existing;
    set_existing.insert(atoms_bucket.begin(), atoms_bucket.end());
    set_existing.insert(atoms.begin(), atoms.end());
    vector<Atom_t*> replaced;
    BOOST_FOREACH (Atom_t& a, atoms_storage)
    {
        set_storage.insert(&a);
        // add any new atoms to the bucket
        if (!set_existing.count(&a))  atoms_bucket.push_back(&a);
        // apply formula, store pointers to atoms that need to be updated
        vector<string>::iterator ei;
        ei = find(fmexpanded.begin(), fmexpanded.end(), a.element);
        if (ei == fmexpanded.end())  replaced.push_back(&a);
        else  fmexpanded.erase(ei);
    }
    assert(fmexpanded.size() == replaced.size());
    BOOST_FOREACH (Atom_t* pa, replaced)
    {
        pa->element = fmexpanded.front();
        fmexpanded.erase(fmexpanded.begin());
    }
    // remove any atom pointers that are not it the storage
    vector<Atom_t*>::iterator ai, ag;
    for (ai = ag = atoms_bucket.begin(); ai != atoms_bucket.end(); ++ai)
    {
        if (set_storage.count(*ai))   (*ag++) = *ai;
    }
    atoms_bucket.erase(ag, atoms_bucket.end());
    for (ai = ag = atoms.begin(); ai != atoms.end(); ++ai)
    {
        if (set_storage.count(*ai))   (*ag++) = *ai;
    }
    atoms.erase(ag, atoms.end());
    this->fetchAtomRadii();
    this->recalculate();
    this->CheckIntegrity();
}


ChemicalFormula Molecule::getChemicalFormula() const
{
    ChemicalFormula rv;
    BOOST_FOREACH (const Atom_t& a, atoms_storage)
    {
        if (rv.empty() || rv.back().first != a.element)
        {
            ChemicalFormula::value_type elcnt(a.element, 0);
            rv.push_back(elcnt);
        }
        rv.back().second += 1;
    }
    return rv;
}


void Molecule::setAtomRadiiTable(const string& s)
{
    AtomRadiiTable table;
    table.fromString(s);
    this->setAtomRadiiTable(table);
}


void Molecule::setAtomRadiiTable(const AtomRadiiTable& radiitable)
{
    if (&(this->getAtomRadiiTable()) == &radiitable)  return;
    this->_atom_radii_table.reset(new AtomRadiiTable(radiitable));
    this->fetchAtomRadii();
    this->recalculateOverlap();
}


const AtomRadiiTable& Molecule::getAtomRadiiTable() const
{
    static AtomRadiiTable emptytable;
    const AtomRadiiTable& rv = (_atom_radii_table.get()) ?
        (*_atom_radii_table) : emptytable;
    return rv;
}


void Molecule::setSamePairRadius(double samepairradius)
{
    this->_samepairradius = samepairradius;
}


const double& Molecule::getSamePairRadius() const
{
    return this->_samepairradius;
}


int Molecule::getMaxAtomCount() const
{
    return atoms_storage.size();
}


void Molecule::Shift(const R3::Vector& drc)
{
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
        Atom_t* pa = seq.ptr();
        pa->r += drc;
    }
}


void Molecule::Center()
{
    R3::Vector avg_rc(0.0, 0.0, 0.0);
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
        Atom_t* pa = seq.ptr();
        avg_rc += pa->r;
    }
    avg_rc /= countAtoms();
    Shift(-avg_rc);
}


AtomPtr Molecule::getNearestAtom(const R3::Vector& rc) const
{
    R3::Vector rij(0.0, 0.0, 0.0);
    double mindistance = DOUBLE_MAX;
    int index = -1;
    for (AtomSequenceIndex seq(this); !seq.finished(); seq.next())
    {
        rij = seq.ptr()->r - rc;
        if (R3::norm(rij) >= mindistance)   continue;
        mindistance = R3::norm(rij);
        index = seq.idx();
    }
    AtomPtr rv;
    if (index >= 0)  rv.reset(new Atom_t(this->getAtom(index)));
    return rv;
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
    atoms.erase(atoms.begin() + aidx);
    atoms_bucket.push_back(pa);
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
    // move all atoms to the bucket
    atoms_bucket.insert(atoms_bucket.end(), atoms.begin(), atoms.end());
    assert(atoms_bucket.size() == atoms_storage.size());
    atoms.clear();
    free_pmx_slots.clear();
    ResetBadness();
    ResetOverlap();
    CheckIntegrity();
}


void Molecule::AddAt(const string& smbl, double rx0, double ry0, double rz0)
{
    vector<Atom_t*>::iterator ai;
    for (ai = atoms_bucket.begin(); ai != atoms_bucket.end(); ++ai)
    {
        if ((*ai)->element == smbl)  break;
    }
    if (ai == atoms_bucket.end())
    {
        ostringstream emsg;
        emsg << "Cannot add '" << smbl << "', element not available.";
        throw invalid_argument(emsg.str());
    }
    this->AddInternalAt(*ai, rx0, ry0, rz0);
}


void Molecule::AddAt(const string& smbl, const R3::Vector& rc)
{
    this->AddAt(smbl, rc[0], rc[1], rc[2]);
}


void Molecule::Add(const Atom_t& a)
{
    this->AddAt(a.element, a.r);
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


namespace {
bool pAtom_is_fixed(const Atom_t* pa)
{
    return pa->fixed;
}
}   // namespace


int Molecule::NFixed() const
{
    return count_if(atoms.begin(), atoms.end(), pAtom_is_fixed);
}


void Molecule::filter_good_atoms(AtomArray& vta,
        double evolve_range, double hi_abad)
{
    if (countAtoms() == getMaxAtomCount())
    {
        ostringstream emsg;
        emsg << "E: Molecule too large in filter_good_atoms()";
        throw InvalidMolecule(emsg.str());
    }
    typedef AtomArray::iterator VAit;
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


void Molecule::filter_bucket_atoms(AtomArray& vta)
{
    set<Atom_t*> bucket(atoms_bucket.begin(), atoms_bucket.end());
    assert(bucket.size() == atoms_bucket.size());
    AtomArray::iterator asrc = vta.begin();
    AtomArray::iterator adst = vta.begin();
    for (; asrc != vta.end(); ++asrc)
    {
        if (bucket.count(asrc->mstorage_ptr))
        {
            *adst = *asrc;
            ++adst;
        }
    }
    vta.erase(adst, vta.end());
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

// GSL helper functions for RelaxExternalAtom
namespace {

template <class V> void copyGSLvector(const gsl_vector* src, V& dest)
{
    for (size_t i = 0; i < src->size; i++)
    {
        dest[i] = gsl_vector_get(src, i);
    }
}


template <class V> void copyGSLvector(const V& src, gsl_vector* dest)
{
    for (size_t i = 0; i < dest->size; i++)
    {
        gsl_vector_set(dest, i, src[i]);
    }
}


struct rxa_fg
{
    double f;
    R3::Vector g;
};

struct rxa_compare_R3Vector
{
    bool operator()(const R3::Vector& v0, const R3::Vector& v1) const
    {
        return lexicographical_compare(
                v0.data(), v0.data() + R3::Ndim,
                v1.data(), v1.data() + R3::Ndim);
    }
};

const Molecule* rxa_molecule = NULL;
AtomCost* rxa_atomcost = NULL;
AtomCost* rxa_atomoverlap = NULL;

class rxa_CacheType : public map<R3::Vector,rxa_fg,rxa_compare_R3Vector>
{
    public:

        const mapped_type* lookup(const R3::Vector& x)
        {
            const mapped_type* rv = NULL;
            if (this->empty())  return rv;
            const_iterator kv = this->lower_bound(x);
            if (kv != this->end() &&
                R3::norm(x - kv->first) < eps_distance)
            {
                rv = &(kv->second);
            }
            else if (kv != this->begin() &&
                     R3::norm(x - (--kv)->first) < eps_distance)
            {
                rv = &(kv->second);
            }
            return rv;
        }
};

rxa_CacheType rxa_cache;

void rxa_fdf(const gsl_vector* x, void* params, double* f, gsl_vector* g);

double rxa_f(const gsl_vector* x, void* params)
{
    double rv;
    gsl_vector* g = NULL;
    rxa_fdf(x, params, &rv, g);
    return rv;
}


void rxa_df(const gsl_vector* x, void* params, gsl_vector* g)
{
    double f;
    rxa_fdf(x, params, &f, g);
}


void rxa_fdf(const gsl_vector* x, void* params, double* f, gsl_vector* g)
{
    Atom_t* pa = static_cast<Atom_t*>(params);
    Atom_t& ta = *pa;
    copyGSLvector(x, ta.r);
    // try to get the result from the cache
    const rxa_fg* pfg = rxa_cache.lookup(ta.r);
    if (!pfg)
    {
        // calculate and store in the cache
        ta.ResetBadness(rxa_atomcost->eval(ta, AtomCost::GRADIENT));
        ta.ResetOverlap(rxa_atomoverlap->eval(ta, AtomCost::GRADIENT));
        rxa_fg fg;
        fg.f = ta.costShare(rxa_molecule->pairsPerAtomInc());
        fg.g = 0.0, 0.0, 0.0;
        fg.g += rxa_atomcost->gradient();
        fg.g += rxa_molecule->pairsPerAtomInc() * rxa_atomoverlap->gradient();
        rxa_cache[ta.r] = fg;
        pfg = rxa_cache.lookup(ta.r);
        assert(R3::norm(fg.g - pfg->g) == 0.0);
        assert(fg.f == pfg->f);
    }
    // here pfg should be non NULL
    assert(pfg);
    *f = pfg->f;
    if (g)  copyGSLvector(pfg->g, g);
}


void rxa_useMolecule(const Molecule* mol)
{
    rxa_molecule = mol;
    rxa_atomcost = mol ? mol->getAtomCostCalculator() : NULL;
    rxa_atomoverlap = mol ? mol->getAtomOverlapCalculator() : NULL;
    rxa_cache.clear();
}

}   // namespace


void Molecule::RelaxAtom(const int cidx)
{
    if (cidx < 0 || cidx >= countAtoms())
    {
        ostringstream emsg;
        emsg << "Molecule::RelaxAtom(const int) invalid index " << cidx << '.';
        throw range_error(emsg.str());
    }
    // RelaxAtom should never get called on a fixed atom
    Atom_t* pa = atoms[cidx];
    assert(!pa->fixed);
    Pop(cidx);
    RelaxExternalAtom(pa);
    AddInternal(pa);
}


void Molecule::RelaxExternalAtom(Atom_t* pa)
{
    // Configuration of the GSL multidimensional minimizer
    // For details, see info pages
    //     '(gsl-ref)Multidimensional Minimization'
    //     '(gsl-ref)Initializing the Multidimensional Minimizer'
    const int maximum_iterations = 500;
    const gsl_multimin_fdfminimizer_type* minimizer_type;
    minimizer_type = gsl_multimin_fdfminimizer_vector_bfgs2;
    const double minimizer_step = 1.0e-4;
    const double minimizer_tol = 0.1;
    // do relaxation on a copy of pa
    Atom_t rta(*pa);
    // loop while badness is improved
    gsl_multimin_function_fdf fdfmin;
    rxa_useMolecule(this);
    fdfmin.f = &rxa_f;
    fdfmin.df = &rxa_df;
    fdfmin.fdf = &rxa_fdf;
    // FIXME handling of ndim
    fdfmin.n = 3;
    fdfmin.params = &rta;
    // copy rta coordinates to vector x
    gsl_vector* x = gsl_vector_alloc(3);
    copyGSLvector(rta.r, x);
    // determine initial_cost before relaxation
    double initial_cost = rxa_f(x, fdfmin.params);
    // allocate minimizer
    gsl_multimin_fdfminimizer* minimizer;
    minimizer = gsl_multimin_fdfminimizer_alloc(minimizer_type, fdfmin.n);
    gsl_multimin_fdfminimizer_set(minimizer,
            &fdfmin, x, minimizer_step, minimizer_tol);
    // iterate unless initial_cost is close to zero
    int iter, status;
    bool skipiter = isNearZeroRoundOff(initial_cost);
    for (iter = 0; iter < maximum_iterations && !skipiter; ++iter)
    {
        status = gsl_multimin_fdfminimizer_iterate(minimizer);
        if (status != GSL_SUCCESS)  break;
        status = gsl_multifit_test_delta(minimizer->dx,
                minimizer->x, eps_distance, eps_distance);
        if (status == GSL_SUCCESS)
        {
            // one last iteration to improve the precision
            gsl_multimin_fdfminimizer_iterate(minimizer);
            break;
        }
        if (status != GSL_CONTINUE) break;
    }
    // copy new coordinates to rta
    // update relaxed atom if cost improved
    double relaxed_cost = minimizer->f;
    if (eps_lt(relaxed_cost, initial_cost))
    {
        copyGSLvector(minimizer->x, pa->r);
    }
    // release minimizer objects
    gsl_multimin_fdfminimizer_free(minimizer);
    gsl_vector_free(x);
    rxa_useMolecule(NULL);
}


void Molecule::AddInternalAt(Atom_t* pa, double rx0, double ry0, double rz0)
{
    pa->r = rx0, ry0, rz0;
    AddInternal(pa);
}


void Molecule::AddInternalAt(Atom_t* pa, const R3::Vector& rc)
{
    pa->r = rc;
    AddInternal(pa);
}


void Molecule::AddInternal(Atom_t* pa)
{
    vector<Atom_t*>::iterator ai;
    ai = find(atoms_bucket.begin(), atoms_bucket.end(), pa);
    assert(ai != atoms_bucket.end());
    // reset cost related attributes
    pa->ResetBadness();
    pa->ResetOverlap();
    pa->pmxidx = getPairMatrixIndex();
    // create new pairs while summing up the costs
    addNewAtomPairs(pa);
    atoms_bucket.erase(ai);
    atoms.push_back(pa);
    if (full())     reassignPairs();
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
    if (isNearZeroRoundOff(this->Badness()))  this->ResetBadness();
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
    // add overlap contributions
    this->applyOverlapContributions(pa, ADD);
}


void Molecule::removeAtomPairs(Atom_t* pa)
{
    // remove associated pair costs
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
        if (pa == seq.ptr())    continue;
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
            if (pa == seq.ptr())        continue;
            // return any used distances
            int idx0 = pa->pmxidx;
            int idx1 = seq.ptr()->pmxidx;
            double& udst = pmx_used_distances(idx0, idx1);
            if (udst > 0.0)     this->_distance_table->return_back(udst);
            udst = 0.0;
        }
    }
    if (isNearZeroRoundOff(this->Badness()))  this->ResetBadness();
    // remove overlap contributions
    this->applyOverlapContributions(pa, REMOVE);
}


Atom_t* Molecule::pickAtomFromBucket() const
{
    assert(!atoms_bucket.empty());
    int idx = randomInt(atoms_bucket.size());
    return atoms_bucket[idx];
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
        Atom_t ad1front(pickAtomFromBucket(), nr);
        ad1front.ttp = LINEAR;
        vta.push_back(ad1front);
        ++push_count;
        // check opposite direction when it makes sense
        // this accounts for extra trial
        if (lattice_direction)
        {
            ++nt;
            nr = anch.B0 - direction*radius;
            Atom_t ad1back(pickAtomFromBucket(), nr);
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
            int ijk = minIndex(ald)[0];
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
                Atom_t ad2(pickAtomFromBucket(), P);
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
                Atom_t ad3(pickAtomFromBucket(), P4);
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
                Atom_t ad3(pickAtomFromBucket(), P4);
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
            Atom_t ad3top(pickAtomFromBucket(), P4);
            ad3top.ttp = SPATIAL;
            vta.push_back(ad3top);
            ++push_count;
            // and bottom one, which makes an extra trial
            ++nt;
            P4 = yP4*uvj - zP4*uvk + vT;
            Atom_t ad3bottom(pickAtomFromBucket(), P4);
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
    assert(!atoms_bucket.empty());
    // containter for test atoms
    AtomArray vta;
    // evolution is trivial for empty or 1-atom molecule
    switch (countAtoms())
    {
        case 0:
            AddInternalAt(pickAtomFromBucket(), 0.0, 0.0, 0.0);
            acc[LINEAR] = 1;
            tot[LINEAR] = 1;
            return acc_tot;
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
    typedef AtomArray::iterator VAit;
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
        filter_bucket_atoms(vta);
        filter_good_atoms(vta, evolve_range, hi_abad);
        // finished when no test atoms left
        if (vta.empty())   break;
        // calculate fitness of test atoms
        valarray<double> vtafit(vta.size());
        // calculate fitness as reciprocal value of badness
        // fill the vtafit array with badness
        double* pfit = &vtafit[0];
        double ppa = this->pairsPerAtomInc();
        for (VAit ai = vta.begin(); ai != vta.end(); ++ai, ++pfit)
        {
            *pfit = ai->costShare(ppa);
        }
        // then get the reciprocal value
        double* ftnfirst = &(vtafit[0]);
        double* ftnlast = &(vtafit[vtafit.size()]);
        transform(ftnfirst, ftnlast, ftnfirst, convertCostToFitness);
        // vtafit is ready here
        int idx = randomWeighedInt(vtafit.size(), &vtafit[0]);
        AddInternalAt(vta[idx].mstorage_ptr, vta[idx].r);
        acc[vta[idx].ttp]++;
        hi_abad = vta[idx].Badness() + evolve_range;
        vta.erase(vta.begin()+idx);
        if (true)
        {
            int worst_overlap_idx = max_element(atoms.begin(), atoms.end(),
                    comp_pAtom_FreeOverlap) - atoms.begin();
            this->MinimizeSiteOverlap(worst_overlap_idx);
        }
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


void Molecule::Degenerate(int Npop, DegenerateFlags flags)
{
    Npop = min(countAtoms(), Npop);
    if (Npop == 0)  return;
    // build array of atom badnesses
    double freebad[countAtoms()];
    int freeidx[countAtoms()];
    int Nfree = 0;
    double ppa = this->pairsPerAtom();
    for (int i = 0; i != countAtoms(); ++i)
    {
        Atom_t* pai = atoms[i];
        if ( pai->fixed )  continue;
        freebad[Nfree] = pai->costShare(ppa);
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
    if (this->countAtoms() && flags != FAST)
    {
        int worst_overlap_idx = max_element(atoms.begin(), atoms.end(),
                comp_pAtom_FreeOverlap) - atoms.begin();
        this->MinimizeSiteOverlap(worst_overlap_idx);
    }
    if (demoterelax && countAtoms() > 1 && flags != FAST)
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


double Molecule::getContactRadius(const Atom_t& a0, const Atom_t& a1) const
{
    double rv = (a0.radius < 0 || a1.radius < 0) ? 0.0 :
        (this->getSamePairRadius() >= 0 && a0.element == a1.element) ?
        (2 * this->getSamePairRadius()) : (a0.radius + a1.radius);
    return rv;
}


void Molecule::FlipSites(int idx0, int idx1)
{
    this->checkAtomIndex(idx0);
    this->checkAtomIndex(idx1);
    Atom_t* pa0 = this->atoms[idx0];
    Atom_t* pa1 = this->atoms[idx1];
    // short circuit when overlap cost does not change
    if (pa0->radius == pa1->radius)  return;
    double radius0 = pa0->radius;
    double radius1 = pa1->radius;
    this->applyOverlapContributions(pa0, REMOVE);
    pa0->radius = -DOUBLE_MAX;
    this->applyOverlapContributions(pa1, REMOVE);
    swap(pa0->element, pa1->element);
    pa1->radius = radius0;
    this->applyOverlapContributions(pa1, ADD);
    pa0->radius = radius1;
    this->applyOverlapContributions(pa0, ADD);
}


void Molecule::DownhillOverlapMinimization()
{
    if (this->getMaxAtomRadius() <= 0.0)  return;
    auto_ptr<Molecule> m1(this->copy());
    while (true)
    {
        struct { double overlap; int i, j; } best = {this->Overlap(), 0, 0};
        for (int i = 0; i < this->countAtoms(); ++i)
        {
            for (int j = i + 1; j < this->countAtoms(); ++j)
            {
                m1->FlipSites(i, j);
                if (m1->Overlap() < best.overlap)
                {
                    best.overlap = m1->Overlap();
                    best.i = i;
                    best.j = j;
                }
                m1->FlipSites(i, j);
                assert(eps_eq(this->Overlap(), m1->Overlap()));
            }
        }
        // get out if overlap did not improve
        if (!eps_lt(best.overlap, this->Overlap()))    break;
        assert(best.i != best.j);
        // perform this flip
        this->FlipSites(best.i, best.j);
        m1->FlipSites(best.i, best.j);
#ifndef NDEBUG
        this->recalculateOverlap();
        assert(eps_eq(this->Overlap(), m1->Overlap()));
#endif
    }
}


void Molecule::MinimizeSiteOverlap(int idx0)
{
    if (this->getMaxAtomRadius() <= 0.0)  return;
    // no need to check atom index, it will be checked in FlipSites
    double initial_overlap = this->Overlap();
    struct { double overlap; int idx1; } best = {initial_overlap, idx0};
    for (int idx1 = 0; idx1 < this->countAtoms(); ++idx1)
    {
        this->FlipSites(idx0, idx1);
        if (this->Overlap() < best.overlap)
        {
            best.overlap = this->Overlap();
            best.idx1 = idx1;
        }
        this->FlipSites(idx0, idx1);
    }
    assert(eps_eq(initial_overlap, this->Overlap()));
    // perform the best flip
    this->FlipSites(idx0, best.idx1);
    assert(!eps_gt(this->Overlap(), initial_overlap));
}


string Molecule::PickElementFromBucket() const
{
    const Atom_t* pa = this->pickAtomFromBucket();
    return pa->element;
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
            if (udst > 0.0)     this->_distance_table->push_back(udst);
            udst = 0.0;
        }
    }
    sort(this->_distance_table->begin(), this->_distance_table->end());
}

// Molecule IO Functions -----------------------------------------------------

void Molecule::ReadFile(const string& filename)
{
    namespace python = boost::python;
    try {
        initializePython();
        python::object stru = this->newDiffPyStructure();
        python::call_method<void>(stru.ptr(), "read", filename);
        this->setFromDiffPyStructure(stru);
    }
    catch (python::error_already_set) {
        if (PyErr_Occurred())   PyErr_Print();
        const char* emsg = "Cannot read structure.";
        throw IOError(emsg);
    }
    this->CheckIntegrity();
}


void Molecule::WriteFile(const string& filename, string title)
{
    // invalid structure format can throw exception,
    // check if write operator works first
    ostringstream output;
    this->WriteStream(output, title);
    // test if filename is writeable
    ofstream fid(filename.c_str(), ios_base::out|ios_base::ate);
    if (!fid)
    {
        ostringstream emsg;
        emsg << "WriteFile(): unable to write to '" << filename << "'";
        throw IOError(emsg.str());
    }
    fid.close();
    // write via temporary file
    string writefile = filename + "XXXXXX";
    mktempofstream(fid, writefile);
    fid << output.str();
    fid.close();
    rename(writefile.c_str(), filename.c_str());
}


void Molecule::setOutputFormat(const std::string& format)
{
    Molecule::output_format = format;
}


void Molecule::WriteStream(ostream& fid, string title) const
{
    namespace python = boost::python;
    try {
        python::object stru = this->convertToDiffPyStructure();
        stru.attr("title") = title;
        string s;
        s = python::call_method<string>(stru.ptr(),
                "writeStr", Molecule::output_format);
        fid << s;
    }
    catch (python::error_already_set) {
        if (PyErr_Occurred())   PyErr_Print();
        const char* emsg = "Cannot output structure.";
        throw IOError(emsg);
    }
}


boost::python::object Molecule::newDiffPyStructure() const
{
    namespace python = boost::python;
    python::object mstru = python::import("diffpy.structure");
    python::object stru = mstru.attr("Structure")();
    return stru;
}


boost::python::object Molecule::convertToDiffPyStructure() const
{
    namespace python = boost::python;
    initializePython();
    python::object stru = this->newDiffPyStructure();
    for (AtomSequence seq(this); !seq.finished(); seq.next())
    {
        const Atom_t& ai = seq.ref();
        stru.attr("addNewAtom")(ai.element);
        python::object alast = stru.attr("getLastAtom")();
        python::object xyz_cartn;
        xyz_cartn = python::make_tuple(ai.r[0], ai.r[1], ai.r[2]);
        alast.attr("xyz_cartn") = xyz_cartn;
        alast.attr("cost") = ai.Badness();
    }
    // xcfg format can save atom cost as an auxiliary property
    if (Molecule::output_format == "xcfg")
    {
        python::list auxiliaries;
        auxiliaries.append("cost");
        python::dict xcfg;
        xcfg["auxiliaries"] = auxiliaries;
        stru.attr("xcfg") = xcfg;
    }
    return stru;
}


void Molecule::setFromDiffPyStructure(boost::python::object stru)
{
    namespace python = boost::python;
    int num_atoms = python::len(stru);
    this->Clear();
    atoms_storage.clear();
    atoms_bucket.clear();
    for (int i = 0; i != num_atoms; ++i)
    {
        python::object ai;
        python::object xyz_cartn;
        double x, y, z;
        ai = stru[i];
        xyz_cartn = ai.attr("xyz_cartn");
        x = python::extract<double>(xyz_cartn[0]);
        y = python::extract<double>(xyz_cartn[1]);
        z = python::extract<double>(xyz_cartn[2]);
        string smbl = python::extract<string>(ai.attr("element"));
        atoms_storage.push_back(Atom_t(smbl, x, y, z));
        Atom_t* pa = &(atoms_storage.back());
        atoms_bucket.push_back(pa);
    }
    this->fetchAtomRadii();
    list<Atom_t>::iterator ai = atoms_storage.begin();
    for (; ai != atoms_storage.end(); ++ai)
    {
        Atom_t* pa = &(*ai);
        this->AddInternal(pa);
    }
}


void Molecule::recalculateOverlap() const
{
    this->ResetOverlap();
    AtomCost* atomoverlap = this->getAtomOverlapCalculator();
    BOOST_FOREACH (Atom_t* pa, this->atoms)
    {
        double aoself = atomoverlap->eval(pa, AtomCost::SELFCOST);
        double aohalf = atomoverlap->eval(pa) / 2.0;
        pa->ResetOverlap(aoself + aohalf);
        this->IncOverlap(aoself + aohalf);
    }
    BOOST_FOREACH (Atom_t* pa, this->atoms_bucket)
    {
        pa->ResetOverlap(0.0);
    }
}


void Molecule::applyOverlapContributions(Atom_t* pa, AddRemove sign)
{
    // add overlap contributions
    AtomCost* atomoverlap = getAtomOverlapCalculator();
    // self contribution
    double aoself = atomoverlap->eval(pa, AtomCost::SELFCOST);
    pa->IncOverlap(sign * aoself);
    this->IncOverlap(sign * aoself);
    // cross contributions
    atomoverlap->eval(pa);
    for (AtomSequenceIndex seq(this); !seq.finished(); seq.next())
    {
        double aohalf = sign * atomoverlap->partialCosts()[seq.idx()] / 2.0;
        seq.ptr()->IncOverlap(aohalf);
        pa->IncOverlap(aohalf);
    }
    this->IncOverlap(sign * atomoverlap->totalCost());
    if (sign == REMOVE && isNearZeroRoundOff(this->Overlap()))
    {
        this->ResetOverlap();
    }
}


void Molecule::fetchAtomRadii()
{
    const AtomRadiiTable& radiitable = this->getAtomRadiiTable();
    BOOST_FOREACH (Atom_t& a, atoms_storage)
    {
        a.radius = (radiitable.empty() || a.element.empty()) ?
            0.0 : radiitable.lookup(a.element);
    }
}


void Molecule::checkAtomIndex(int idx)
{
    if (idx < 0 || idx >= this->countAtoms())
    {
        ostringstream emsg;
        emsg << "Invalid atom index " << idx << ".";
        throw range_error(emsg.str());
    }
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


void Molecule::CheckIntegrity() const
{
#ifndef NDEBUG
    set<const Atom_t*> set_atoms(atoms.begin(), atoms.end());
    set<const Atom_t*> set_bucket(atoms_bucket.begin(), atoms_bucket.end());
    assert(set_atoms.size() + set_bucket.size() == atoms_storage.size());
    set<const Atom_t*> set_storage;
    BOOST_FOREACH (const Atom_t& a, atoms_storage)
    {
        set_storage.insert(&a);
    }
    assert(set_storage.size() == atoms_storage.size());
    BOOST_FOREACH (const Atom_t* pa, set_atoms)
    {
        assert(set_storage.count(pa));
    }
    BOOST_FOREACH (const Atom_t* pa, set_bucket)
    {
        assert(set_storage.count(pa));
    }
#endif  // NDEBUG
}


R3::Vector Molecule::rxaCheckGradient(const Atom_t* pa) const
{
    double c;
    R3::Vector rv;
    this->rxaCheckEval(pa, &c, &rv);
    return rv;
}


double Molecule::rxaCheckCost(const Atom_t* pa) const
{
    double c;
    R3::Vector rv;
    this->rxaCheckEval(pa, &c, &rv);
    return c;
}


void Molecule::rxaCheckEval(const Atom_t* pa,
        double* pcost, R3::Vector* pg) const
{
    rxa_useMolecule(this);
    Atom_t rta = *pa;
    void* params = &rta;
    // copy rta coordinates to vector x
    gsl_vector* x = gsl_vector_alloc(3);
    gsl_vector_set(x, 0, rta.r[0]);
    gsl_vector_set(x, 1, rta.r[1]);
    gsl_vector_set(x, 2, rta.r[2]);
    gsl_vector* g = gsl_vector_alloc(3);
    // calculate both cost and gradient
    rxa_fdf(x, params, pcost, g);
    copyGSLvector(g, *pg);
    gsl_vector_free(g);
    gsl_vector_free(x);
    rxa_useMolecule(NULL);
}

// Non-member Operators ------------------------------------------------------

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


ostream& operator<<(ostream& fid, const Molecule& M)
{
    M.WriteStream(fid);
    return fid;
}


// End of file
