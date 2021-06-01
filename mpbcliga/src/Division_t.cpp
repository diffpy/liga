/***********************************************************************
* Short Title: one division of the liga system
*
* Comments:
*
* <license text>
***********************************************************************/

#include <cassert>

#include "Division_t.hpp"
#include "Molecule.hpp"
#include "LigaUtils.hpp"

using namespace std;
using namespace NS_LIGA;

////////////////////////////////////////////////////////////////////////
// local helper routines
////////////////////////////////////////////////////////////////////////

namespace {

typedef Division_t::PMOL PMOL;

inline bool comp_PMOL_Badness( const PMOL lhs, const PMOL rhs)
{
    return lhs->Badness() < rhs->Badness();
}

}   // namespace


////////////////////////////////////////////////////////////////////////
// class Division_t
////////////////////////////////////////////////////////////////////////

// constructors and destructor

Division_t::Division_t(size_t fullsize, size_t level) :
    vector<PMOL>(),
    _fullsize(fullsize), _level(level), _trials(0.0)
{
    fill(acc_triang, acc_triang + NTGTYPES, 0);
    fill(tot_triang, tot_triang + NTGTYPES, 0);
    fill(est_triang, est_triang + NTGTYPES, 0);
}

Division_t::Division_t(const Division_t& src) :
    vector<PMOL>(src),
    _fullsize(src._fullsize), _level(src._level), _trials(0.0)
{
    fill(acc_triang, acc_triang + NTGTYPES, 0);
    fill(tot_triang, tot_triang + NTGTYPES, 0);
    fill(est_triang, est_triang + NTGTYPES, 0);
}

Division_t::~Division_t()
{
    for (iterator ii = begin(); ii != end(); ++ii)  delete *ii;
}

// operators

Division_t& Division_t::operator= (const Division_t& src)
{
    static_cast< vector<PMOL> >(*this) = src;
    _fullsize = src._fullsize;
    _level = src._level;
    copy(src.acc_triang, src.acc_triang + NTGTYPES, acc_triang);
    copy(src.tot_triang, src.tot_triang + NTGTYPES, tot_triang);
    copy(src.est_triang, src.est_triang + NTGTYPES, est_triang);
    return *this;
}

// public methods

int Division_t::find_winner()
{
    // evaluate molecule fitness
    valarray<double> vmfit(size());
    double *pd = &vmfit[0];
    for (iterator mi = begin(); mi != end(); ++mi, ++pd)
        *pd = (*mi)->cost();
    // then get the reciprocal value
    vmfit = costToFitness(vmfit);
    double *mfit = &vmfit[0];
    int idx = randomWeighedInt(size(), mfit);
    return idx;
}

int Division_t::find_looser()
{
    // evaluate molecule fitness
    valarray<double> vmbad(size());
    double *pd = &vmbad[0];
    for (iterator mi = begin(); mi != end(); ++mi, ++pd)
        *pd = (*mi)->cost();
    double *mbad = &vmbad[0];
    int idx = randomWeighedInt(size(), mbad);
    return idx;
}

int Division_t::find_best()
{
    int idx = min_element(begin(), end(), comp_PMOL_Badness) - begin();
    return idx;
}

Division_t::PMOL& Division_t::best()
{
    return at(find_best());
}

double Division_t::averageCost() const
{
    double total = 0.0;
    for (const_iterator ii = begin(); ii != end(); ++ii)
    {
        total += (*ii)->cost();
    }
    return size() ? total / size() : 0.0;
}

const int* Division_t::estimateTriangulations()
{
    const double pdef[NTGTYPES] = { 2.0/18, 4.0/18, 12.0/18 };
    valarray<double> pbtg(pdef, NTGTYPES);
    for (size_t i = 0; i != NTGTYPES; ++i)
    {
        if (!tot_triang[i])     continue;
        // probability of successful triangulation is given
        // by beta distribution
        double a = acc_triang[i] + 1;
        double b = tot_triang[i] - acc_triang[i] + 1;
        pbtg[i] = randomBeta(a, b);
    }
    // nasty hack to fix triangulation probabilities for Crystal object
    // Triangulation statistics should be moved to Molecule and Crystal.
    // size_t nd = min(ndim, level());
    size_t nd = ndim;
    if (!this->empty() && this->front()->type() == MOLECULE)
    {
        nd = min(ndim, this->level());
    }
    switch (nd)
    {
        case 0:
            pbtg[LINEAR] = 0.0;
        case 1:
            pbtg[PLANAR] = 0.0;
        case 2:
            pbtg[SPATIAL] = 0.0;
        default:
            double ptot = pbtg.sum();
            if (ptot)   pbtg /= ptot;
    }
    for (size_t i = 0; i != NTGTYPES; ++i)
    {
        est_triang[i] = int(ceil(pbtg[i]*_trials));
    }
    return est_triang;
}

void Division_t::noteTriangulations(const pair<int*,int*>& acc_tot)
{
    const int* acc = acc_tot.first;
    const int* tot = acc_tot.second;
    for (int i = 0; i < NTGTYPES; ++i)
    {
        acc_triang[i] += acc[i];
        tot_triang[i] += tot[i];
        assert(acc_triang[i] <= tot_triang[i]);
    }
}

// class data

size_t Division_t::ndim = 3;

// End of file
