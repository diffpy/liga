/***********************************************************************
*
* Liga Algorithm    for structure determination from pair distances
*                   Pavol Juhas
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* class AtomOverlapCrystal
*
* Comments: OverlapCost calculation for a Molecule
*
***********************************************************************/

#include <cassert>
#include <vector>

#include "AtomOverlapCrystal.hpp"
#include "AtomSequence.hpp"
#include "Crystal.hpp"
#include "Lattice.hpp"
#include "LigaUtils.hpp"
#include "Counter.hpp"

using namespace std;

// Constructor ---------------------------------------------------------------

AtomOverlapCrystal::AtomOverlapCrystal(const Crystal* crst) :
    AtomCostCrystal(crst)
{ }

// Public Methods ------------------------------------------------------------

void AtomOverlapCrystal::resetFor(const Molecule* clust)
{
    this->AtomCost::resetFor(clust);
    this->arg_cluster = static_cast<const Crystal*>(clust);
    assert(this->arg_cluster == this->AtomCost::arg_cluster);
    assert(!this->use_distances);
    this->_rmax = 2 * arg_cluster->getMaxAtomRadius();
    pair<double,double> rext = arg_cluster->getRExtent(0.0, this->_rmax);
    _sph.reset(new PointsInSphere(rext.first, rext.second,
                arg_cluster->getLattice()) );
}


const pair<double,double>&
AtomOverlapCrystal::pairDistanceDifference(const double& d) const
{
    static pair<double,double> rv(0.0, 1.0);
    assert(1.0 == rv.second);
    // force zero overlap when disabled by negative atom radius
    double r0r1 = arg_cluster->getContactRadius(*arg_atom, *crst_atom);
    rv.first = (d < r0r1) ? (r0r1 - d) : 0.0;
    return rv;
}


// End of file
