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
* $Id$
*
***********************************************************************/

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
    pair<double,double> rext = arg_cluster->getRExtent();
    _sph.reset(new PointsInSphere(rext.first, rext.second,
                arg_cluster->getLattice()) );
}


double AtomOverlapCrystal::pairDistanceDifference(const double& d) const
{
    double r0r1 = this->arg_atom->radius + this->crst_atom->radius;
    double dd = (d < r0r1) ? (r0r1 - d) : 0.0;
    return dd;
}


// End of file