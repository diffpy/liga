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
* class AtomOverlap
*
* Comments: OverlapCost calculation for a Molecule
*
***********************************************************************/

#include <vector>

#include "AtomOverlap.hpp"
#include "AtomSequence.hpp"
#include "Molecule.hpp"
#include "LigaUtils.hpp"
#include "Counter.hpp"

using namespace std;

// Constructor ---------------------------------------------------------------

AtomOverlap::AtomOverlap(const Molecule* m) : AtomCost(m)
{ }

// Public Methods ------------------------------------------------------------

void AtomOverlap::resetFor(const Molecule* m)
{
    arg_cluster = m;
    noCutoff();
    use_distances = false;
}


double AtomOverlap::eval(const Atom_t& a, int flags)
{
    return this->eval(&a, flags);
}


double AtomOverlap::eval(const Atom_t* pa, int flags)
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
    fill(partial_costs.begin(), partial_costs.end(), 0.0);
    if (this->_selfcost_flag)  return this->totalCost();
    for (AtomSequenceIndex seq(arg_cluster); !seq.finished(); seq.next())
    {
        // do not calculate own overlap
        if (seq.ptr() == arg_atom)  continue;
        // calculation
        double d = R3::distance(arg_atom->r, seq.ptr()->r);
        // force zero overlap when disabled by negative atom radius
        double r0r1 = arg_cluster->getContactRadius(*arg_atom, seq.ref());
        double dd = (d < r0r1) ? (r0r1 - d) : 0.0;
        double pcost = this->penaltyScaled(dd, 1.0);
        partial_costs[seq.idx()] = pcost;
        total_cost += pcost;
        if (this->_gradient_flag && d > NS_LIGA::eps_distance)
        {
            static R3::Vector g_dd_xyz;
            g_dd_xyz = (-1.0/d) * (arg_atom->r - seq.ptr()->r);
            double g_pcost_dd = penalty_gradient(dd, 1.0) * this->getScale();
            this->_gradient += g_pcost_dd * g_dd_xyz;
        }
    }
    this->_gradient_cached = this->_gradient_flag;
    return total_cost;
}

// End of file
