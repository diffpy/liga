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

#ifndef ATOMOVERLAP_HPP_INCLUDED
#define ATOMOVERLAP_HPP_INCLUDED

#include "AtomCost.hpp"

class AtomOverlap : public AtomCost
{
    public:

        // constructor
        AtomOverlap(const Molecule* m);

        // public methods - overloaded
        virtual void resetFor(const Molecule* clust);
        double eval(const Atom_t& a, int flags=NONE);
        virtual double eval(const Atom_t* pa, int flags=NONE);

};  // class AtomOverlap

#endif  // ATOMOVERLAP_HPP_INCLUDED
