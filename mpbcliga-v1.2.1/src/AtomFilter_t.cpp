/***********************************************************************
* Short Title: definitions of AtomFilter_t and derived classes
*
* Comments:
*
* <license text>
***********************************************************************/

#include <list>
#include "AtomFilter_t.hpp"
#include "LigaUtils.hpp"
#include "AtomSequence.hpp"
#include "Molecule.hpp"

using namespace NS_LIGA;
using namespace std;


////////////////////////////////////////////////////////////////////////
// class BondAngleFilter_t - derived from AtomFilter_t
////////////////////////////////////////////////////////////////////////

// constructor

BondAngleFilter_t::BondAngleFilter_t(double _max_blen)
{
    max_blen = _max_blen;
    setBondAngleRange(0.0, DOUBLE_MAX);
}

// public methods

bool BondAngleFilter_t::Check(Atom_t* pta, Molecule* pm)
{
    // first neighbors are closer than max_blen
    list<Atom_t*> first_neighbors;
    list<double>  first_distances;
    // second neighbors are closer than 2*max_blen
    list<Atom_t*> second_neighbors;
    list<double>  second_distances;
    // find first and second neighbors and corresponding distances
    typedef vector<Atom_t*>::iterator VPAit;
    for (AtomSequence seq(pm); !seq.finished(); seq.next())
    {
        double blen = R3::distance(pta->r, seq.ptr()->r);
        if (0.0 < blen && blen < 2*max_blen)
        {
            second_neighbors.push_back(seq.ptr());
            second_distances.push_back(blen);
            if (blen < max_blen)
            {
                first_neighbors.push_back(seq.ptr());
                first_distances.push_back(blen);
            }
        }
    }
    // check all new bond angles formed by test atom
    typedef list<Atom_t*>::iterator LPAit;
    typedef list<double>::iterator LDit;
    LPAit pfn1i, pfn2i;
    LDit dfn1i, dfn2i;
    for (   pfn1i = first_neighbors.begin(), dfn1i = first_distances.begin();
            pfn1i != first_neighbors.end();  ++pfn1i, ++dfn1i )
    {
        // check test atom angles, ta_angle
        pfn2i = pfn1i; ++pfn2i;
        dfn2i = dfn1i; ++dfn2i;
        for (; pfn2i != first_neighbors.end(); ++pfn2i, ++dfn2i)
        {
            double dneib = R3::distance((*pfn1i)->r, (*pfn2i)->r);
            double ta_angle = 180.0 / M_PI * acos(
                    ( (*dfn1i)*(*dfn1i) + (*dfn2i)*(*dfn2i) - dneib*dneib ) /
                    (2.0*(*dfn1i)*(*dfn2i))  );
            if (ta_angle < lo_bangle || ta_angle > hi_bangle)
            {
                return false;
            }
        }
        // check neighbor atom angles, na_angle
        LPAit psni = second_neighbors.begin();
        LDit dsni = second_distances.begin();
        for (; psni != second_neighbors.end(); ++psni, ++dsni)
        {
            double dneib = R3::distance((*pfn1i)->r, (*psni)->r);
            if (dneib == 0.0 || !(dneib < max_blen))
            {
                continue;
            }
            double na_angle = 180.0 / M_PI * acos(
                    ( (*dfn1i)*(*dfn1i) + dneib*dneib - (*dsni)*(*dsni) ) /
                    (2.0*(*dfn1i)*dneib)  );
            if (na_angle < lo_bangle || na_angle > hi_bangle)
            {
                return false;
            }
        }
    }
    return true;
}

void BondAngleFilter_t::setBondAngleRange(double lo, double hi)
{
    lo_bangle = lo;
    hi_bangle = hi;
}

////////////////////////////////////////////////////////////////////////
// class LoneAtomFilter_t - derived from AtomFilter_t
////////////////////////////////////////////////////////////////////////

// constructor

LoneAtomFilter_t::LoneAtomFilter_t(double maxbondlength)
{
    mmaxbondlength = maxbondlength;
}

// methods

bool LoneAtomFilter_t::Check(Atom_t* pta, Molecule* pm)
{
    AtomPtr apnear = pm->getNearestAtom(pta->r);
    // apnear holds no pointer for empty molecule.
    // the first atom is never alone.
    bool isnotalone = true;
    if (apnear.get())
    {
        double dnear = R3::distance(pta->r, apnear->r);
        isnotalone = eps_gt(dnear, 0.0) && (dnear < mmaxbondlength);
    }
    return isnotalone;
}

// End of file
