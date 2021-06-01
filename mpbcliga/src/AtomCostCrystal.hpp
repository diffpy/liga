/***********************************************************************
* Short Title: AtomCostCrystal - functor for cost calculation with PBC.
*
* Comments:
*
* <license text>
***********************************************************************/

#ifndef ATOMCOSTCRYSTAL_HPP_INCLUDED
#define ATOMCOSTCRYSTAL_HPP_INCLUDED

#include <vector>
#include <memory>
#include "AtomCost.hpp"
#include "R3linalg.hpp"
#include "PointsInSphere.hpp"

class Crystal;
class Atom_t;

class AtomCostCrystal : public AtomCost
{
    public:

        // constructor
        AtomCostCrystal(const Crystal* m);

        // destructor
        virtual ~AtomCostCrystal() { }

        // public methods - overloaded
        virtual void resetFor(const Molecule* clust);
        double eval(const Atom_t& a, int flags=NONE);
        virtual double eval(const Atom_t* pa, int flags=NONE);
        int totalPairCount() const;
        const std::vector<int>& pairCounts() const;

        // public methods - specific
        std::pair<double,int> pairCostCount(const R3::Vector& cv);

    protected:

        // data - results
        int total_pair_count;
        std::vector<int> pair_counts;

        // data - arguments
        const Crystal* arg_cluster;
        R3::Vector arg_rcuc;    // cartesian positions offset to unit cell
        const Atom_t* crst_atom;

        // data - for intermediate cost evaluation
        // maximum r for PDF range
        double _rmax;
        // lattice points sequencer
        std::auto_ptr<PointsInSphere> _sph;

        // methods
        virtual const std::pair<double,double>&
            pairDistanceDifference(const double& d) const;
        void resizeArrays();

};  // class AtomCostCrystal

#endif  // ATOMCOSTCRYSTAL_HPP_INCLUDED
