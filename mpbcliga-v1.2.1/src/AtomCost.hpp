/***********************************************************************
* Short Title: AtomCost - functor for external atom cost calculation
*
* Comments: base class for AtomCostCrystal
*
* <license text>
***********************************************************************/

#ifndef ATOMCOST_HPP_INCLUDED
#define ATOMCOST_HPP_INCLUDED

#include <vector>
#include "R3linalg.hpp"

class Molecule;
class Atom_t;

////////////////////////////////////////////////////////////////////////
// Declarations
////////////////////////////////////////////////////////////////////////

// functions

double penalty(const double& dd, const double& desd);
double penalty_gradient(const double& dd, const double& desd);

// classes

class AtomCost
{
    public:

        enum EvalFlag { NONE, GRADIENT=1, SELFCOST=2 };

        // constructor
        AtomCost(const Molecule* m);

        // destructor
        virtual ~AtomCost() { }

        // public methods
        virtual void resetFor(const Molecule* m);
        double eval(const Atom_t& pa, int flags=NONE);
        virtual double eval(const Atom_t* pa, int flags=NONE);
        const R3::Vector& gradient();
        double lowest() const;
        double cutoff() const;
        void setCutoff(double cf);
        double cutoffRange() const;
        void setCutoffRange(double cutrng);
        void noCutoff();
        double totalCost() const;
        const std::vector<double>& partialCosts() const;
        const std::vector<double>& targetDistances() const;
        const std::vector<int>& usedTargetDistanceIndices() const;
        const std::vector<int>& usedTargetAtomIndices() const;
        void setScale(double);
        const double& getScale() const;
        double penaltyScaled(const double& dd, const double& desd) const;

    protected:

        // data - arguments
        const Molecule* arg_cluster;
        const Atom_t* arg_atom;

        // data - results
        bool use_distances;
        std::vector<bool> useflag;
        bool apply_cutoff;
        double lowest_cost;
        double cutoff_cost;
        double cutoff_range;
        double total_cost;
        std::vector<double> partial_costs;
        std::vector<double> target_distances;
        std::vector<int> useflag_indices;
        std::vector<int> useatom_indices;

        // optimizer specific data
        bool _selfcost_flag;
        bool _gradient_flag;
        bool _gradient_cached;
        R3::Vector _gradient;

        // protected methods
        const std::vector<Atom_t*>& getClusterAtoms() const;
        virtual void resizeArrays();
        void resetUseFlags();
        void resetGradient();
        size_t nearDistanceIndex(const double& d) const;
        double nearDistance(const double& d) const;

    private:

        // data - penalty configuration
        double mscale;

};  // class AtomCost

#endif  // ATOMCOST_HPP_INCLUDED
