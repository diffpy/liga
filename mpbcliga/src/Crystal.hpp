/***********************************************************************
* Short Title: derived class from Molecule for periodic crystal structure
*
* Comments:
*
* <license text>
***********************************************************************/

#ifndef CRYSTAL_T_HPP_INCLUDED
#define CRYSTAL_T_HPP_INCLUDED

#include <boost/shared_ptr.hpp>
#include "Molecule.hpp"

class Lattice;
class AtomCostCrystal;

class Crystal : public Molecule
{
    public:

        // friends
        friend class AtomCostCrystal;

        // data

        // constructors
        Crystal();
        Crystal(const Crystal&);

        // operators
        virtual Crystal& operator=(const Molecule&);
        Crystal& operator=(const Crystal&);
        virtual Molecule* copy() const;     // create a copy

        // methods - class registration and type info
        virtual StructureType type() const  { return CRYSTAL; }
        virtual std::string typeStr() const { return "crystal"; }

        // methods
        virtual void setDistanceTable(const DistanceTable&);
        virtual void setDistReuse(bool);

        void setLattice(const Lattice&);
        const Lattice& getLattice() const;

        void setRmax(double rmax);
        double getRmax() const;
        std::pair<double,double> getRExtent(double rlo, double rhi) const;

        virtual double cost() const;
        virtual int countPairs() const;
        virtual void recalculate() const;
        virtual AtomCost* getAtomCostCalculator() const;
        virtual AtomCost* getAtomOverlapCalculator() const;

        virtual void Shift(const R3::Vector& drc);

        virtual AtomPtr getNearestAtom(const R3::Vector& rc) const;
        virtual void Clear();
        virtual const std::pair<int*,int*>& Evolve(const int* est_triang);
        virtual void Degenerate(int Npop, DegenerateFlags flags=NONE);

    protected:

        // data
        boost::shared_ptr<DistanceTable> _full_distance_table;
        mutable SymmetricMatrix<int> pmx_pair_counts;

        // methods
        virtual void AddInternal(Atom_t* pa);  // add atom from the storage
        virtual void addNewAtomPairs(Atom_t* pa);
        virtual void removeAtomPairs(Atom_t* pa);
        virtual const TriangulationAnchor&
            getLineAnchor(const RandomWeighedGenerator& rwg);
        virtual const TriangulationAnchor&
            getPlaneAnchor(const RandomWeighedGenerator& rwg);
        virtual const TriangulationAnchor&
            getPyramidAnchor(const RandomWeighedGenerator& rwg);
        virtual void resizePairMatrices(int sz);
        virtual boost::python::object newDiffPyStructure() const;
        virtual void setFromDiffPyStructure(boost::python::object);

    private:

        // crystal specific data
        boost::shared_ptr<const Lattice> _lattice;
        double _rmin;
        double _rmax;
        mutable double _cost;
        mutable int _count_pairs;
        mutable bool _cost_data_cached;
        double _lattice_max_ucd;

        // class methods
        boost::shared_ptr<Lattice> getDefaultLattice();

        // methods
        void init();
        void uncacheCostData() const;
        void cropDistanceTable();
        const R3::Vector&
            anyOffsetAtomSite(const RandomWeighedGenerator& rwg) const;
        R3::Vector ucvCartesianAdjusted(const R3::Vector& cv) const;
        void shiftToOrigin();
};


#endif  // CRYSTAL_T_HPP_INCLUDED
