/***********************************************************************
* Short Title: object declarations for Atom_t class
*
* Comments:
*
* <license text>
***********************************************************************/

#ifndef ATOM_T_HPP_INCLUDED
#define ATOM_T_HPP_INCLUDED

#include <string>
#include <boost/shared_ptr.hpp>
#include "R3linalg.hpp"

// may be later moved to dedicated triangulation class
enum triangulation_type { LINEAR, PLANAR, SPATIAL, NTGTYPES };

////////////////////////////////////////////////////////////////////////
// class Atom_t
////////////////////////////////////////////////////////////////////////

class Atom_t
{
    public:

        // friends
        friend class Molecule;
        friend class Crystal;

        // constructors
        Atom_t(const std::string& elsmbl, double rx, double ry, double rz);
        template <class V> Atom_t(const std::string& elsmbl, const V& rxyz);
        template <class V> Atom_t(Atom_t* asrc, const V& rxyz);

        // data
        std::string element;
        R3::Vector r;
        bool fixed;
        triangulation_type ttp;
        double radius;

        // methods
        const double& Badness() const;
        double FreeBadness() const;
        void IncBadness(const double& db);
        void DecBadness(const double& db);
        void ResetBadness(double b=0.0);
        const double& Overlap() const;
        double FreeOverlap() const;
        void IncOverlap(const double& doverlap);
        void DecOverlap(const double& doverlap);
        void ResetOverlap(double overlap=0.0);
        double costShare(double pairsperatom);


    private:

        // data
        double _badness;
        double _overlap;
        mutable int pmxidx;     // pair matrix index
        Atom_t* mstorage_ptr;

        // methods
        void init(double rx, double ry, double rz);
};


typedef boost::shared_ptr<Atom_t> AtomPtr;

// declarations of non-member operators and functions

bool operator==(const Atom_t& a1, const Atom_t& a2);

////////////////////////////////////////////////////////////////////////
// template and inline functions
////////////////////////////////////////////////////////////////////////

// template constructor

template <class V>
Atom_t::Atom_t(const std::string& elsmbl, const V& rxyz) :
    element(elsmbl), mstorage_ptr(NULL)
{
    this->init(double(rxyz[0]), double(rxyz[1]), double(rxyz[2]));
}


template <class V>
Atom_t::Atom_t(Atom_t* asrc, const V& rxyz) :
    element(asrc->element), mstorage_ptr(asrc)
{
    this->init(double(rxyz[0]), double(rxyz[1]), double(rxyz[2]));
}




#endif  // ATOM_T_HPP_INCLUDED
