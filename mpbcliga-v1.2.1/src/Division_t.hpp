/*****************************************************************************
* Short Title: one division of the liga system, declarations
*
* Comments:
*
* <license text>
*****************************************************************************/

#ifndef DIVISION_T_HPP_INCLUDED
#define DIVISION_T_HPP_INCLUDED

#include <vector>
#include "Atom_t.hpp"   // for NTGTYPES
#include "Random.hpp"

class Molecule;

class Division_t : public std::vector<Molecule*>
{
    public:

        // friends
        friend class Liga_t;

        // types
        typedef Molecule* PMOL;

        // constructors and destructor
        Division_t(size_t fullsize, size_t level);
        Division_t(const Division_t& src);
        ~Division_t();

        // operators
        Division_t& operator= (const Division_t& src);

        // public methods
        int find_winner();
        int find_looser();
        int find_best();
        PMOL& best();
        bool full() const;
        size_t fullsize() const;
        size_t level() const;
        void assignTrials(double t);
        double trials();
        double averageCost() const;
        const int* estimateTriangulations();
        void noteTriangulations(const std::pair<int*,int*>& acc_tot);

    private:

        // class data
        static size_t ndim;

        // data members
        size_t _fullsize;
        size_t _level;
        double _trials;
        long long acc_triang[NTGTYPES];
        long long tot_triang[NTGTYPES];
        int est_triang[NTGTYPES];

};

// Inline Definitions --------------------------------------------------------

inline bool Division_t::full() const
{
    return !(size() < _fullsize);
}


inline size_t Division_t::fullsize() const
{
    return _fullsize;
}


inline size_t Division_t::level() const
{
    return _level;
}


inline void Division_t::assignTrials(double t)
{
    _trials = t;
}


inline double Division_t::trials()
{
    return _trials;
}

#endif  // DIVISION_T_HPP_INCLUDED
