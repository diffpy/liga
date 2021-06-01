/***********************************************************************
* Short Title: definitions for Atom_t class
*
* Comments:
*
* <license text>
***********************************************************************/

#include "LigaUtils.hpp"
#include "Atom_t.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////
// class Atom_t
////////////////////////////////////////////////////////////////////////

// constructors

Atom_t::Atom_t(const string& elsmbl, double rx, double ry, double rz) :
    element(elsmbl), mstorage_ptr(NULL)
{
    this->init(rx, ry, rz);
}

// public methods

const double& Atom_t::Badness() const
{
    return this->_badness;
}

double Atom_t::FreeBadness() const
{
    double rv = this->fixed ? 0.0 : this->Badness();
    return rv;
}

void Atom_t::IncBadness(const double& db)
{
    this->_badness += db;
    // Take care of round-offs, but only if they are very small.
    if (db < 0.0 && isNearZeroRoundOff(this->_badness))
    {
        this->ResetBadness(0.0);
    }
}

void Atom_t::DecBadness(const double& db)
{
    this->IncBadness(-db);
}

void Atom_t::ResetBadness(double b)
{
    this->_badness = b;
}

const double& Atom_t::Overlap() const
{
    return this->_overlap;
}

double Atom_t::FreeOverlap() const
{
    double rv = this->fixed ? 0.0 : this->Overlap();
    return rv;
}

void Atom_t::IncOverlap(const double& doverlap)
{
    this->_overlap += doverlap;
    // Take care of round-offs, but only if they are very small.
    if (doverlap < 0.0 && isNearZeroRoundOff(this->_overlap))
    {
        this->ResetOverlap(0.0);
    }
}

void Atom_t::DecOverlap(const double& doverlap)
{
    this->IncOverlap(-doverlap);
}

void Atom_t::ResetOverlap(double overlap)
{
    this->_overlap = overlap;
}


double Atom_t::costShare(double pairsperatom)
{
    double c = pairsperatom * this->Overlap() + this->Badness();
    return c;
}

// private methods

void Atom_t::init(double rx, double ry, double rz)
{
    this->r = rx, ry, rz;
    this->fixed = false;
    this->ttp = LINEAR;
    this->radius = 0.0;
    this->ResetBadness(0.0);
    this->ResetOverlap(0.0);
}

// non-member operators and functions

bool operator==(const Atom_t& a1, const Atom_t& a2)
{
    bool rv = (a1.element == a2.element &&
            equal(a1.r.data(), a1.r.data() + 3, a2.r.data()));
    return rv;
}

// End of file
