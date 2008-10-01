/***********************************************************************
* Short Title: class DistanceTable - declaration
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef DISTANCETABLE_HPP_INCLUDED
#define DISTANCETABLE_HPP_INCLUDED

#include <vector>
#include <iostream>

class DistanceTable : public std::vector<double>
{
    public:

        // friends
	friend std::istream& operator>>(std::istream&, DistanceTable&);

        // constructors
        DistanceTable();
        DistanceTable(const double* v, size_t sz);
        DistanceTable(const std::vector<double>&);
        DistanceTable(const DistanceTable&);

        // methods
        DistanceTable& operator= (const std::vector<double>&);
        DistanceTable& operator= (const DistanceTable&);
        const_iterator find_nearest(const double& d) const;
        iterator return_back(const double&);
        int estNumAtoms() const;
        int countUnique() const;
        std::vector<double> unique() const;
        double getResolution() const;
        void setResolution(double res);
        double maxDistance() const;

    private:

        // data
        mutable int count_unique;
        double resolution;

        // methods
        void init();
        void readPWAFormat(std::istream&);
        void readSimpleFormat(std::istream&);
};

// non-member operators
std::istream& operator>>(std::istream&, DistanceTable&);

#endif	// DISTANCETABLE_HPP_INCLUDED
