/*****************************************************************************
* Short Title: class DistanceTable - declaration
*
* Comments:
*
* <license text>
*****************************************************************************/

#ifndef DISTANCETABLE_HPP_INCLUDED
#define DISTANCETABLE_HPP_INCLUDED

#include <vector>
#include <iostream>
#include <boost/unordered_map.hpp>

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
        const double& getesd(const double& d) const;
        void setESDs(const std::vector<double>& esds);
        void clearESDs();
        bool hasESDs() const;
        int estNumAtoms() const;
        int countUnique() const;
        DistanceTable unique() const;
        double getResolution() const;
        void setResolution(double res);
        double maxDistance() const;
        double maxDistanceRepr() const;

    private:

        // data
        mutable int mcount_unique;
        mutable double mmaxdistancerepr;
        double mresolution;
        boost::unordered_map<double,double> mesd;

        // methods
        void init();
        void readESDFormat(std::istream&);
        void readPWAFormat(std::istream&);
        void readSimpleFormat(std::istream&);
};

// non-member operators
std::istream& operator>>(std::istream&, DistanceTable&);

#endif  // DISTANCETABLE_HPP_INCLUDED
