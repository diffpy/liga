/*****************************************************************************
*
* Liga Algorithm    for structure determination from pair distances
*                   Pavol Juhas
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* AtomRadiiTable
*
* Comments: Storage of empirical atomic radii
*
*****************************************************************************/

#ifndef ATOMRADIITABLE_HPP_INCLUDED
#define ATOMRADIITABLE_HPP_INCLUDED

#include <map>
#include <string>


class AtomRadiiTable : public std::map<std::string,double>
{
    public:

        /// fast value lookup, which does not change the table.
        double lookup(const std::string& smbl) const;
        /// initialize from a string in (A1:r1, A2:r2, ...) format
        void fromString(std::string);
        /// convert to a string in (A1:r1,A2:r2,...) format
        std::string toString(std::string separator=",") const;
};


#endif  // ATOMRADIITABLE_HPP_INCLUDED
