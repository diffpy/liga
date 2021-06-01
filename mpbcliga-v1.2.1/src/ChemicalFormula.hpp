/***********************************************************************
*
* Liga Algorithm    for structure determination from pair distances
*                   Pavol Juhas
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* ChemicalFormula
*
* Comments: Storage and manipulation of ChemicalFormula.
*
***********************************************************************/

#ifndef CHEMICALFORMULA_HPP_INCLUDED
#define CHEMICALFORMULA_HPP_INCLUDED

#include <utility>
#include <vector>
#include <string>

class ChemicalFormula : public std::vector< std::pair<std::string,int> >
{
    public:

        // Constructors
        ChemicalFormula();
        ChemicalFormula(const std::string&);

        // Methods

        /// total count of elements
        int countElements() const;
        /// expand to a vector of elements
        std::vector<std::string> expand() const;
        /// initialize from a string
        void fromString(const std::string&);
        /// convert to a string
        std::string toString(std::string separator="") const;

};

#endif  // CHEMICALFORMULA_HPP_INCLUDED
