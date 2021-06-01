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

#include <cstdlib>
#include <cctype>
#include <sstream>
#include <cassert>
#include <map>

#include "ChemicalFormula.hpp"

using namespace std;

// Constructors --------------------------------------------------------------

ChemicalFormula::ChemicalFormula()
{
}


ChemicalFormula::ChemicalFormula(const string& formula)
{
    this->fromString(formula);
}

// Public Methods ------------------------------------------------------------

int ChemicalFormula::countElements() const
{
    int cnt = 0;
    for (const_iterator sc = this->begin(); sc != this->end(); ++sc)
    {
        cnt += max(0, sc->second);
    }
    return cnt;
}


vector<string> ChemicalFormula::expand() const
{
    vector<string> rv;
    for (const_iterator sc = this->begin(); sc != this->end(); ++sc)
    {
        int c = max(0, sc->second);
        rv.insert(rv.end(), c, sc->first);
    }
    assert(int(rv.size()) == this->countElements());
    return rv;
}


void ChemicalFormula::fromString(const std::string& s)
{
    ChemicalFormula chfm;
    for (string::const_iterator c = s.begin(); c != s.end(); ++c)
    {
        if (isblank(*c))  continue;
        if (chfm.empty() || isupper(*c))
        {
            chfm.push_back(ChemicalFormula::value_type("", 0));
        }
        string& smblcnt = chfm.back().first;
        smblcnt += *c;
    }
    for (iterator sc = chfm.begin(); sc != chfm.end(); ++sc)
    {
        string& smblcnt = sc->first;
        string::size_type pcnt;
        pcnt = smblcnt.find_last_not_of("0123456789") + 1;
        string scnt = smblcnt.substr(pcnt);
        smblcnt.erase(pcnt);
        sc->second = scnt.empty() ? 1 : atoi(scnt.c_str());
    }
    // assign new formula when everything went well
    this->swap(chfm);
}


string ChemicalFormula::toString(string separator) const
{
    vector<string> smbl_order;
    map<string,int> smbl_count;
    for (const_iterator sc = this->begin(); sc != this->end(); ++sc)
    {
        if (!smbl_count.count(sc->first) && sc->second > 0)
        {
            smbl_order.push_back(sc->first);
        }
        smbl_count[sc->first] += sc->second;
    }
    ostringstream rv;
    vector<string>::const_iterator smb;
    for (smb = smbl_order.begin(); smb != smbl_order.end(); ++smb)
    {
        if (smb != smbl_order.begin())  rv << separator;
        rv << *smb;
        int cnt = smbl_count[*smb];
        if (cnt > 1)  rv << cnt;
    }
    return rv.str();
}

// End of file
