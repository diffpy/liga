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

#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <cctype>
#include <sstream>

#include <boost/foreach.hpp>

#include "AtomRadiiTable.hpp"
#include "StringUtils.hpp"

using namespace std;

// Public Methods ------------------------------------------------------------

double AtomRadiiTable::lookup(const string& smbl) const
{
    const_iterator ti = this->find(smbl);
    if (ti == this->end())
    {
        ostringstream emsg;
        emsg << "Undefined radius for element '" << smbl << "'.";
        throw invalid_argument(emsg.str());
    }
    return ti->second;
}


void AtomRadiiTable::fromString(string s)
{
    // remove any blank characters
    s.erase(remove_if(s.begin(), s.end(), isblank), s.end());
    // replace commas with space so we can use the split function
    replace(s.begin(), s.end(), ',', ' ');
    vector<string> words;
    split(s, words);
    // perform parsing on a new table and only assign when everything works
    AtomRadiiTable tnew;
    BOOST_FOREACH (string w, words)
    {
        string::size_type p = w.find(':');
        if (p == string::npos)
        {
            ostringstream emsg;
            emsg << "Invalid radius specification, missing ':' in '" <<
                w << "'.";
            throw invalid_argument(emsg.str());
        }
        string smbl = w.substr(0, p);
        double value;
        istringstream sv(w.substr(p + 1));
        if (!(sv >> value))
        {
            ostringstream emsg;
            emsg << "Invalid floating point number in '" << w << "'.";
            throw invalid_argument(emsg.str());
        }
        tnew[smbl] = value;
    }
    // everything worked up to here, we can do the assignment
    *this = tnew;
}


string AtomRadiiTable::toString(string separator) const
{
    ostringstream rv;
    for (const_iterator ti = this->begin(); ti != this->end(); ++ti)
    {
        if (ti != this->begin())  rv << separator;
        rv << ti->first << ':' << ti->second;
    }
    return rv.str();
}

// End of file
