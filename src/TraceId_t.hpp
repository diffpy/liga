/***********************************************************************
*
* Liga Algorithm    for structure determination from pair distances
*                   Pavol Juhas
*                   (c) 2007 trustees of the Michigan State University
*                   All rights reserved.
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* struct TraceId_t
*
* Comments: Data storage for a single trace point of Liga competitor.
*
***********************************************************************/

#ifndef TRACEID_T_HPP_INCLUDED
#define TRACEID_T_HPP_INCLUDED

struct TraceId_t
{
    int season;                 // season when molecule was modified
    int level;                  // liga level played when change occurred
    size_t mol_natoms;          // molecule size after change
    double mol_norm_badness;    // molecule cost after change
    long mol_id;                // unique identifier of affected molecule
};

#endif  // TRACEID_T_HPP_INCLUDED
