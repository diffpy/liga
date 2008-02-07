/***********************************************************************
* Short Title: functions for obtaining subversion information
*
* Comments: subversion information definitions
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef VERSION_HPP_INCLUDED
#define VERSION_HPP_INCLUDED

#include <string>

namespace NS_VERSION {

using namespace std;

int getRevisionNumber();
const string& getRevision();
const string& getDate();
const string& getAuthor();
const string& getId();

}   // namespace NS_VERSION

#endif	// VERSION_HPP_INCLUDED
