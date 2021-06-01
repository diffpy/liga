/***********************************************************************
* Short Title: functions for obtaining subversion information
*
* Comments: subversion information definitions.
*     DO NOT edit Version.cpp, change Version.tpl instead.
*
* <license text>
***********************************************************************/

#include "Version.hpp"

namespace NS_VERSION {

const string& getVersion()
{
    static const string srev = "${version}";
    return srev;
}

const string& getCommit()
{
    static const string srev = "${commit}";
    return srev;
}

const string& getDate()
{
    static const string sdate = "${date}";
    return sdate;
}

const string& getAuthor()
{
    static const string sauthor = "${author}";
    return sauthor;
}

const string& getId()
{
    static string sid;
    if (sid.empty())
    {
        sid = getVersion() + ' ' + getCommit().substr(0, 7) + ' ' +
            getDate() + ' ' + getAuthor();
    }
    return sid;
}

}   // namespace NS_VERSION

// vim:ft=cpp:
// End of Version.cpp
