/***********************************************************************
* Short Title: useful functions for string manipulation
*
* Comments:
*
* <license text>
***********************************************************************/

#ifndef STRINGUTILS_HPP_INCLUDED
#define STRINGUTILS_HPP_INCLUDED

#include <string>
#include <sstream>

template<typename StrType, typename Sequence>
std::string join(StrType sep, Sequence words)
{
    std::string joined;
    typename Sequence::iterator w = words.begin();
    typename Sequence::iterator last = words.end();
    if (w == last)      return joined;
    for (joined = *(w++); w != last; ++w)
    {
        joined += sep;
        joined += *w;
    }
    return joined;
}

template<typename StrType, typename Iterator>
std::string join(StrType sep, Iterator first, Iterator last)
{
    std::string joined;
    Iterator w = first;
    if (w == last)      return joined;
    for (joined = *(w++); w != last; ++w)
    {
        joined += sep;
        joined += *w;
    }
    return joined;
}

template<typename Container>
void split(const std::string& s, Container& words)
{
    std::istringstream istrs(s);
    std::string w;
    words.clear();
    while (istrs >> w)  words.push_back(w);
}


#endif  // STRINGUTILS_HPP_INCLUDED
