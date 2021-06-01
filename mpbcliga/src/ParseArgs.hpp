/*****************************************************************************
* Short Title: program arguments parser
*
* Comments:
*
* <license text>
*****************************************************************************/

#ifndef PARSEARGS_HPP_INCLUDED
#define PARSEARGS_HPP_INCLUDED
#include <iostream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <stdexcept>
#include <getopt.h>

class ParseArgsError : public std::runtime_error
{
    public:
        ParseArgsError (std::string what_arg="") :
            std::runtime_error(what_arg)
        { }
};

class ParseArgs
{
public:
    // constructor
    ParseArgs(int nargc, char * const nargv[]);
    ParseArgs(int nargc, char * const nargv[], const char *optstring);
    ParseArgs(int nargc, char * const nargv[], const char *optstring,
            const struct option *longopts);
    int argc;
    char * const * argv;
    const char* optstring;
    const struct option *longopts;
    std::map<std::string,std::string> opts;
    std::map<std::string,std::string> pars;
    std::vector<std::string> args;
    std::string cmd_h, cmd_t;
    void Parse();
    void Dump();
    void ReadPars(const std::string& parfile);
    std::istream& ReadPars(std::istream& fid=std::cin);
    template<typename T> T GetPar(std::string par);
    std::vector<int> ExpandRangePar(std::string par);
    template<typename T> T GetPar(std::string par, T defval);
    template<typename T> std::vector<T> GetParVec(std::string par);
    void ValidatePars(const std::list<std::string>& validpars);
    void defParameterAlias(const std::string& a, const std::string& par);
    std::string expandParameterAlias(const std::string& p);
    const std::list<std::string>& usedParameterAliases() const;
    inline bool ispar(const std::string& p) { return pars.count(p); }
    inline bool isopt(const std::string& o) { return opts.count(o); }
private:
    void init();
    void do_getopt();
    void do_getopt_long();
    void arg_or_par(const char *s);
    void force_valid_par(const std::string&p) const;
    std::set<std::string> cmdl_par;
    std::map<std::string,std::string> par_alias;
    std::list<std::string> used_par_aliases;
};

template<typename T> T ParseArgs::GetPar(std::string par)
{
    using namespace std;
    this->force_valid_par(par);
    T val;
    istringstream iss(pars[par]);
    iss >> val || (iss.clear(), iss >> boolalpha >> val);
    if (!iss)
    {
        ostringstream oss;
        oss << "invalid value for parameter '" << par << "'";
        throw ParseArgsError(oss.str());
    }
    return val;
}

template<typename T> T ParseArgs::GetPar(std::string par, T defval)
{
    return pars.count(par) ? GetPar<T>(par) : defval;
}

template<typename T> std::vector<T> ParseArgs::GetParVec(std::string par)
{
    using namespace std;
    // replace all commas in par with <space>
    string values(pars[par]);
    for (   string::size_type pcomma = values.find(',');
            pcomma != string::npos; pcomma = values.find(',', pcomma) )
    {
        values[pcomma] = ' ';
    }
    vector<T> v;
    istringstream iss(values);
    T val;
    while (iss >> val)
        v.push_back(val);
    return v;
}

#endif
