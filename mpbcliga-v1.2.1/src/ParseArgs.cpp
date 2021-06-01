/*****************************************************************************
* Short Title: program arguments parser
*
* Comments:
*
* <license text>
*****************************************************************************/

#include <algorithm>
#include <fstream>
#include <bitset>

#include "Exceptions.hpp"
#include "ParseArgs.hpp"

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// class ParseArgs
//////////////////////////////////////////////////////////////////////////////

ParseArgs::ParseArgs(int nargc, char* const nargv[]) :
    argc(nargc), argv(nargv)
{
    optstring = NULL;
    longopts = NULL;
    init();
}

ParseArgs::ParseArgs(int nargc, char* const nargv[], const char *noptstring) :
    argc(nargc), argv(nargv), optstring(noptstring)
{
    longopts = NULL;
    init();
}

ParseArgs::ParseArgs(int nargc, char* const nargv[], const char *noptstring,
            const struct option *nlongopts) :
    argc(nargc), argv(nargv), optstring(noptstring), longopts(nlongopts)
{
    init();
}

void ParseArgs::Parse()
{
    if (longopts && optstring)
        do_getopt_long();
    else if (optstring)
        do_getopt();
    else
    {
        for (int i = 1; i < argc; ++i)
            arg_or_par(argv[i]);
    }
}

void ParseArgs::Dump()
{
    cout << "cmd_h = '" << cmd_h << "'" << endl;
    cout << "cmd_t = '" << cmd_t << "'" << endl;
    typedef map<string,string>::iterator MSSit;
    for (MSSit ii = opts.begin(); ii != opts.end(); ++ii)
        cout << "opts[" << ii->first << "] = '" << ii->second << "'" << endl;
    for (MSSit ii = pars.begin(); ii != pars.end(); ++ii)
        cout << "pars[" << ii->first << "] = '" << ii->second << "'" << endl;
    for (size_t i = 0; i < args.size(); ++i)
        cout << "args[" << i << "] = '" << args[i] << "'" << endl;
}

vector<int> ParseArgs::ExpandRangePar(string par)
{
    if (!pars.count(par))
    {
        ostringstream emsg_ostream;
        emsg_ostream << "parameter '" << par << "' is not defined";
        throw ParseArgsError(emsg_ostream.str());
    }
    // replace all commas in par with <space>
    string ranges(pars[par]);
    for (   string::size_type pcomma = ranges.find(',');
            pcomma != string::npos; pcomma = ranges.find(',', pcomma) )
    {
        ranges[pcomma] = ' ';
    }
    vector<int> indices;
    istringstream iss(ranges);
    string word;
    while (iss >> word)
    {
        // prepare error message
        ostringstream emsg_ostream;
        emsg_ostream << par << ": invalid range value '"
            << word << "'";
        // process word
        string::size_type prange = word.find("..");
        if (prange == string::npos)
        {
            istringstream word_istream(word);
            int idx;
            if ( !(word_istream >> idx) )
            {
                throw ParseArgsError(emsg_ostream.str());
            }
            indices.push_back(idx);
        }
        else
        {
            word.replace(prange, 2, " ");
            istringstream word_istream(word);
            int start, stop;
            if ( !(word_istream >> start >> stop) )
            {
                throw ParseArgsError(emsg_ostream.str());
            }
            for (int idx = start; idx < stop+1; ++idx)
            {
                indices.push_back(idx);
            }
        }
    }
    return indices;
}

void ParseArgs::ValidatePars(const list<string>& validpars)
{
    list<string> invalidpars;
    typedef map<string,string>::iterator MSSit;
    typedef list<string>::const_iterator LSit;
    for (MSSit ii = pars.begin(); ii != pars.end(); ++ii)
    {
        LSit ivp = find(validpars.begin(), validpars.end(), ii->first);
        if (ivp == validpars.end())
            invalidpars.push_back(ii->first);
    }
    if (invalidpars.size())
    {
        ostringstream emsg_ostream;
        if (invalidpars.size() == 1)
            emsg_ostream << "invalid parameter -- ";
        else
            emsg_ostream << "invalid parameters -- ";
        LSit ii = invalidpars.begin();
        emsg_ostream << *ii;
        for (++ii; ii != invalidpars.end(); ++ii)
            emsg_ostream << ", " << *ii;
        throw ParseArgsError(emsg_ostream.str());
    }
}

void ParseArgs::defParameterAlias(const string& a, const string& par)
{
    par_alias[a] = par;
}

string ParseArgs::expandParameterAlias(const string& p)
{
    string pe = p;
    if (par_alias.count(p))
    {
        pe = par_alias[p];
        if (!count(used_par_aliases.begin(), used_par_aliases.end(), p))
        {
            used_par_aliases.push_back(p);
        }
    }
    return pe;
}

const list<string>& ParseArgs::usedParameterAliases() const
{
    return used_par_aliases;
}

void ParseArgs::ReadPars(const string& parfile)
{
    // open parfile for reading
    ifstream fid(parfile.c_str());
    if (!fid)
    {
        ostringstream emsg_ostream;
        emsg_ostream << "Unable to read '" << parfile << "'";
        throw IOError(emsg_ostream.str());
    }
    try {
        ReadPars(fid);
    }
    catch (ParseArgsError(e)) {
        ostringstream emsg;
        emsg << "invalid syntax in parameter file\n" <<
            parfile << ":" << e.what();
        throw ParseArgsError(emsg.str());
    }
    fid.close();
}

bool check_backslash(string& s)
{
    string::reverse_iterator rii;
    // count backslashes at the end of string
    int nbs = 0;
    for (rii = s.rbegin(); rii != s.rend() && *rii == '\\'; ++rii)
        ++nbs;
    bool line_continues = (nbs % 2 != 0);
    if (nbs > 0)
    {
        string::size_type pend;
        pend = s.size() - nbs/2 - (line_continues);
        s.erase(pend);
    }
    return line_continues;
}

istream& ParseArgs::ReadPars(istream& fid)
{
    const char* blank = " \t\r\n";
    // Parameters reading will finish after long hash line and
    // blank line at first non-parameter line.
    enum TermFlags { HASHLINE, BLANKLINE, NUMCOND };
    bitset<NUMCOND> term_conditions(0);
    string pline, fline;
    bool line_continues = false;
    for (int nr = 1; !fid.eof() && getline(fid, fline); ++nr)
    {
        line_continues = check_backslash(fline);
        pline = line_continues ? (pline + fline) : fline;
        if (line_continues)    continue;
        // here we should have complete pline
        string::size_type lb, le;
        lb = pline.find_first_not_of(blank);
        le = pline.find_last_not_of(blank)+1;
        pline = (lb == string::npos) ? "" : pline.substr(lb, le-lb);
        // ship lines that start with hash mark
        if (pline.length() && pline[0] == '#')
        {
            // check for terminating hashline
            string::size_type hashcount = pline.find_first_not_of('#');
            term_conditions[HASHLINE] = (hashcount > 64);
            continue;
        }
        // skip blank line
        if (pline.empty())
        {
            term_conditions[BLANKLINE] = term_conditions[HASHLINE];
            continue;
        }
        string::size_type eq, pe, vb;
        eq = pline.find('=');
        if (eq == string::npos)
        {
            // get out when term_conditions are all true
            if (term_conditions.count() == term_conditions.size())  break;
            ostringstream emsg_ostream;
            emsg_ostream << nr << ": missing equal symbol";
            throw ParseArgsError(emsg_ostream.str());
        }
        // here we have something that looks like parameter, keep reading
        term_conditions.reset();
        pe = pline.find_last_not_of(string(blank) + "=", eq);
        pe = (pe == string::npos) ? eq : pe+1;
        string par0 = pline.substr(0, pe);
        string par = expandParameterAlias(par0);
        bool ispar = isalpha(par[0]);
        for (   string::iterator ii = par.begin();
                ii != par.begin() && ispar; ++ii)
        {
            ispar = isalnum(*ii) || *ii == '_';
        }
        if (!ispar)
        {
            ostringstream emsg_ostream;
            emsg_ostream << nr << ": invalid parameter name '" << par << "'";
            throw ParseArgsError(emsg_ostream.str());
        }
        if (cmdl_par.count(par))    continue;
        vb = pline.find_first_not_of(blank, eq+1);
        pars[par] = (vb == string::npos) ? "" : pline.substr(vb);
    }
    return fid;
}

void ParseArgs::init()
{
    cmd_h = argv[0];
    string::size_type pslash = cmd_h.find_last_of('/');
    if (pslash != string::npos)
    {
        cmd_t = cmd_h.substr(pslash+1);
        cmd_h.erase(pslash);
    }
    else
        cmd_h.swap(cmd_t);
}

void ParseArgs::do_getopt()
{
    while (true)
    {
        int c = getopt(argc, argv, optstring);
        if (c == -1)
            break;
        else if (c == '?')
            throw ParseArgsError();
        else
            opts[string(1,c)] = (optarg) ? optarg : "";
    }
    for (; optind < argc; ++optind)
    {
        arg_or_par(argv[optind]);
    }
}

void ParseArgs::do_getopt_long()
{
    while (true)
    {
        int option_index = 0;
        int c = getopt_long(argc, argv, optstring,
                longopts, &option_index);
        if (c == -1)
            break;
        string o;
        switch (c)
        {
            case '?':
                throw ParseArgsError();
                break;
            case 0:
            case '-':
                o = longopts[option_index].name;
                opts[o] = (optarg) ? optarg : "";
                break;
            default:
                o = c;
                opts[o] = (optarg) ? optarg : "";
        }
    }
    for (; optind < argc; ++optind)
    {
        arg_or_par(argv[optind]);
    }
}

void ParseArgs::arg_or_par(const char *s)
{
    string a(s);
    string::size_type peq = a.find('=');
    if (peq == string::npos)
    {
        args.push_back(a);
        return;
    }
    // here it looks like a parameter
    string pname = a.substr(0, peq);
    bool ispar = true;
    for (string::iterator p = pname.begin(); ispar && p != pname.end(); ++p)
    {
        ispar = isalnum(*p) || *p == '_';
    }
    if (ispar)
    {
        string parname = expandParameterAlias(pname);
        pars[parname] = a.substr(peq + 1);
        cmdl_par.insert(parname);
    }
    else
    {
        args.push_back(a);
    }
}


void ParseArgs::force_valid_par(const string& par) const
{
    if (!this->pars.count(par))
    {
        ostringstream emsg;
        emsg << "parameter '" << par << "' is not defined";
        throw ParseArgsError(emsg.str());
    }
}

// End of file
