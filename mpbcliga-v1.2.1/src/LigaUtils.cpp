#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include "LigaUtils.hpp"
#include "Exceptions.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////////

// valarray operations

double vdnorm(const valarray<double>& v)
{
    double s2 = 0.0;
    for (const double *pv = &v[0]; pv != &v[v.size()]; ++pv)
        s2 += (*pv)*(*pv);
    return sqrt(s2);
}

double vddot(const valarray<double>& v1, const valarray<double>& v2)
{
    // assuming v1.size() == v2.size()
    double dotprod = 0.0;
    const double *pv1 = &v1[0], *pv2 = &v2[0];
    for (; pv1 != &v1[v1.size()]; ++pv1, ++pv2)
        dotprod += (*pv1)*(*pv2);
    return dotprod;
}

valarray<double> vdcross(const valarray<double>& v1, const valarray<double>& v2)
{
    if (v1.size() != 3 || v2.size() != 3)
    {
        throw runtime_error("vdcross(): invalid valarray size");
    }
    valarray<double> cross(3);
    cross[0] = v1[1]*v2[2] - v1[2]*v2[1];
    cross[1] = v1[2]*v2[0] - v1[0]*v2[2];
    cross[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return cross;
}

// file utilities

ofstream& mktempofstream(ofstream& out, string& writefile)
{
    const int max_attempts = 10;
    string::size_type pX = writefile.find("XXXXXX");
    if (pX == string::npos)
    {
        throw IOError("mktempofstream(): XXXXXX missing in filename");
    }
    const string alnum = "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    for (int i = 0; ; ++i)
    {
        // create random name
        for (int j = 0; j < 6; ++j)
        {
            writefile[pX + j] = alnum[rand() % alnum.size()];
        }
        ifstream fid(writefile.c_str());
        if (!fid)   break;
        else        fid.close();
        if (i >= max_attempts)
        {
            const char* emsg = "mktempofstream(): max_attempts exceeded.";
            throw IOError(emsg);
        }
    }
    // we should have a good name here
    out.open(writefile.c_str());
    if (!out)
    {
        ostringstream emsg;
        emsg << "mktempofstream(): unable to open " << writefile <<
            "for writing.";
        throw IOError(emsg.str());
    }
    return out;
}

// read lines that do not start with number
bool read_header(istream& fid, string& header)
{
    double x;
    string line;
    istringstream istrs;
    header.clear();
    for (   istream::pos_type p = fid.tellg();
            !fid.eof() && getline(fid, line);
            p = fid.tellg()
        )
    {
        istrs.clear();
        istrs.str(line);
        if (istrs >> x)
        {
            fid.seekg(p);
            break;
        }
        else
        {
            header.append(line + '\n');
        }
    }
    return !(fid.rdstate() & ios::badbit);
}

bool read_header(istream& fid)
{
    string dummy;
    return read_header(fid, dummy);
}

// End of file
