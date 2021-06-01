/*****************************************************************************
* Short Title: class DistanceTable - definitions
*
* Comments:
*
* <license text>
*****************************************************************************/

#include <sstream>
#include <cmath>
#include <cassert>

#include "DistanceTable.hpp"
#include "Exceptions.hpp"
#include "LigaUtils.hpp"
#include "StringUtils.hpp"

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// class DistanceTable
//////////////////////////////////////////////////////////////////////////////

// Constructors --------------------------------------------------------------

DistanceTable::DistanceTable() : vector<double>()
{
    init();
}


DistanceTable::DistanceTable(const double* v, size_t sz) :
    vector<double>(v, v + sz)
{
    init();
}


DistanceTable::DistanceTable(const vector<double>& v) : vector<double>(v)
{
    init();
}


DistanceTable::DistanceTable(const DistanceTable& d0)
{
    *this = d0;
}

// Operators -----------------------------------------------------------------

DistanceTable& DistanceTable::operator= (const vector<double>& v)
{
    vector<double>& this_vector = *this;
    this_vector = v;
    init();
    return *this;
}


DistanceTable& DistanceTable::operator= (const DistanceTable& d0)
{
    if (this == &d0)    return *this;
    assign(d0.begin(), d0.end());
    mcount_unique = d0.mcount_unique;
    mresolution = d0.mresolution;
    mesd = d0.mesd;
    return *this;
}

// Public Methods ------------------------------------------------------------

DistanceTable::const_iterator
DistanceTable::find_nearest(const double& dfind) const
{
    const_iterator ii = lower_bound(begin(), end(), dfind);
    bool shiftleft = (ii == end() && !empty()) ||
            (ii != begin() && (dfind - *(ii-1)) < (*ii - dfind));
    if (shiftleft)  --ii;
    return ii;
}


DistanceTable::iterator DistanceTable::return_back(const double& dback)
{
    iterator ii = lower_bound(begin(), end(), dback);
    return insert(ii, dback);
}


const double& DistanceTable::getesd(const double& d) const
{
    static const double one = 1.0;
    return mesd.empty() ? one : mesd.at(d);
}


void DistanceTable::setESDs(const std::vector<double>& esds)
{
    if (this->size() != esds.size())
    {
        const char* emsg = "Incompatible length of the ESD vector.";
        throw std::invalid_argument(emsg);
    }
    this->clearESDs();
    for (size_t i = 0; i != this->size(); ++i)
    {
        mesd[this->at(i)] = esds[i];
    }
    return;
}


void DistanceTable::clearESDs()
{
    mesd.clear();
}


bool DistanceTable::hasESDs() const
{
    return !mesd.empty();
}


int DistanceTable::estNumAtoms() const
{
    double xn = ( 1.0 + sqrt(1.0 + 8.0*this->size()) )/2.0;
    return int(xn);
}


int DistanceTable::countUnique() const
{
    if (mcount_unique >= 0)  return mcount_unique;
    mcount_unique = empty() ? 0 : 1;
    const_iterator ilo = begin();
    const_iterator ihi = (ilo == end()) ? ilo : ilo + 1;
    for (; ihi != end(); ++ilo, ++ihi)
    {
        if ( (*ihi - *ilo) > mresolution )   ++mcount_unique;
    }
    return mcount_unique;
}


DistanceTable DistanceTable::unique() const
{
    DistanceTable dtu;
    dtu.reserve(this->countUnique());
    vector<double> uesds;
    uesds.reserve(this->countUnique());
    for (const_iterator di = begin(); di != end(); ++di)
    {
        if (dtu.empty() || (*di - dtu.back()) > mresolution)
        {
            dtu.push_back(*di);
            uesds.push_back(this->getesd(*di));
        }
    }
    if (this->hasESDs())  dtu.setESDs(uesds);
    return dtu;
}


double DistanceTable::getResolution() const
{
    return mresolution;
}


void DistanceTable::setResolution(double res)
{
    mresolution = res;
    mcount_unique = -1;
}


double DistanceTable::maxDistance() const
{
    double mxd;
    mxd = this->empty() ? 0.0 : this->back();
    return mxd;
}


double DistanceTable::maxDistanceRepr() const
{
    double dr = mmaxdistancerepr - this->maxDistance();
    if (dr < 0.0 || dr > 1e-6)
    {
        for (int digits = int(1e8); digits >= 1; digits /= 10)
        {
            double mxdroundup = ceil(this->maxDistance() * digits) / digits;
            ostringstream omxd;
            omxd << mxdroundup;
            istringstream imxd(omxd.str());
            imxd >> mmaxdistancerepr;
            dr = mmaxdistancerepr - this->maxDistance();
            if (dr >= 0.0)  break;
        }
        assert(dr >= 0.0);
    }
    return mmaxdistancerepr;
}

// Private Methods -----------------------------------------------------------

void DistanceTable::init()
{
    mcount_unique = -1;
    mresolution = NS_LIGA::eps_cost;
    mmaxdistancerepr = -1;
    this->clearESDs();
    if (empty())    return;
    // check for any negative element
    size_t minidx = min_element(begin(), end()) - begin();
    if (at(minidx) < 0.0)
    {
        ostringstream emsg;
        emsg << "DistanceTable::init() negative distance at " << minidx;
        throw InvalidDistanceTable(emsg.str());
    }
    // sort and count unique distance
    sort(begin(), end());
}


void DistanceTable::readESDFormat(istream& fid)
{
    vector<double> positions;
    vector<double> esds;
    string line;
    vector<string> words;
    while (!fid.eof() && getline(fid, line))
    {
        split(line, words);
        if (words.empty() || words[0][0] == '#')    continue;
        istringstream istrs(line);
        double p, e;
        if (istrs >> p >> e)
        {
            positions.push_back(p);
            esds.push_back(e);
        }
        else
        {
            ostringstream emsg;
            int line_start = int(fid.tellg()) - int(line.size());
            emsg << "Invalid ESD format at file position " << line_start;
            throw IOError(emsg.str());
        }
    }
    this->swap(positions);
    this->init();
    this->setESDs(esds);
}


void DistanceTable::readPWAFormat(istream& fid)
{
    vector<double> positions;
    vector<double> widths;
    vector<double> areas;
    string line;
    vector<string> words;
    while (!fid.eof() && getline(fid, line))
    {
        split(line, words);
        if (words.empty() || words[0][0] == '#')    continue;
        istringstream istrs(line);
        double p, w, a;
        if (istrs >> p >> w >> a)
        {
            positions.push_back(p);
            widths.push_back(w);
            areas.push_back(a);
        }
        else
        {
            ostringstream emsg;
            int line_start = int(fid.tellg()) - int(line.size());
            emsg << "Invalid PWA format at file position " << line_start;
            throw IOError(emsg.str());
        }
    }
    // widths and areas are ignored
    this->swap(positions);
    this->init();
}


void DistanceTable::readSimpleFormat(istream& fid)
{
    vector<double> positions;
    double x;
    while (fid >> x)    positions.push_back(x);
    this->swap(positions);
    this->init();
}

// Non-member Operators ------------------------------------------------------

istream& operator>>(istream& fid, DistanceTable& dstbl)
{
    // skip any lines that do not start with number
    istream::pos_type data_start = fid.tellg();
    string line;
    vector<string> esd_marker;
    split("#L position esd", esd_marker);
    vector<string> pwa_marker;
    split("#L position fwhm area", pwa_marker);
    vector<string> words;
    for (; !fid.eof() && getline(fid, line); data_start = fid.tellg())
    {
        split(line, words);
        // check if we found marker of ESD file
        if (words == esd_marker)
        {
            dstbl.readESDFormat(fid);
            break;
        }
        // check if we found marker of PWA file
        if (words == pwa_marker)
        {
            dstbl.readPWAFormat(fid);
            break;
        }
        // otherwise check if line starts with a number
        if (!words.empty())
        {
            // check if the fi
            istringstream wrd0(words[0]);
            double x;
            if (wrd0 >> x && wrd0.eof())
            {
                fid.seekg(data_start);
                dstbl.readSimpleFormat(fid);
                break;
            }
        }
    }
    return fid;
}

// End of file
