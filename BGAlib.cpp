/***********************************************************************
* Short Title: object definitions for Biosphere Genetic Algorithm
*
* Comments:
*
* $Id$
*
* <license text>
***********************************************************************/

#include <stdexcept>
#include <limits>
#include <utility>
#include <map>
#include <gsl/gsl_randist.h>
#include "BGAlib.hpp"

// random number generator
gsl_rng * BGA::rng = gsl_rng_alloc(gsl_rng_default);

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

inline bool read_header(istream& fid)
{
    string dummy;
    return read_header(fid, dummy);
}

// read as many numbers as possible
template<class T> bool read_data(istream& fid, vector<T>& v)
{
    // prepare v
    T x;
    while (fid >> x)
    {
	v.push_back(x);
    }
    return !(fid.rdstate() & ios::badbit);
}

// read as many numbers as possible
template<class T> typename list<T>::iterator list_at(list<T>& lst, int n)
{
    typename list<T>::iterator ii;
    if (n <= lst.size()/2)
    {
	ii = lst.begin();
	advance(ii, n);
    }
    else
    {
	ii = lst.end();
	advance(ii, n-(int)lst.size());
    }
    return ii;
}



////////////////////////////////////////////////////////////////////////
// SandSphere definitions
////////////////////////////////////////////////////////////////////////

SandSphere::SandSphere(int GridMax, const vector<double>& vt) :
    gridmax(GridMax)
{
    init(vt);
}

SandSphere::SandSphere(int GridMax, int s, const double *pt) :
    gridmax(GridMax)
{
    vector<double> vt(s);
    for (int i = 0; i < s; ++i) { vt[i] = pt[i]; }
    init(vt);
}

SandSphere::SandSphere(int GridMax, const char *file) :
    gridmax(GridMax)
{
    // open file for reading
    ifstream fid(file);
    if (!fid)
    {
	cerr << "E: unable to read '" << file << "'" << endl;
	throw IOError();
    }
    // read values to vt
    vector<double> vt;
    bool result = read_header(fid) && read_data(fid, vt);
    // check if everything was read
    if ( !result || !fid.eof() )
    {
	fid.clear();
	cerr << "E: " << file << ':' << fid.tellg() <<
	    ": error reading SandSphere" << endl;
	throw IOError();
    }
    fid.close();
    init(vt);
}

void SandSphere::init(const vector<double>& t)
{
    NDist = t.size();
    if (NDist == 0)
    {
	cerr << "E: target distance table is empty" << endl;
	throw InvalidDistanceTable();
    }
    // calculate and check NAtoms
    double xNAtoms = 0.5 + sqrt(1 + 8.0*NDist)/2.0;
    NAtoms = int(xNAtoms);
    if (double(NAtoms) != xNAtoms)
    {
	cerr << "E: incorrect length of target distance table, NAtoms=" <<
		xNAtoms << '\n';
	throw InvalidDistanceTable();
    }
    // fill in and check distance valarray d
    d.resize(NDist);
    copy(t.begin(), t.end(), &d[0]);
    sort(&d[0], &d[d.size()]);
    if (d[0] <= 0)
    {
	cerr << "E: non-positive entry in target distance table, " <<
	    "d[0]=" << d[0] << '\n';
	throw InvalidDistanceTable();
    }
    // calculate grid parameters
    dmax = d[d.size()-1];
    delta = dmax/gridmax;
    // calculate d2, d2lo, d2hi
    d2.resize(NDist);
    d2lo.resize(NDist);
    d2hi.resize(NDist);
    d /= delta;
    d2 = d*d;
    // take care of grid tolerance:
    SetGridTol(defGridTol);
}

void SandSphere::SetGridTol(double t)
{
    vGridTol = t;
    for (int i = 0; i < NDist; ++i)
    {
	d2lo[i] = (d[i] < vGridTol) ? 0.0 : pow(d[i] - vGridTol, 2);
	d2hi[i] = pow(d[i] + vGridTol, 2);
    }
    for ( list<Molecule *>::const_iterator i = molecules.begin();
	    i != molecules.end(); ++i )
    {
	(*i)->Recalculate();
    }
}

double SandSphere::GridTol()
{
    return vGridTol;
}


////////////////////////////////////////////////////////////////////////
// Atom_t definitions
////////////////////////////////////////////////////////////////////////

Atom_t::Atom_t(double rx0, double ry0, double rz0, double bad0) :
    rx(rx0), ry(ry0), rz(rz0), badness(bad0)
{
    badness_sum = badness;
    age = 1;
}

double Atom_t::Badness() const
{
    return badness;
}

double Atom_t::AvgBadness() const
{
    return (age != 0) ? 1.0*badness_sum/age : 0.0;
}

double Atom_t::IncBadness(double db)
{
    badness += db;
    badness_sum += badness;
    age++;
    return badness;
}

double Atom_t::DecBadness(double db)
{
    badness -= db;
    badness_sum += badness;
    age++;
    return badness;
}

double Atom_t::ResetBadness(double b)
{
    badness = badness_sum = b;
    age = 1;
    return badness;
}

bool operator==(const Atom_t& a1, const Atom_t& a2)
{
    return a1.rx == a2.rx && a1.ry == a2.ry && a1.rz == a2.rz;
}

double dist2(const Atom_t& a1, const Atom_t& a2)
{
    double dr[3] = { (a1.rx - a2.rx), (a1.ry - a2.ry), (a1.rz - a2.rz) };
    return dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
}

////////////////////////////////////////////////////////////////////////
// Pair_t definitions
////////////////////////////////////////////////////////////////////////

bool comp_find_in_array(const int& lhs, const pair<double,double*>& rhs)
{
    const double value = rhs.first;
    const double *a = rhs.second;
    return a[lhs] < value;
}

Pair_t::Pair_t(Molecule *pM, Atom_t& a1, Atom_t& a2) :
    owner(pM), atom1(&a1), atom2(&a2)
{
    d2 = dist2(*atom1, *atom2);
    d = sqrt(d2);
    badness = 0;
    vector<double>::iterator dit = owner->dTarget.find_nearest(d);
    dUsed = *dit;
    badness += pow(d-dUsed, 2);
    owner->dTarget.erase(dit);
    atom1->IncBadness(badness);
    atom2->IncBadness(badness);
    owner->badness += 2*badness;
}

Pair_t::Pair_t(const Pair_t& pair0)
{
    cerr << "call to Pair_t() copy constructor" << endl;
    throw runtime_error("call to Pair_t() copy constructor");
}

Pair_t& Pair_t::operator=(const Pair_t&)
{
    cerr << "call to Pair_t operator=()" << endl;
    throw runtime_error("call to Pair_t operator=()");
}

Pair_t::~Pair_t()
{
    atom1->DecBadness(badness);
    atom2->DecBadness(badness);
    owner->badness -= 2*badness;
    owner->dTarget.return_back(dUsed);
}


////////////////////////////////////////////////////////////////////////
// DistanceTable definitions
////////////////////////////////////////////////////////////////////////

DistanceTable::DistanceTable() : vector<double>()
{
    init();
}

DistanceTable::DistanceTable(const double* v, size_t s) : vector<double>()
{
    resize(s);
    copy(v, v+s, begin());
    init();
}

template <class InputIterator>
DistanceTable::DistanceTable(InputIterator first, InputIterator last) :
    vector<double>(first, last)
{
    init();
}

DistanceTable::DistanceTable(const char* file) : vector<double>()
{
    // open file for reading
    ifstream fid(file);
    if (!fid)
    {
	cerr << "E: unable to read '" << file << "'" << endl;
	throw IOError();
    }
    bool result = read_header(fid) && read_data(fid, *this);
    // check if everything was read
    if ( !result || !fid.eof() )
    {
	fid.clear();
	cerr << "E: " << file << ':' << fid.tellg() <<
	    ": error reading DistanceTable" << endl;
	throw IOError();
    }
    fid.close();
    init();
}

DistanceTable::DistanceTable(const vector<double>& v) : vector<double>(v)
{
    init();
}

DistanceTable& DistanceTable::operator= (const DistanceTable& DT)
{
    *this = vector<double>(DT);
    NAtoms = DT.NAtoms;
    max_d = DT.max_d;
}

DistanceTable& DistanceTable::operator= (const vector<double>& v)
{
    *this = v;
    init();
}

vector<double>::iterator DistanceTable::find_nearest(const double& dfind)
{
    iterator ii = lower_bound(begin(), end(), dfind);
    if (    ( ii == end() && size() != 0 ) ||
	    ( ii != begin() && (dfind - *(ii-1)) < (*ii - dfind) )
       )
	--ii;
    return ii;
}

vector<double>::iterator DistanceTable::return_back(const double& dback)
{
    iterator ii = lower_bound(begin(), end(), dback);
    return insert(ii, dback);
}

void DistanceTable::init()
{
    if (size() == 0)
    {
	cerr << "E: target distance table is empty" << endl;
	throw InvalidDistanceTable();
    }
    // calculate and check NAtoms
    double xNAtoms = 0.5 + sqrt(1 + 8.0*size())/2.0;
    NAtoms = int(xNAtoms);
    if (double(NAtoms) != xNAtoms)
    {
	cerr << "E: incorrect length of target distance table, NAtoms=" <<
		xNAtoms << '\n';
	throw InvalidDistanceTable();
    }
    // sort and check values
    sort(begin(), end());
    if (at(0) <= 0)
    {
	cerr << "E: non-positive entry in DistanceTable, " <<
	    "d[0]=" << at(0) << endl;
	throw InvalidDistanceTable();
    }
    max_d = *rend();
}


////////////////////////////////////////////////////////////////////////
// Molecule definitions
////////////////////////////////////////////////////////////////////////

Molecule::Molecule(const DistanceTable& dtab) : dTarget(dtab)
{
    init();
}

Molecule::Molecule(const DistanceTable& dtab,
	const int s, const double *px, const double *py, const double *pz
	) : dTarget(dtab)
{
    init();
    for (int i = 0; i < s; ++i)
    {
	Add(px[i], py[i], pz[i]);
    }
}

Molecule::Molecule(const DistanceTable& dtab,
	const vector<double>& vx, const vector<double>& vy,
	const vector<double>& vz
	) : dTarget(dtab)
{
    init();
    if (vx.size() != vy.size() || vx.size() != vz.size())
    {
	cerr << "E: invalid coordinate vectors" << endl;
	throw InvalidMolecule();
    }
    for (int i = 0; i < vx.size(); ++i)
    {
	Add(vx[i], vy[i], vz[i]);
    }
}

Molecule::Molecule(const Molecule& M)
{
    init();
    *this  = M;
}

Molecule& Molecule::operator=(const Molecule& M)
{
    if (this == &M) return *this;
    // Clear() must be the first statement
    Clear();
    // eventually this can be speed up by Pair_t copying
    atoms = M.atoms;
    // restore dTarget to the original state
    dTarget = M.dTarget;
    typedef list<Pair_t*>::const_iterator LPPcit;
    for (LPPcit ii = M.pairs.begin(); ii != M.pairs.end(); ++ii)
    {
	dTarget.push_back((*ii)->dUsed);
    }
    sort(dTarget.begin(), dTarget.end());
    // and calculate everything from scratch
    Recalculate();
    // IO helpers
    output_format = M.output_format;
    opened_file = M.opened_file;
    return *this;
}

void Molecule::init()
{
    max_abad = 0;
    badness = 0;
    // default output format
    OutFmtGrid();
}

Molecule::~Molecule()
{
    // we must call Clear() to delete Pair_t objects
    Clear();
}


//////////////////////////////////////////////////////////////////////////
//// Molecule badness/fitness evaluation
//////////////////////////////////////////////////////////////////////////

// recalculate everything from scratch
void Molecule::Recalculate()
{
    if (NAtoms() > max_NAtoms())
    {
	cerr << "E: molecule too large in Recalculate()" << endl;
	throw InvalidMolecule();
    }
    // destroy all pairs
    typedef list<Pair_t*>::iterator LPPit;
    for (LPPit ii = pairs.begin(); ii != pairs.end(); ++ii)
    {
	delete *ii;
    }
    pairs.clear();
    // molecule parameters
    max_abad = -1;
    badness = 0;
    // reset all atoms
    typedef list<Atom_t>::iterator LAit;
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	ai->ResetBadness();
    }
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	LAit aj = ai;
	for (++aj; aj != atoms.end(); ++aj)
	{
	    pairs.push_back(new Pair_t(this, *ai, *aj));
	}
    }
}

//double Molecule::ABadness(int i) const
//{
//    if (!cached) calc_db();
//    return abad[i];
//}
//
//double Molecule::AFitness(int i) const
//{
//    // this will update abadMax if necessary
//    double badness_i = ABadness(i);
//    return abadMax - badness_i;
//}

double Molecule::Badness() const
{
    return badness;
}

double Molecule::Fitness() const
{
    // call to MaxABadness() will update max_abad if necessary
    return NAtoms()*MaxABadness() - Badness();
}

bool comp_Atom_Badness(const Atom_t& lhs, const Atom_t& rhs)
{
    return lhs.Badness() < rhs.Badness();
}

bool comp_Atom_rx(const Atom_t& lhs, const Atom_t& rhs)
{
    return lhs.rx < rhs.rx;
}

bool comp_Atom_ry(const Atom_t& lhs, const Atom_t& rhs)
{
    return lhs.ry < rhs.ry;
}

bool comp_Atom_rz(const Atom_t& lhs, const Atom_t& rhs)
{
    return lhs.rz < rhs.rz;
}

double Molecule::MaxABadness() const
{
    if (max_abad < 0)
    {
	//pj: check this
	max_abad = (NAtoms() == 0) ?  0 : max( (double)NAtoms()-1, 
		max_element(atoms.begin(), atoms.end(),
		    comp_Atom_Badness)-> Badness() );
    }
    return max_abad;
}

//////////////////////////////////////////////////////////////////////////
// Molecule operators
//////////////////////////////////////////////////////////////////////////

Molecule& Molecule::Shift(double dx, double dy, double dz)
{
    for (list<Atom_t>::iterator ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	ai->rx += dx;
	ai->ry += dy;
	ai->rz += dz;
    }
    return *this;
}

Molecule& Molecule::Center()
{
    double avg_rx = 0.0, avg_ry = 0.0, avg_rz = 0.0;
    typedef list<Atom_t>::iterator LAit;
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	avg_rx += ai->rx;
	avg_ry += ai->ry;
	avg_rz += ai->rz;
    }
    avg_rx /= NAtoms();
    avg_ry /= NAtoms();
    avg_rz /= NAtoms();
    Shift(-avg_rx, -avg_ry, -avg_rz);
    return *this;
}

//Molecule& Molecule::Rotate(double phi, double h0, double k0)
//{
//    // define rotation matrix
//    double Rz[2][2] = {
//	{ cos(phi),	-sin(phi)  },
//	{ sin(phi), 	 cos(phi)  }
//    };
//    double hshift, kshift;
//    for (int i = 0; i < NAtoms; ++i)
//    {
//	// shift coordinates origin to the center of rotation
//	double hshift = h[i] - h0;
//	double kshift = k[i] - k0;
//	// rotate
//	double hrot = Rz[0][0]*hshift + Rz[0][1]*kshift;
//	double krot = Rz[1][0]*hshift + Rz[1][1]*kshift;
//	// shift back
//	h[i] = (int) round(hrot + h0);
//	k[i] = (int) round(krot + k0);
//    }
//    return *this;
//}
//
//Molecule& Molecule::Part(const Molecule& M, const int cidx)
//{
//    h.resize(1, M.h[cidx]);
//    k.resize(1, M.k[cidx]);
//    fix_size();
//    return *this;
//}
//
//Molecule& Molecule::Part(const Molecule& M, const list<int>& cidx)
//{
//    vector<int> *h_new, *k_new;
//    if (&M == this)
//    {
//	h_new = new vector<int>;
//	k_new = new vector<int>;
//    }
//    else
//    {
//	h_new = &h;
//	k_new = &k;
//    }
//    h_new->resize(cidx.size());
//    k_new->resize(cidx.size());
//    int idx = 0;
//    for ( list<int>::const_iterator li = cidx.begin();
//	    li != cidx.end(); ++li )
//    {
//	if (*li < 0 || *li >= M.NAtoms)
//	{
//	    throw range_error("in Molecule::Part()");
//	}
//	(*h_new)[idx] = M.h[*li];
//	(*k_new)[idx] = M.k[*li];
//	++idx;
//    }
//    if (&M == this)
//    {
//	h = *h_new;
//	k = *k_new;
//	delete h_new, k_new;
//    }
//    fix_size();
//    return *this;
//}

Molecule& Molecule::Pop(list<Atom_t>::iterator ai)
{
    // delete all pairs that refer to *ai
    Atom_t *ap = &(*ai);
    typedef list<Pair_t*>::iterator LPPit;
    for (LPPit ii = pairs.begin(); ii != pairs.end(); )
    {
	Pair_t* pp = *ii;
	if (pp->atom1 == ap || pp->atom2 == ap)
	{
	    delete *ii;
	    ii = pairs.erase(ii);
	}
	else
	{
	    ++ii;
	}
    }
    atoms.erase(ai);
    // uncache max atom badness
    max_abad = -1;
    return *this;
}

Molecule& Molecule::Pop(const int cidx)
{
    if (cidx < 0 || cidx >= NAtoms())
    {
	throw range_error("in Molecule::Pop(list<int>)");
    }
    list<Atom_t>::iterator apop = list_at(atoms, cidx);
    Pop(apop);
    return *this;
}

Molecule& Molecule::Pop(const list<int>& cidx)
{
    typedef list<Atom_t>::iterator LAit;
    list<LAit> atoms2pop;
    // build a list of iterators to atoms to be popped
    for ( list<int>::const_iterator ii = cidx.begin();
	    ii != cidx.end(); ++ii )
    {
	if (*ii < 0 || *ii >= NAtoms())
	{
	    cerr << "index out of range in Molecule::Pop(list<int>)" << endl;
	    throw range_error("in Molecule::Pop(list<int>)");
	}
	atoms2pop.push_back( list_at(atoms, *ii) );
    }
    // now pop those atoms
    for ( list<LAit>::iterator ii = atoms2pop.begin(); 
	    ii != atoms2pop.end(); ++ii )
    {
	Pop(*ii);
    }
    return *this;
}

Molecule& Molecule::Clear()
{
    // pairs must be destroyed before atoms;
    typedef list<Pair_t*>::iterator LPPit;
    for (LPPit ii = pairs.begin(); ii != pairs.end(); ++ii)
    {
	delete *ii;
    }
    pairs.clear();
    atoms.clear();
    // uncache max atom badness
    max_abad = -1;
    badness = 0;
    return *this;
}

Molecule& Molecule::Add(Molecule& M)
{
    for ( list<Atom_t>::iterator ai = M.atoms.begin();
	    ai != M.atoms.end(); ++ai )
    {
	Add(*ai);
    }
    return *this;
}

Molecule& Molecule::Add(double rx0, double ry0, double rz0)
{
    Add(Atom_t(rx0, ry0, rz0));
    return *this;
}

Molecule& Molecule::Add(Atom_t atom)
{
    if (NAtoms() == max_NAtoms())
    {
	cerr << "E: molecule too large in Add()" << endl;
	throw InvalidMolecule();
    }
    list<Atom_t>::iterator this_atom, ai;
    this_atom = atoms.insert(atoms.end(), atom);
    this_atom->ResetBadness();
    // uncache max atom badness
    max_abad = -1;
    for (ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	if (ai == this_atom)  continue;
	pairs.push_back(new Pair_t(this, *ai, *this_atom));
    }
    return *this;
}

void Molecule::calc_test_badness(Atom_t& ta)
{
    if (NAtoms() == max_NAtoms())
    {
	cerr << "E: Molecule too large in calc_test_badness()" << endl;
	throw InvalidMolecule();
    }
    double tbad = 0;
    typedef vector<double>::iterator VDit;
    list<double> used_distances;
    list<VDit> erased_positions;
    typedef list<Atom_t>::iterator LAit;
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	double td = dist(*ai, ta);
	VDit dnear = dTarget.find_nearest(td);
	tbad += pow(td - *dnear, 2);
	used_distances.push_back(*dnear);
	erased_positions.push_back(dTarget.erase(dnear));
    }
    // restore dTarget to the original state
    list<double>::reverse_iterator rud = used_distances.rbegin();
    list<VDit>::reverse_iterator rpos = erased_positions.rbegin();
    for (; rud != used_distances.rend(); ++rud, ++rpos)
    {
	dTarget.insert(*rpos, *rud);
    }
    ta.ResetBadness(tbad);
}

int Molecule::push_good_distances(
	vector<Atom_t>& vta, double *afit, int ntrials
	)
{
    // prepare discrete generator
    gsl_ran_discrete_t *table = gsl_ran_discrete_preproc(NAtoms(), afit);
    for (int nt = 0; nt < ntrials; ++nt)
    {
	// pick a random direction
	double phi = 2*M_PI*gsl_rng_uniform(BGA::rng);
	double z = 2*gsl_rng_uniform(BGA::rng) - 1.0;
	double w = sqrt(1 - z*z);
	double rdir[3];
	rdir[0] = w*cos(phi);
	rdir[1] = w*sin(phi);
	rdir[2] = z;
	// pick one atom and free distance
	int aidx = gsl_ran_discrete(BGA::rng, table);
	int didx = gsl_rng_uniform_int(BGA::rng, dTarget.size());
	double radius = dTarget[didx];
	Atom_t& a1 = *list_at(atoms, aidx);
	double nrx = a1.rx + rdir[0]*radius;
	double nry = a1.ry + rdir[1]*radius;
	double nrz = a1.rz + rdir[2]*radius;
	Atom_t ad1(nrx, nrz, nrz);
	calc_test_badness(ad1);
	ad1.IncBadness(a1.Badness());
	vta.push_back(ad1);
    }
    gsl_ran_discrete_free(table);
    return ntrials;
}

int Molecule::push_good_triangles(
	vector<Atom_t>& vta, double *afit, int ntrials
	)
{
    // generate randomly oriented triangles
    if (NAtoms() == max_NAtoms())
    {
	cerr << "E: molecule too large for finding a new position" << endl;
	throw InvalidMolecule();
    }
    else if (NAtoms() < 2)
    {
	cerr << "E: molecule too small, triangulation not possible" << endl;
	throw InvalidMolecule();
    }
    gsl_ran_discrete_t *table = gsl_ran_discrete_preproc(NAtoms(), afit);
    int push_count = 0;
    for (int nt = 0; nt < ntrials; ++nt)
    {
	// pick first atom and free distance
	list<int> aidx = random_wt_choose(2, afit, NAtoms());
	list<int>::iterator aidxit = aidx.begin();
	Atom_t& a1 = *list_at(atoms, *(aidxit++)); 
	Atom_t& a2 = *list_at(atoms, *(aidxit++)); 
	int idf1 = gsl_rng_uniform_int(BGA::rng, dTarget.size());
	int idf2 = gsl_rng_uniform_int(BGA::rng, dTarget.size()-1) + 1;
	idf2 = (idf2 + idf1) % dTarget.size();
	double r13 = dTarget[idf1];
	double r23 = dTarget[idf2];
	double r12 = dist(a1, a2);
	// is triangle base reasonably large?
	if (r12 < 1.0)
	    continue;
	double xlong = (r13*r13 + r12*r12 - r23*r23) / (2.0*r12);
	double xperp2 = r13*r13 - xlong*xlong;
	double xperp;
	if (xperp2 > 0.0)
	    xperp = sqrt(xperp2);
	else if (xperp2 > -0.25)
	    xperp = 0.0;
	else
	    continue;
	// direction along triangle base:
	double longdir[3] = {
	    (a2.rx - a1.rx)/r12,
	    (a2.ry - a1.ry)/r12,
	    (a2.rz - a1.rz)/r12 };
	// generate random direction in the plane perpendicular to longdir
	double pdir1[3] = { -longdir[1],  longdir[0],  0.0 };
	if (pdir1[0] == 0 && pdir1[1] == 0)
	    pdir1[0] = 1.0;
	double norm_pdir1 =
	    sqrt(pdir1[0]*pdir1[0] + pdir1[1]*pdir1[1] + pdir1[2]*pdir1[2]);
	pdir1[0] /= norm_pdir1; pdir1[1] /= norm_pdir1; pdir1[2] /= norm_pdir1;
	double pdir2[3] = {
	    longdir[1]*pdir1[2]-longdir[2]*pdir1[1],
	    longdir[2]*pdir1[0]-longdir[0]*pdir1[2],
	    longdir[0]*pdir1[1]-longdir[1]*pdir1[0] };
	double phi = 2*M_PI*gsl_rng_uniform(BGA::rng);
	double perpdir[3] = {
	    cos(phi)*pdir1[0]+sin(phi)*pdir2[0],
	    cos(phi)*pdir1[1]+sin(phi)*pdir2[1],
	    cos(phi)*pdir1[2]+sin(phi)*pdir2[2] };
	// prepare new atom
	double nrx = a1.rx + xlong*longdir[0] + xperp*perpdir[0];
	double nry = a1.ry + xlong*longdir[1] + xperp*perpdir[1];
	double nrz = a1.rz + xlong*longdir[2] + xperp*perpdir[2];
	Atom_t ad2(nrx, nry, nrz);
	// is ad2 already in vta?
	vector<Atom_t>::reverse_iterator rvtadupend = 
	    vta.rbegin()+min((int)vta.size(), 100);
	if (find(vta.rbegin(), rvtadupend, ad2) != rvtadupend)
	    continue;
	calc_test_badness(ad2);
	ad2.IncBadness(a1.Badness() + a2.Badness());
	vta.push_back(ad2);
	push_count++;
    }
    gsl_ran_discrete_free(table);
    return push_count;
}

int Molecule::push_good_pyramids(
	vector<Atom_t>& vta, double *afit, int ntrials
	)
{
    if (NAtoms() == max_NAtoms())
    {
	cerr << "E: molecule too large for finding a new position" << endl;
	throw InvalidMolecule();
    }
    else if (NAtoms() < 3)
    {
	cerr << "E: molecule too small, cannot construct pyramid" << endl;
	throw InvalidMolecule();
    }
    int push_count = 0;
    for (int nt = 0; nt < ntrials;)
    {
	// pick 3 base atoms
	list<int> aidx = random_wt_choose(3, afit, NAtoms());
	list<int>::iterator aidxit = aidx.begin();
	Atom_t& a1 = *list_at(atoms, *(aidxit++)); 
	Atom_t& a2 = *list_at(atoms, *(aidxit++)); 
	Atom_t& a3 = *list_at(atoms, *(aidxit++)); 
	double base_badness = a1.Badness()+a2.Badness()+a3.Badness();
	// pick 3 vertex distances
	list<int> didx = random_choose_few(3, dTarget.size());
	list<int>::iterator didxit = didx.begin();
	double dv1 = dTarget[*(didxit++)];
	double dv2 = dTarget[*(didxit++)];
	double dv3 = dTarget[*(didxit++)];
	double dvperm[6][3] = {
	    {dv1, dv2, dv3},
	    {dv1, dv3, dv2},
	    {dv2, dv1, dv3},
	    {dv2, dv3, dv1},
	    {dv3, dv1, dv2},
	    {dv3, dv2, dv1} };
	for (int iperm = 0; iperm < 6; ++iperm)
	{
	    ++nt;
	    double r14 = dvperm[iperm][0];
	    double r24 = dvperm[iperm][1];
	    double r34 = dvperm[iperm][2];
	    // uvi is a unit vector in a1a2 direction
	    double uvi_val[3] = { a2.rx-a1.rx, a2.ry-a1.ry, a2.rz-a1.rz };
	    valarray<double> uvi(uvi_val, 3);
	    double r12 = vdnorm(uvi);
	    if (r12 < 1.0)
		continue;
	    uvi /= r12;
	    // v13 is a1a3 vector
	    double v13_val[3] = { a3.rx-a1.rx, a3.ry-a1.ry, a3.rz-a1.rz };
	    valarray<double> v13(v13_val, 3);
	    // uvj lies in a1a2a3 plane and is perpendicular to uvi
	    valarray<double> uvj = v13;
	    uvj -= uvi*vddot(uvi, uvj);
	    double nm_uvj = vdnorm(uvj);
	    if (nm_uvj < 1.0)
		continue;
	    uvj /= nm_uvj;
	    // uvk is a unit vector perpendicular to a1a2a3 plane
	    valarray<double> uvk = vdcross(uvi, uvj);
	    double xP1 = -0.5/(r12)*(r12*r12 + r14*r14 - r24*r24);
	    // Pn are coordinates in pyramid coordinate system
	    valarray<double> P1(3);
	    P1[0] = xP1; P1[1] = P1[2] = 0.0;
	    // vT is translation from pyramid to normal system
	    valarray<double> vT(3);
	    vT[0] = a1.rx - xP1*uvi[0];
	    vT[1] = a1.ry - xP1*uvi[1];
	    vT[2] = a1.rz - xP1*uvi[2];
	    // obtain coordinates of P3
	    valarray<double> P3 = P1;
	    P3[0] += vddot(uvi, v13);
	    P3[1] += vddot(uvj, v13);
	    P3[2] = 0.0;
	    double xP3 = P3[0];
	    double yP3 = P3[1];
	    // find pyramid vertices
	    valarray<double> P4(0.0, 3);
	    double h2 = r14*r14 - xP1*xP1;
	    // prepare stop marker for duplicate search
	    vector<Atom_t>::reverse_iterator rvtadupend;
	    // does P4 belong to a1a2 line?
	    if (fabs(h2) < 0.25)
	    {
		// is vertex on a1a2
		if (fabs(vdnorm(P3)) - r14 > 0.25)
		    continue;
		P4 = vT;
		Atom_t ad3(P4[0], P4[1], P4[2]);
		// is ad3 already in vta?
		rvtadupend = vta.rbegin()+min((int)vta.size(), 100);
		if (find(vta.rbegin(), rvtadupend, ad3) != rvtadupend)
		    continue;
		calc_test_badness(ad3);
		ad3.IncBadness(base_badness);
		vta.push_back(ad3);
		push_count++;
		continue;
	    }
	    else if (h2 < 0)
	    {
		continue;
	    }
	    double yP4 = 0.5/(yP3)*(h2 + xP3*xP3 + yP3*yP3 - r34*r34);
	    double z2P4 = h2 - yP4*yP4;
	    // does P4 belong to a1a2a3 plane?
	    if (fabs(z2P4) < 0.25)
	    {
		P4 = yP4*uvj + vT;
		Atom_t ad3(P4[0], P4[1], P4[2]);
		// is ad3 already in vta?
		rvtadupend = vta.rbegin()+min((int)vta.size(), 100);
		if (find(vta.rbegin(), rvtadupend, ad3) != rvtadupend)
		    continue;
		calc_test_badness(ad3);
		ad3.IncBadness(base_badness);
		vta.push_back(ad3);
		push_count++;
		continue;
	    }
	    else if (z2P4 < 0)
	    {
		continue;
	    }
	    // here we can construct 2 pyramids
	    double zP4 = sqrt(z2P4);
	    // top one
	    P4 = yP4*uvj + zP4*uvk + vT;
	    Atom_t ad3top(P4[0], P4[1], P4[2]);
	    // is ad3top already in vta?
	    rvtadupend = vta.rbegin()+min((int)vta.size(), 100);
	    if (find(vta.rbegin(), rvtadupend, ad3top) != rvtadupend)
		continue;
	    calc_test_badness(ad3top);
	    ad3top.IncBadness(base_badness);
	    vta.push_back(ad3top);
	    push_count++;
	    // and bottom one
	    P4 = yP4*uvj - zP4*uvk + vT;
	    Atom_t ad3bottom(P4[0], P4[1], P4[2]);
	    // is ad3bottom already in vta?
	    rvtadupend = vta.rbegin()+min((int)vta.size(), 100);
	    if (find(vta.rbegin(), rvtadupend, ad3bottom) != rvtadupend)
		continue;
	    calc_test_badness(ad3bottom);
	    ad3bottom.IncBadness(base_badness);
	    vta.push_back(ad3bottom);
	    push_count++;
	}
    }
    return push_count;
}


Molecule& Molecule::Evolve(int ntd1, int ntd2, int ntd3)
{
    if (NAtoms() == max_NAtoms())
    {
	cerr << "E: full-sized molecule cannot Evolve()" << endl;
	throw InvalidMolecule();
    }
    vector<Atom_t> vta;
    // evolution is trivial for empty or 1-atom molecule
    switch (NAtoms())
    {
	case 0:
	    Add(0, 0, 0);
	    return *this;
	case 1:
	    double afit1 = 1.0;
	    push_good_distances(vta, &afit1, 1);
	    Add(vta[0]);
	    Center();
	    return *this;
    }
    // otherwise we need to build array of atom fitnesses
    double afit[NAtoms()];
    double *pf = afit;
    typedef list<Atom_t>::iterator LAit;
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai, ++pf)
    {
	*pf = MaxABadness() - ai->Badness();
    }
    // and push appropriate numbers of test atoms
    // here NAtoms() >= 2
    push_good_distances(vta, afit, ntd1);
    push_good_triangles(vta, afit, ntd2);
    if (NAtoms() > 2)  push_good_pyramids(vta, afit, ntd3);
    // select among the minimal elements of vta
    double min_badness = min_element(vta.begin(), vta.end(),
	    comp_Atom_Badness)-> Badness();
    typedef vector<Atom_t>::iterator VAit;
    vector<VAit> min_iterators;
    for (VAit ai = vta.begin(); ai != vta.end(); ++ai)
    {
	if (ai->Badness() == min_badness)
	    min_iterators.push_back(ai);
    }
    int itidx = gsl_rng_uniform_int(BGA::rng, min_iterators.size());
    Atom_t& best = *min_iterators[itidx];
    Add(best);
    if (NAtoms() < 40)    Center();
    return *this;
}

Molecule& Molecule::Degenerate(int Npop)
{
    Npop = min(NAtoms(), Npop);
    if (Npop == 0)  return *this;
    // build array of atom badnesses
    double abad[NAtoms()];
    double *pb = abad;
    typedef list<Atom_t>::iterator LAit;
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai, ++pb)
    {
	*pb = ai->Badness();
    }
    // generate list of atoms to pop
    list<int> ipop = random_wt_choose(Npop, abad, NAtoms());
    Pop(ipop);
    if (NAtoms() < 40)    Center();
    return *this;
}

list<int> random_choose_few(int K, int Np)
{
    list<int> lst;
    if (K > Np)
    {
	cerr << "random_wt_choose(): too many items to choose" << endl;
	throw range_error("too many items to choose");
    }
    // check trivial case
    else if (K == 0)
    {
	return lst;
    }
    int N = Np;
    map<int,int> tr;
    for (int i = 0; i < K; ++i, --N)
    {
	int k = gsl_rng_uniform_int(BGA::rng, N);
	for (int c = 0;  tr.count(k) == 1;  k = tr[k], ++c)
	{
	    if (!(c < Np))
	    {
		cerr << "random_choose_few(): too many translations" << endl;
		throw runtime_error("too many translations");
	    }
	}
	lst.push_back(k);
	tr[k] = N-1;
    }
    return lst;
}

list<int> random_wt_choose(int K, const double *p, int Np)
{
    list<int> lst;
    if (K > Np)
    {
	cerr << "random_wt_choose(): too many items to choose" << endl;
	throw range_error("too many items to choose");
    }
    // check trivial case
    else if (K == 0)
    {
	return lst;
    }
    if ( p+Np != find_if(p, p+Np, bind2nd(less<double>(),0.0)) )
    {
	cerr << "random_wt_choose(): negative choice probability" << endl;
	throw runtime_error("negative choice probability");
    }
    // now we need to do some real work
    double prob[Np];
    copy(p, p+Np, prob);
    // integer encoding
    int val[Np];
    for (int i = 0; i != Np; ++i)  val[i] = i;
    // cumulative probability
    double cumprob[Np];
    int Nprob = Np;
    // main loop
    for (int i = 0, Nprob = Np; i != K; ++i, --Nprob)
    {
	// calculate cumulative probability
	partial_sum(prob, prob+Nprob, cumprob);
	// if all probabilities are 0.0, set them to equal value
	if (cumprob[Nprob-1] == 0.0)
	{
	    for (int j = 0; j != Nprob; ++j)
	    {
		prob[j] = 1.0;
		cumprob[j] = (j+1.0)/Nprob;
	    }
	}
	// otherwise we can normalize cumprob
	else
	{
	    for (double *pcp = cumprob; pcp != cumprob+Nprob; ++pcp)
	    {
		*pcp /= cumprob[Nprob-1];
	    }
	}
	// now let's do binary search on cumprob:
	double r = gsl_rng_uniform(BGA::rng);
	double *pcp = upper_bound(cumprob, cumprob+Nprob, r);
	int idx = pcp - cumprob;
	lst.push_back(val[idx]);
	// overwrite this element with the last number
	prob[idx] = prob[Nprob-1];
	val[idx] = val[Nprob-1];
    }
    return lst;
}

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
	cerr << "vdcross(): invalid valarray size" << endl;
	throw runtime_error("vdcross(): invalid valarray size");
    }
    valarray<double> cross(3);
    cross[0] = v1[1]*v2[2] - v1[2]*v2[1];
    cross[1] = v1[2]*v2[0] - v1[0]*v2[2];
    cross[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return cross;
}

////////////////////////////////////////////////////////////////////////
// Molecule IO functions
////////////////////////////////////////////////////////////////////////

Molecule::ParseHeader::ParseHeader(const string& s) : header(s)
{
    // parse format
    string fmt;
    // initialize members:
    state =
	read_token("BGA molecule format", fmt) &&
	read_token("NAtoms", NAtoms);
    if (!state)
    {
	return;
    }
    if (fmt == "grid")
	format = GRID;
    else if (fmt == "xyz")
	format = XYZ;
    else if (fmt == "atomeye")
	format = ATOMEYE;
    else
    {
	state = false;
	return;
    }
}

template<class T> bool Molecule::ParseHeader::read_token(
	const char *token, T& value
	)
{
    const char *fieldsep = ":= ";
    int ltoken = strlen(token);
    string::size_type sp;
    const string::size_type npos = string::npos;
    if (
	    npos == (sp = header.find(token)) ||
	    npos == (sp = header.find_first_not_of(fieldsep, sp+ltoken))
       )
    {
	return false;
    }
    istringstream istrs(header.substr(sp));
    bool result = (istrs >> value);
    return result;
}

istream& Molecule::ReadGrid(istream& fid)
{
    // read values to integer vector vhkl
    string header;
    vector<int> vhkl;
    bool result = read_header(fid, header) && read_data(fid, vhkl);
    if (!result) return fid;
    // parse header
    int vhkl_NAtoms = vhkl.size()/3;
    ParseHeader ph(header);
    if (ph)
    {
	if ( vhkl_NAtoms != ph.NAtoms )
	{
	    cerr << "E: " << opened_file << ": expected " << ph.NAtoms <<
		" atoms, read " << vhkl_NAtoms << endl;
	    throw IOError();
	}
    }
    // check if all coordinates have been read
    if ( vhkl.size() % 3 )
    {
	cerr << "E: " << opened_file << ": incomplete data" << endl;
	throw IOError();
    }
    Clear();
    for (int i = 0; i < vhkl.size(); i += 3)
    {
	int h = vhkl[i+0];
	int k = vhkl[i+1];
	int l = vhkl[i+2];
	Add(Atom_t(h, k, l));
    }
    return fid;
}

bool Molecule::ReadGrid(const char* file)
{
    // open file for reading
    ifstream fid(file);
    if (!fid)
    {
	cerr << "E: unable to read '" << file << "'" << endl;
	throw IOError();
    }
    opened_file = file;
    bool result = ReadGrid(fid);
    opened_file.clear();
    fid.close();
    return result;
}

istream& Molecule::ReadXYZ(istream& fid)
{
    // read values to integer vector vxyz
    string header;
    vector<double> vxyz;
    bool result = read_header(fid, header) && read_data(fid, vxyz);
    if (!result) return fid;
    int vxyz_NAtoms = vxyz.size()/3;
    // check how many numbers were read
    ParseHeader ph(header);
    if (ph && vxyz_NAtoms != ph.NAtoms)
    {
	cerr << "E: " << opened_file << ": expected " << ph.NAtoms <<
	    " atoms, read " << vxyz_NAtoms << endl;
	throw IOError();
    }
    if ( vxyz.size() % 3 )
    {
	cerr << "E: " << opened_file << ": incomplete data" << endl;
	throw IOError();
    }
    Clear();
    for (int i = 0; i < vxyz.size(); i += 3)
    {
	Add(Atom_t(vxyz[i+0], vxyz[i+1], vxyz[i+2]));
    }
    return fid;
}

bool Molecule::ReadXYZ(const char* file)
{
    // open file for reading
    ifstream fid(file);
    if (!fid)
    {
	cerr << "E: unable to read '" << file << "'" << endl;
	throw IOError();
    }
    opened_file = file;
    bool result = ReadXYZ(fid);
    opened_file.clear();
    fid.close();
    return result;
}

bool write_file(const char* file, Molecule& M)
{
    // open file for writing
    ofstream fid(file);
    if (!fid)
    {
	cerr << "E: unable to write to '" << file << "'" << endl;
	throw IOError();
    }
    bool result = (fid << M);
    fid.close();
    return result;
}

bool Molecule::WriteGrid(const char* file)
{
    file_fmt_type org_ofmt = output_format;
    OutFmtGrid();
    bool result = write_file(file, *this);
    output_format = org_ofmt;
    return result;
}

bool Molecule::WriteXYZ(const char* file)
{
    file_fmt_type org_ofmt = output_format;
    OutFmtXYZ();
    bool result = write_file(file, *this);
    output_format = org_ofmt;
    return result;
}

bool Molecule::WriteAtomEye(const char* file)
{
    file_fmt_type org_ofmt = output_format;
    OutFmtAtomEye();
    bool result = write_file(file, *this);
    output_format = org_ofmt;
    return result;
}

Molecule& Molecule::OutFmtGrid()
{
    output_format = GRID;
    return *this;
}

Molecule& Molecule::OutFmtXYZ()
{
    output_format = XYZ;
    return *this;
}

Molecule& Molecule::OutFmtAtomEye()
{
    output_format = ATOMEYE;
    return *this;
}

istream& operator>>(istream& fid, Molecule& M)
{
    string header;
    istream::pos_type p = fid.tellg();
    if( !read_header(fid, header) )
    {
	fid.setstate(ios_base::failbit);
	return fid;
    }
    fid.seekg(p);
    Molecule::ParseHeader ph(header);
    if (!ph)
    {
	fid.setstate(ios_base::failbit);
	return fid;
    }
    bool result;
    switch (ph.format)
    {
	case M.GRID:
	    result = M.ReadGrid(fid);
	    break;
	case M.XYZ:
	    result = M.ReadXYZ(fid);
	    break;
	case M.ATOMEYE:
	    throw runtime_error("reading of atomeye files not implemented");
	    break;
    }
    if (!result)
    {
	fid.setstate(ios_base::failbit);
    }
    return fid;
}

ostream& operator<<(ostream& fid, Molecule& M)
{
    typedef list<Atom_t>::iterator LAit;
		LAit afirst = M.atoms.begin();
		LAit alast = M.atoms.end();
    switch (M.output_format)
    {
	case M.GRID:
	    fid << "# BGA molecule format = grid" << endl;
	    fid << "# NAtoms = " << M.NAtoms() << endl;
	    for (LAit ai = afirst; ai != alast; ++ai)
	    {
		fid << ai->rx << '\t' << ai->ry << '\t' << ai->rz << endl;
	    }
	    break;
	case M.XYZ:
	    fid << "# BGA molecule format = xy" << endl;
	    fid << "# NAtoms = " << M.NAtoms() << endl;
	    for (LAit ai = afirst; ai != alast; ++ai)
	    {
		fid <<
		    ai->rx << '\t' <<
		    ai->ry << '\t' <<
		    ai->rz << endl;
	    }
	    break;
	case M.ATOMEYE:
	    double xyz_lo = 0.0;
	    double xyz_hi = 1.0;
	    double xyz_range = xyz_hi - xyz_lo;
	    if (M.NAtoms() > 0)
	    {
		const double scale = 1.01;
		double xyz_extremes[8] = {
		    -M.dTarget.max_d,
		    scale * min_element(afirst, alast, comp_Atom_rx)->rx,
		    scale * min_element(afirst, alast, comp_Atom_ry)->ry,
		    scale * min_element(afirst, alast, comp_Atom_rz)->rz,
		    +M.dTarget.max_d,
		    scale * max_element(afirst, alast, comp_Atom_rx)->rx,
		    scale * max_element(afirst, alast, comp_Atom_ry)->ry,
		    scale * max_element(afirst, alast, comp_Atom_rz)->rz,
		};
		xyz_lo = *min_element(xyz_extremes, xyz_extremes+8);
		xyz_hi = *max_element(xyz_extremes, xyz_extremes+8);
		xyz_range = xyz_hi - xyz_lo;
	    }
	    fid << "# BGA molecule format = atomeye" << endl;
	    fid << "# NAtoms = " << M.NAtoms() << endl;
	    fid << "Number of particles = " << M.NAtoms() << endl;
	    fid << "A = 1.0 Angstrom (basic length-scale)" << endl;
	    fid << "H0(1,1) = " << xyz_range << " A" << endl;
	    fid << "H0(1,2) = 0 A" << endl;
	    fid << "H0(1,3) = 0 A" << endl;
	    fid << "H0(2,1) = 0 A" << endl;
	    fid << "H0(2,2) = " << xyz_range << " A" << endl;
	    fid << "H0(2,3) = 0 A" << endl;
	    fid << "H0(3,1) = 0 A" << endl;
	    fid << "H0(3,2) = 0 A" << endl;
	    fid << "H0(3,3) = " << xyz_range << " A" << endl;
	    fid << ".NO_VELOCITY." << endl;
	    // 4 entries: x, y, z, Uiso
	    fid << "entry_count = 4" << endl;
	    fid << "auxiliary[0] = abad [au]" << endl;
	    fid << endl;
	    // pj: now it only works for a single Carbon atom in the molecule
	    fid << "12.0111" << endl;
	    fid << "C" << endl;
	    for (LAit ai = afirst; ai != alast; ++ai)
	    {
		fid <<
		    (ai->rx - xyz_lo) / xyz_range << " " <<
		    (ai->ry - xyz_lo) / xyz_range << " " <<
		    (ai->rz - xyz_lo) / xyz_range << " " <<
		    ai->Badness() << endl;
	    }
	    break;
    }
    return fid;
}

void Molecule::PrintBadness()
{
    // call to MBadness() will update abadMax if necessary
    cout << "MBadness() = " << Badness() << endl;
    cout << "ABadness() =";
    double mab = max_element(atoms.begin(), atoms.end(),
		    comp_Atom_Badness) -> Badness();
    bool marked = false;
    typedef list<Atom_t>::iterator LAit;
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	cout << ' ' << ai->Badness();
	if (!marked && ai->Badness() == mab)
	{
	    cout << '*';
	    marked = true;
	}
    }
    cout << endl;
}

void Molecule::PrintFitness()
{
    cout << "MFitness() = " << Fitness() << endl;
    cout << "AFitness() =";
    double mab = max_element(atoms.begin(), atoms.end(),
		    comp_Atom_Badness) -> Badness();
    bool marked = false;
    typedef list<Atom_t>::iterator LAit;
    for (LAit ai = atoms.begin(); ai != atoms.end(); ++ai)
    {
	cout << ' ' << MaxABadness() - ai->Badness();
	if (!marked && ai->Badness() == mab)
	{
	    cout << '*';
	    marked = true;
	}
    }
    cout << endl;
}
