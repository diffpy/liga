/*****************************************************************************
* Short Title: run parameters for mpbcliga application
*
* Comments:
*
* <license text>
*****************************************************************************/

#include <cstdlib>
#include <csignal>
#include <fstream>

#include "RunPar_t.hpp"
#include "EmbedPython.hpp"
#include "Random.hpp"
#include "TrialDistributor.hpp"
#include "StringUtils.hpp"
#include "Molecule.hpp"
#include "Crystal.hpp"
#include "Lattice.hpp"
#include "Liga_t.hpp"
#include "AtomFilter_t.hpp"
#include "LigaUtils.hpp"
#include "Exceptions.hpp"
#include "Counter.hpp"
#include "Version.hpp"

using namespace std;
using namespace NS_LIGA;

// Local Helper Functions ----------------------------------------------------

namespace {


void SIGABRT_handler(int signum)
{
    char hostname[255];
    gethostname(hostname, 255);
    cerr << RunPar_t().getAppName() << " on " << hostname <<
        " stopped due to abort signal.  Open in gdb using\n"
        << "kill -CONT " << getpid() << " && gdb --pid=" << getpid() << endl;
    raise(SIGSTOP);
    sleep(1);
}

}   // namespace

//////////////////////////////////////////////////////////////////////////////
// class RunPar_t
//////////////////////////////////////////////////////////////////////////////

// Public Methods ------------------------------------------------------------

void RunPar_t::processArguments(int argc, char* const argv[])
{
    const char* short_options = "p:hV";
    // parameters and options
    option long_options[] = {
        {"parfile", 1, 0, 'p'},
        {"help", 0, 0, 'h'},
        {"version", 0, 0, 'V'},
        {"db-abortstop", 0, 0, '-'},
        {0, 0, 0, 0}
    };
    args.reset(new ParseArgs(argc, argv, short_options, long_options));
    // define aliases for renamed parameters
    args->defParameterAlias("tol_bad", "tolcost");
    args->defParameterAlias("evolve_frac", "promotefrac");
    args->defParameterAlias("evolve_relax", "promoterelax");
    args->defParameterAlias("degenerate_relax", "demoterelax");
    args->defParameterAlias("trials_sharing", "trialsharing");
    args->defParameterAlias("max_dist", "maxbondlength");
    args->Parse();
    if (args->isopt("h") || argc == 1)
    {
        print_help();
        exit(EXIT_SUCCESS);
    }
    else if (args->isopt("V"))
    {
        cout << version_string() << endl;
        exit(EXIT_SUCCESS);
    }
    else if (args->isopt("db-abortstop"))
    {
        signal(SIGABRT, SIGABRT_handler);
    }
    if (args->isopt("p"))
    {
        args->ReadPars(args->opts["p"].c_str());
    }
    args->ValidatePars(validpars());
    // assign run parameters
    this->process_cmdline_args();
    // distfile
    if (!args->ispar("distfile"))
    {
        string emsg = "Distance file not defined.";
        throw ParseArgsError(emsg);
    }
    distfile = args->pars["distfile"];
    DistanceTable dtab;
    ifstream dstfid(distfile.c_str());
    if (!dstfid)
    {
        ostringstream emsg;
        emsg << "Unable to read '" << distfile << "'";
        throw IOError(emsg.str());
    }
    dstfid >> dtab;
    // figure out if we have Molecule or Crystal
    // crystal
    crystal = args->GetPar<bool>("crystal", true);
    // create empty molecule
    if (this->crystal)
    {
        this->mol.reset(new Crystal);
    }
    else
    {
        this->mol.reset(new Molecule);
    }
    // assign distance table
    mol->setDistanceTable(dtab);
    // distreuse
    distreuse = args->GetPar<bool>("distreuse", this->crystal);
    mol->setDistReuse(distreuse);
    // inistru
    if (args->ispar("inistru"))
    {
        this->inistru = args->pars["inistru"];
        this->mol->ReadFile(inistru);
    }
    // latpar
    // must be processed after reading inistru, which may define lattice
    latpar.resize(6);
    if (this->crystal)
    {
        Crystal& crst = dynamic_cast<Crystal&>(*mol);
        const Lattice& L = crst.getLattice();
        latpar[0] = L.a();
        latpar[1] = L.b();
        latpar[2] = L.c();
        latpar[3] = L.alpha();
        latpar[4] = L.beta();
        latpar[5] = L.gamma();
    }
    if (args->ispar("latpar"))
    {
        if (!this->crystal)
        {
            const char* emsg = "latpar has no sense when crystal=false.";
            throw ParseArgsError(emsg);
        }
        latpar = args->GetParVec<double>("latpar");
        if (latpar.size() != 6)
        {
            const char* emsg = "latpar must define 6 lattice parameters.";
            throw ParseArgsError(emsg);
        }
        Lattice L(latpar[0], latpar[1], latpar[2],
                latpar[3], latpar[4], latpar[5]);
        Crystal& crst = dynamic_cast<Crystal&>(*mol);
        crst.setLattice(L);
    }
    // rmax
    if (args->ispar("rmax"))
    {
        if (!this->crystal)
        {
            const char* emsg = "rmax has no sense when crystal=false.";
            throw ParseArgsError(emsg);
        }
        this->rmax = args->GetPar<double>("rmax");
        Crystal& crst = dynamic_cast<Crystal&>(*mol);
        crst.setRmax(rmax);
    }
    if (crystal)
    {
        Crystal& crst = dynamic_cast<Crystal&>(*mol);
        this->rmax = crst.getRmax();
    }
    // outstru
    if (args->ispar("outstru"))
    {
        outstru = args->pars["outstru"];
    }
    // outfmt
    outfmt = args->GetPar<string>("outfmt", "rawxyz");
    Molecule::setOutputFormat(outfmt);
    // saverate, saveall
    if (args->ispar("outstru"))
    {
        // saverate
        saverate = args->GetPar<int>("saverate", 0);
        // saveall
        saveall = args->GetPar<bool>("saveall", false);
    }
    // frames, framestrace, framesrate
    if (args->ispar("frames"))
    {
        // frames
        frames = args->pars["frames"];
        // framestrace
        if (args->ispar("framestrace"))
        {
            framesrate = 0;
            vector<long> stp = args->GetParVec<long>("framestrace");
            if (stp.size() % 3)
            {
                string emsg = "framestrace must have 3n entries.";
                throw ParseArgsError(emsg);
            }
            // check if seasons in stp are ordered
            for (size_t i = 3; i < stp.size(); i += 3)
            {
                if (stp[i] < stp[i-3])
                {
                    string emsg = "framestrace seasons must be sorted.";
                    throw ParseArgsError(emsg);
                }
            }
            TraceId_t tid;
            for (size_t i = 0; i < stp.size(); i += 3)
            {
                tid.season = stp[i];
                tid.level = stp[i+1];
                tid.mol_id = stp[i+2];
                tid.mol_natoms = 0;
                tid.mol_norm_badness = 0.0;
                framestrace.push_back(tid);
            }
        }
        // framesrate - parse only if framestrace is empty
        if (framestrace.size() == 0)
        {
            framesrate = args->GetPar<int>("framesrate", 0);
        }
    }
    // trace
    trace = args->GetPar<bool>("trace", false);
    // scoopfunction
    if (args->ispar("scoopfunction"))
    {
        this->scoopfunction = args->pars["scoopfunction"];
    }
    // scooprate
    this->scooprate = args->ispar("scoopfunction") ?
        args->GetPar("scooprate", 0) : 0;
    // ncpu
    this->ncpu = args->GetPar("ncpu", 1);
    // verbose
    verbose = Liga_t::getDefaultVerbose();
    if (args->ispar("verbose"))
    {
        fill(verbose.begin(), verbose.end(), false);
        vector<string> flagwords = args->GetParVec<string>("verbose");
        vector<string>::iterator w;
        try {
            for (w = flagwords.begin(); w != flagwords.end(); ++w)
            {
                Liga_t::setVerboseVector(verbose, *w, true);
            }
        }
        catch (invalid_argument(e)) {
            throw ParseArgsError(e.what());
        }
    }
    // liga parameters
    // ndim
    ndim = args->GetPar<size_t>("ndim", 3);
    if (ndim < 1 || ndim > 3)
    {
        string emsg = "ndim value must be 1, 2 or 3.";
        throw ParseArgsError(emsg);
    }
    // costweights
    this->costweights.assign(2, 1.0);
    if (args->ispar("costweights"))
    {
        vector<double> cwts = args->GetParVec<double>("costweights");
        this->costweights = cwts;
    }
    if (this->costweights.size() != 2)
    {
        const char* emsg = "costweights must have 2 elements.";
        throw ParseArgsError(emsg);
    }
    mol->setAtomCostScale(this->costweights[0]);
    mol->setAtomOverlapScale(this->costweights[1]);
    // tolcost
    tolcost = args->GetPar<double>("tolcost", 1.0e-4);
    Molecule::tol_nbad = tolcost;
    // natoms - apply only when formula is not set
    if (args->ispar("natoms") && !args->ispar("formula"))
    {
        this->natoms = args->GetPar<int>("natoms");
        this->formula.clear();
        this->formula.push_back(
                ChemicalFormula::value_type("C", this->natoms));
        this->mol->setChemicalFormula(this->formula);
    }
    // formula must be set after distreuse
    if (args->ispar("formula"))
    {
        this->formula.fromString(args->pars["formula"]);
        this->mol->setChemicalFormula(this->formula);
    }
    this->formula = mol->getChemicalFormula();
    // radii
    if (args->ispar("radii"))
    {
        this->radii.fromString(args->pars["radii"]);
        this->mol->setAtomRadiiTable(this->radii);
    }
    // samepairradius
    this->samepairradius = args->GetPar<double>("samepairradius", -1.0);
    this->mol->setSamePairRadius(this->samepairradius);
    // fixed_atoms must be set after inistru
    if (args->ispar("fixed_atoms"))
    {
        fixed_atoms = args->ExpandRangePar("fixed_atoms");
        sort(fixed_atoms.begin(), fixed_atoms.end());
        vector<int>::iterator last;
        last = unique(fixed_atoms.begin(), fixed_atoms.end());
        fixed_atoms.erase(last, fixed_atoms.end());
        vector<int>::const_iterator ii;
        try {
            for (ii = fixed_atoms.begin(); ii != fixed_atoms.end(); ++ii)
            {
                mol->Fix(*ii);
            }
        }
        catch (range_error) {
            ostringstream emsg;
            emsg << "fixed_atoms - invalid index " << *ii << '.';
            throw ParseArgsError(emsg.str());
        }
    }
    base_level = mol->NFixed();
    // maxcputime
    maxcputime = args->GetPar<double>("maxcputime", 0.0);
    // maxwalltime
    maxwalltime = args->GetPar<double>("maxwalltime", 0.0);
    // rngseed
    rngseed = args->GetPar<int>("rngseed", 0);
    if (rngseed)
    {
        randomSeed(rngseed);
    }
    // promotefrac
    promotefrac = args->GetPar<double>("promotefrac", 0.1);
    Molecule::promotefrac = promotefrac;
    // promoterelax
    promoterelax = args->GetPar<bool>("promoterelax", false);
    Molecule::promoterelax = promoterelax;
    // demoterelax
    demoterelax = args->GetPar<bool>("demoterelax", false);
    Molecule::demoterelax = demoterelax;
    // ligasize
    ligasize = args->GetPar<int>("ligasize", 10);
    // stopgame
    stopgame = args->GetPar<double>("stopgame", 0.0025);
    // seasontrials
    seasontrials = int( args->GetPar<double>("seasontrials", 16384.0) );
    // trialsharing
    trialsharing = args->GetPar<string>("trialsharing", "success");
    if (!TrialDistributor::isType(trialsharing))
    {
        ostringstream emsg;
        emsg << "trialsharing must be one of (";
        emsg << join(", ", TrialDistributor::getTypes()) << ").";
        throw ParseArgsError(emsg.str());
    }
    // bangle_range
    if (args->ispar("bangle_range"))
    {
        bangle_range = args->GetParVec<double>("bangle_range");
        if (bangle_range.size() != 2 && bangle_range.size() != 3)
        {
            string emsg = "bangle_range must have 2 or 3 arguments.";
            throw ParseArgsError(emsg);
        }
        double max_blen = bangle_range[0];
        if (bangle_range.size() == 2)   bangle_range.push_back(DOUBLE_MAX);
        BondAngleFilter_t* pbaf = new BondAngleFilter_t(max_blen);
        pbaf->setBondAngleRange(bangle_range[1], bangle_range[2]);
        Molecule::atom_filters.push_back(pbaf);
    }
    // maxbondlength
    if (args->ispar("maxbondlength"))
    {
        maxbondlength = args->GetPar<double>("maxbondlength");
        LoneAtomFilter_t* plaf = new LoneAtomFilter_t(maxbondlength);
        Molecule::atom_filters.push_back(plaf);
    }
    // Final Checks
    this->mol->CheckIntegrity();
    // done
    print_pars();
}


const string& RunPar_t::getAppName() const
{
    static string appname =  "mpbcliga";
    return appname;
}


bool RunPar_t::outOfCPUTime() const
{
    bool rv = (maxcputime > 0.0 && Counter::CPUTime() > maxcputime);
    return rv;
}


bool RunPar_t::outOfWallTime() const
{
    bool rv = (maxwalltime > 0.0 && Counter::WallTime() > maxwalltime);
    return rv;
}


boost::python::object RunPar_t::importMapFunction()
{
    namespace python = boost::python;
    // short circuit when already imported
    if (mmapfunctionobj.get())  return *mmapfunctionobj;
    // performap import
    initializePython();
    mmapfunctionobj.reset(new python::object());
    python::object py_main = python::import("__main__");
    python::object py_globals = py_main.attr("__dict__");
    python::dict py_locals;
    python::object sysmod = python::import("sys");
    int hexversion = python::extract<int>(sysmod.attr("hexversion"));
    if (hexversion >= 0x02060000 && this->ncpu > 1)
    {
        python::object mpmod = python::import("multiprocessing");
        python::object pool = mpmod.attr("Pool")(this->ncpu);
        *mmapfunctionobj = pool.attr("map");
    }
    else
    {
        // print warning on attempt to use several CPUs with Python 2.5
        if (this->ncpu > 1)
        {
            cout << "Warning ncpu=" << this->ncpu <<
                " ignored, it needs Python 2.6 or later.\n";
            this->ncpu = 1;
        }
        *mmapfunctionobj = python::eval("map", py_globals, py_locals);
    }
    return *mmapfunctionobj;
}


boost::python::object RunPar_t::importScoopFunction() const
{
    namespace python = boost::python;
    // short circuit when already imported
    if (mscoopfunctionobj.get())  return *mscoopfunctionobj;
    initializePython();
    mscoopfunctionobj.reset(new python::object());
    // retrieve module and functio name from scoopfunction string
    string::size_type pcolon = this->scoopfunction.find(':');
    string scmodname = this->scoopfunction.substr(0, pcolon);
    string scfncname = this->scoopfunction.substr(pcolon + 1);
    python::object scmod = python::import(scmodname.c_str());
    *mscoopfunctionobj = scmod.attr(scfncname.c_str());
    return *mscoopfunctionobj;
}


double RunPar_t::applyScoopFunction(Molecule* mol) const
{
    namespace python = boost::python;
    double rv;
    try {
        python::object scfnc = this->importScoopFunction();
        python::object stru = mol->convertToDiffPyStructure();
        python::object coststru = scfnc(stru);
        // check if the first value can be converted to a double
        rv = python::extract<double>(coststru[0]);
        python::object scstru  = coststru[1];
        mol->setFromDiffPyStructure(scstru);
    }
    catch (python::error_already_set) {
        if (PyErr_Occurred())   PyErr_Print();
        ostringstream emsg;
        emsg << "Error executing scoopfunction '" <<
            this->scoopfunction << "'\n" <<
            "scoopfunction must return a (cost, stru) tuple.";
        throw ParseArgsError(emsg.str());
    }
    return rv;
}


void RunPar_t::checkScoopFunction(const Molecule& molecule) const
{
    if (this->scoopfunction.empty())  return;
    auto_ptr<Molecule> testmol(molecule.copy());
    double L = testmol->getDistanceTable().maxDistance();
    while (testmol->countAtoms() < testmol->getMaxAtomCount())
    {
        string smbl = testmol->PickElementFromBucket();
        double xc = L * randomFloat();
        double yc = L * randomFloat();
        double zc = L * randomFloat();
        testmol->AddAt(smbl, xc, yc, zc);
    }
    this->applyScoopFunction(testmol.get());
}

// Private Methods -----------------------------------------------------------

void RunPar_t::print_help()
{
    const string& cmd_t = args->cmd_t;
    cout <<
"usage: " << cmd_t << " [-p PARFILE] [DISTFILE] [par1=val1 par2=val2...]\n"
"run " << this->getAppName() <<
    " simulation using distances from DISTFILE.  Parameters can\n"
"be set in PARFILE or on the command line, which overrides PARFILE.\n"
"Options:\n"
"  -p, --parfile=FILE    read parameters from FILE\n"
"  -h, --help            display this message\n"
"  -V, --version         show program version\n"
"  --db-abortstop        stop process on SIGABRT, allows to attach gdb\n"
"IO parameters:\n"
"  distfile=FILE         target distance table\n"
"  inistru=FILE          [empty] initial structure in Cartesian coordinates\n"
"  outstru=FILE          where to save the best full molecule\n"
"  outfmt=string         [rawxyz], cif, discus,... - outstru file format\n"
"  saverate=int          [0] rate of intermediate saves of outstru\n"
"  saveall=bool          [false] save best molecules from all levels\n"
"  frames=FILE           save intermediate structures to FILE.season\n"
"  framesrate=int        [0] number of iterations between frame saves\n"
"  framestrace=array     [] triplets of (season, level, id)\n"
"  trace=bool            [false] keep and show trace of the best structure\n"
"  scoopfunction=string  (python.module:function) top-level structure scooping\n"
"                        scoopfunction must return a (cost, stru) tuple.\n"
"  scooprate=int         [0] rate of top-level structure scooping\n"
"  ncpu=int              [1] number of CPUs used for structure scooping\n"
"  verbose=array         [ad,wc,bc,sc] output flags from\n"
"                        (" << joined_verbose_flags() << ")\n" <<
"Liga parameters:\n"
"  ndim={1,2,3}          [3] search in n-dimensional space\n"
"  crystal=bool          [true] assume periodic crystal structure\n"
"  latpar=array          [1,1,1,90,90,90] crystal lattice parameters\n"
"  rmax=double           [dmax] distance cutoff when crystal=true\n"
"  distreuse=bool        [false] keep used distances in distance table\n"
"  costweights=array     [1,1] weights for distance and overlap components\n"
"  tolcost=double        [1E-4] target normalized molecule cost\n"
"  natoms=int            obsolete, equivalent to formula=Cn\n"
"  formula=string        chemical formula, use inistru when not specified\n"
"  radii=string          define atomic radii in (A1:r1, A2:r2,...) format\n"
"  samepairradius=double [-1] optional radius for a pair of equal atoms\n"
"  fixed_atoms=ranges    [] indices of fixed atoms in inistru (start at 0)\n"
"  maxcputime=double     [0] when set, maximum CPU time in seconds\n"
"  maxwalltime=double    [0] when set, maximum wall time in seconds\n"
"  rngseed=int           seed of random number generator\n"
"  promotefrac=double    [0.1] fraction of tolcost threshold of tested atoms\n"
"  promoterelax=bool     [false] relax the worst atom after addition\n"
"  demoterelax=bool      [false] relax the worst atom after removal\n"
"  ligasize=int          [10] number of teams per division\n"
"  stopgame=double       [0.0025] skip division when winner is worse\n"
"  seasontrials=int      [16384] number of atom placements in one season\n"
"  trialsharing=string   [success] sharing method from (" <<
        join(",", TrialDistributor::getTypes()) << ")\n" <<
"Constrains (applied only when set):\n"
"  bangle_range=array    (max_blen, low[, high]) bond angle constraint\n"
"  maxbondlength=double  distance limit for rejecting lone atoms\n"
;
}


void RunPar_t::process_cmdline_args()
{
    if (args->args.size() > 1)
    {
        ostringstream emsg;
        emsg << args->cmd_t << ": several DISTFILE arguments.";
        throw ParseArgsError(emsg.str());
    }
    if (args->args.size() == 1)
    {
        args->pars["distfile"] = args->args[0];
    }
}


string RunPar_t::version_string(string quote)
{
    using namespace std;
    ostringstream oss;
    oss << quote << this->getAppName() << ' ' <<
        NS_VERSION::getId() << '\n' <<
        quote << "Compiled with " << __VERSION__ << '.';
    return oss.str();
}


void RunPar_t::print_pars()
{
    // print out all run parameters
    // intro messages
    string hashsep(72, '#');
    cout << hashsep << '\n';
    cout << "# " << args->cmd_t << '\n';
    cout << version_string("# ") << '\n';
    char hostname[255];
    gethostname(hostname, 255);
    cout << "# " << hostname << '\n';
    time_t cur_time = time(NULL);
    cout << "# " << ctime(&cur_time);
    cout << hashsep << '\n';
    // distfile
    cout << "distfile=" << distfile << '\n';
    // inistru
    if (args->ispar("inistru"))
    {
        cout << "inistru=" << this->inistru << '\n';
    }
    // outstru
    if (args->ispar("outstru"))
    {
        cout << "outstru=" << outstru << '\n';
    }
    // outfmt
    if (args->ispar("outstru") || args->ispar("frames"))
    {
        cout << "outfmt=" << outfmt << '\n';
    }
    // saverate, saveall
    if (args->ispar("outstru"))
    {
        cout << "saverate=" << saverate << '\n';
        cout << "saveall=" << saveall << '\n';
    }
    // frames, framestrace, framesrate
    if (args->ispar("frames"))
    {
        cout << "frames=" << frames << '\n';
        // framestrace
        if (framestrace.size() != 0)
        {
            cout << "framestrace=";
            for (deque<TraceId_t>::iterator ii = framestrace.begin();
                    ii != framestrace.end(); ++ii)
            {
                cout << " \\" << '\n';
                cout << "    " << ii->season;
                cout << ' ' << ii->level;
                cout << ' ' << ii->mol_id;
            }
            cout << '\n';
        }
        // framesrate - print only if framestrace is empty
        else
        {
            cout << "framesrate=" << framesrate << '\n';
        }
    }
    // trace
    cout << "trace=" << trace << '\n';
    // scoopfunction, scooprate, ncpu
    if (args->ispar("scoopfunction"))
    {
        cout << "scoopfunction=" << this->scoopfunction << '\n';
        cout << "scooprate=" << this->scooprate << '\n';
        cout << "ncpu=" << this->ncpu << '\n';
    }
    // verbose
    {
        vector<string> flagwords;
        // scan flags only up to ALL
        using namespace NS_LIGA_VERBOSE_FLAG;
        for (size_t i = 0; i < ALL; ++i)
        {
            if (verbose[i]) flagwords.push_back(Liga_t::verbose_flags[i]);
        }
        cout << "verbose=" << join(",", flagwords) << '\n';
    }
    // liga parameters
    // ndim
    cout << "ndim=" << ndim << '\n';
    // crystal, latpar, rmax
    cout << "crystal=" << this->crystal << '\n';
    if (this->crystal)
    {
        cout << "latpar=" <<
            this->latpar[0] << ',' << this->latpar[1] << ',' <<
            this->latpar[2] << ',' << this->latpar[3] << ',' <<
            this->latpar[4] << ',' << this->latpar[5] << '\n';
        cout << "rmax=" << this->rmax << '\n';
    }
    // distreuse
    cout << "distreuse=" << distreuse << '\n';
    // costweights
    cout << "costweights=";
    for (size_t i = 0; i != this->costweights.size(); ++i)
    {
        cout << (i ? "," : "") << costweights[i];
    }
    cout << '\n';
    // tolcost
    cout << "tolcost=" << tolcost << '\n';
    // formula
    cout << "formula=" << this->formula.toString() << '\n';
    // radii
    if (!this->radii.empty())
    {
        cout << "radii=" << this->radii.toString() << '\n';
    }
    // samepairradius
    if (samepairradius >= 0.0)
    {
        cout << "samepairradius=" << samepairradius << '\n';
    }
    // fixed_atoms
    if (args->ispar("fixed_atoms") && !fixed_atoms.empty())
    {
        cout << "fixed_atoms=";
        for (int start = 0; start < int(fixed_atoms.size()); )
        {
            int range = 0;
            for ( ; start + range < int(fixed_atoms.size()) &&
                    fixed_atoms[start]+range == fixed_atoms[start+range];
                    ++range )
            { }
            if (start > 0)  cout << ',';
            cout << fixed_atoms[start];
            if (range > 1)  cout << ".." << fixed_atoms[start+range-1];
            start += range;
        }
        cout << '\n';
    }
    // maxcputime
    if (maxcputime > 0.0)
    {
        cout << "maxcputime=" << maxcputime << '\n';
    }
    // maxwalltime
    if (maxwalltime > 0.0)
    {
        cout << "maxwalltime=" << maxwalltime << '\n';
    }
    // rngseed
    if (rngseed)
    {
        cout << "rngseed=" << rngseed << '\n';
    }
    // promotefrac, promoterelax, demoterelax
    cout << "promotefrac=" << promotefrac << '\n';
    cout << "promoterelax=" << promoterelax << '\n';
    cout << "demoterelax=" << demoterelax << '\n';
    // ligasize, stopgame, seasontrials, trialsharing
    cout << "ligasize=" << ligasize << '\n';
    cout << "stopgame=" << stopgame << '\n';
    cout << "seasontrials=" << seasontrials << '\n';
    cout << "trialsharing=" << trialsharing << '\n';
    // constraints
    // bangle_range
    if (args->ispar("bangle_range"))
    {
        cout << "bangle_range=" << bangle_range[0];
        for (size_t i = 1; i < bangle_range.size(); ++i)
        {
            cout << ',' << bangle_range[i];
        }
        cout << '\n';
    }
    // maxbondlength
    if (args->ispar("maxbondlength"))
    {
        cout << "maxbondlength=" << maxbondlength << '\n';
    }
    // finish done
    cout << hashsep << '\n' << endl;
}


const list<string>& RunPar_t::validpars() const
{
    const char* parnames[] = {
        "distfile",
        "inistru",
        "outstru",
        "outfmt",
        "saverate",
        "saveall",
        "frames",
        "framesrate",
        "framestrace",
        "scoopfunction",
        "scooprate",
        "ncpu",
        "verbose",
        "ndim",
        "crystal",
        "latpar",
        "rmax",
        "distreuse",
        "costweights",
        "tolcost",
        "natoms",
        "formula",
        "radii",
        "samepairradius",
        "fixed_atoms",
        "maxcputime",
        "maxwalltime",
        "rngseed",
        "promotefrac",
        "promoterelax",
        "demoterelax",
        "ligasize",
        "stopgame",
        "seasontrials",
        "trace",
        "trialsharing",
        "bangle_range",
        "maxbondlength",
        // obsolete ignored parameters
        "seed_clusters",
        "lookout_prob",
    };
    size_t parlen = sizeof(parnames)/sizeof(char*);
    static list<string> parlist(parnames, parnames + parlen);
    return parlist;
}


const string& RunPar_t::joined_verbose_flags() const
{
    static string jvf;
    if (jvf.empty())
    {
        jvf = join(", ", Liga_t::verbose_flags);
    }
    return jvf;
}

// End of file
