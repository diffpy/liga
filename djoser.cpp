/***********************************************************************
* Short Title: molecule reconstruction from distance table
*
* Comments: not particularly inteligent search for molecular configuration
*
* $Id$
***********************************************************************/

#include <limits>
#include <unistd.h>
#include <iomanip>
#include "ioerror.hpp"
#include "ParseArgs.hpp"
#include "BGAlib.hpp"

void print_version()
{
    using namespace std;
    const char *ver =
	"$Id$"
#   if defined(__DATE__) && defined(__TIME__)
	"\ncompiled " __DATE__ " " __TIME__
#   endif
	"\n";
    cout << ver;
}

void print_help(ParseArgs& a)
{
    // /usage:/;/;/-s/.*/"&\\n"/
    // /cou/;/;/s/^\s*"\(.*\)\\n"/\1/ | '[put! ='/*' | /;/put ='*/'
    cout << 
"usage: " << a.cmd_t << "[-p PAR_FILE] [DISTFILE] [par1=val1 par2=val2...]\n"
"run madwalk simulation using distances from DISTFILE.  Parameters can\n"
"be set in PAR_FILE or on the command line, which overrides PAR_FILE.\n"
"Options:\n"
"  -p, --parfile=FILE    read parameters from FILE\n"
"  -h, --help            display this message\n"
"  -v, --version         show program version\n"
"IO parameters:\n"
"  distfile=FILE         target distance table\n"
"  outstru=FILE          where to save the best full molecule\n"
"  inistru=FILE          initial structure [empty box]\n"
"  snapshot=FILE         live molecule structure\n"
"  snaprate=int          [100] number of iterations between snapshot updates\n"
"  frames=FILE           save intermediate structures to FILE.iteration\n"
"  framesrate=int        number of iterations between frame saves\n"
"Walk parameters\n"
"  tol_dd=double         [inf] distance is not used when dd=|d-d0|>=tol_dd\n"
"  tol_bad=double        target value of normalized molecule badness\n"
"  seed=int        	 seed random number generator\n"
"  logsize=int           [10] last steps used for success rate evaluation\n"
"  eprob_max=double      high limit of evolve probability\n"
"  eprob_min=double      low limit of evolve probability\n"
"  bustprob=double       probability of forcing the full structure built\n"
"  evolve_jump=bool      [true] allow additions of several atoms\n"
"  evolve_frac=double    selection badness threshold of tested atoms\n"
"  penalty=string        dd penalty function [pow2], fabs, well\n"
"  dist_trials=int       [10] good distance atoms to try\n"
"  tri_trials=int        [20] godd triangle atoms to try\n"
"  pyr_trials=int        [1000] good pyramid atoms to try\n"
;
}

struct RunPar_t
{
    // IO parameters
    const char *distfile;
    const char *outstru;
    const char *inistru;
    const char *snapshot;
    int snaprate;
    const char *frames;
    int framesrate;
    // Walk parameters
    double tol_dd;
    double tol_bad;
    int seed;
    int logsize;
    double eprob_max;
    double eprob_min;
    double bustprob;
    bool evolve_jump;
    double evolve_frac;
    string penalty;
    int dist_trials;
    int tri_trials;
    int pyr_trials;
};

Molecule process_arguments(RunPar_t& rp, int argc, char *argv[])
{
    char *short_options =
	"p:hv";
    // parameters and options
    option long_options[] = {
	{"parfile", 1, 0, 'p'},
	{"help", 0, 0, 'h'},
	{"version", 0, 0, 'v'},
	{0, 0, 0, 0}
    };
    ParseArgs a(argc, argv, short_options, long_options);
    try {
	a.Parse();
//	a.Dump();
    }
    catch (ParseArgsError) {
	exit(EXIT_FAILURE);
    }
    if (a.opts.count("h") || argc == 1)
    {
	print_help(a);
	exit(EXIT_SUCCESS);
    }
    else if (a.opts.count("v"))
    {
	print_version();
	exit(EXIT_SUCCESS);
    }
    if (a.opts.count("p"))
    {
	try {
	    a.ReadPars(a.opts["p"].c_str());
	}
	catch (IOError(e)) {
	    cerr << e.what() << endl;
	    exit(EXIT_FAILURE);
	}
	catch (ParseArgsError(e)) {
	    cerr << "invalid syntax in parameter file" << endl;
	    cerr << e.what() << endl;
	    exit(EXIT_FAILURE);
	}
    }
    // assign run parameters
    // distfile
    if (a.args.size())
	a.pars["distfile"] = a.args[0];
    if (!a.pars.count("distfile"))
    {
	cerr << "Distance file not defined" << endl;
	exit(EXIT_FAILURE);
    }
    rp.distfile = a.pars["distfile"].c_str();
    DistanceTable* dtab;
    try {
	dtab = new DistanceTable(rp.distfile);
    }
    catch (IOError) {
	exit(EXIT_FAILURE);
    }
    string hashsep(72, '#');
    cout << hashsep << endl;
    cout << "# " << a.cmd_t << ' ' <<
	"$Id$" << endl;
    time_t cur_time = time(NULL);
    cout << "# " << ctime(&cur_time);
    cout << hashsep << endl;
    Molecule mol(*dtab);
    cout << "distfile=" << rp.distfile << endl;
    rp.outstru = a.pars.count("outstru") ? a.pars["outstru"].c_str() : NULL;
    if (rp.outstru) cout << "outstru=" << rp.outstru << endl;
    rp.inistru = a.pars.count("inistru") ? a.pars["inistru"].c_str() : NULL;
    if (rp.inistru) cout << "inistru=" << rp.inistru << endl;
    rp.snapshot = a.pars.count("snapshot") ? a.pars["snapshot"].c_str() : NULL;
    if (rp.snapshot)
    {
	cout << "snapshot=" << rp.snapshot << endl;
	rp.snaprate = a.GetPar<int>("snaprate", 100);
	cout << "snaprate=" << rp.snaprate << endl;
    }
    rp.frames = a.pars.count("frames") ? a.pars["frames"].c_str() : NULL;
    if (rp.frames)
    {
	cout << "frames=" << rp.frames << endl;
	rp.framesrate = a.GetPar<int>("framesrate", 100);
	cout << "framesrate=" << rp.framesrate << endl;
    }
    // Walk parameters
    rp.tol_dd = a.GetPar<double>("tol_dd", 0.1);
    cout << "tol_dd=" << rp.tol_dd << endl;
    mol.tol_dd = rp.tol_dd;
    rp.tol_bad = a.GetPar<double>("tol_bad", 1.0e-4);
    cout << "tol_bad=" << rp.tol_bad << endl;
    rp.seed = a.GetPar<int>("seed", 0);
    if (rp.seed)
    {
	gsl_rng_set(BGA::rng, rp.seed);
	cout << "seed=" << rp.seed << endl;
    }
    rp.logsize = a.GetPar<int>("logsize", 10);
    cout << "logsize=" << rp.logsize << endl;
    rp.eprob_max = a.GetPar<double>("eprob_max", 0.75);
    cout << "eprob_max=" << rp.eprob_max << endl;
    rp.eprob_min = a.GetPar<double>("eprob_min", 0.25);
    cout << "eprob_min=" << rp.eprob_min << endl;
    rp.bustprob = a.GetPar<double>("bustprob", 0.01);
    cout << "bustprob=" << rp.bustprob << endl;
    rp.evolve_jump = a.GetPar<bool>("evolve_jump", true);
    cout << "evolve_jump=" << rp.evolve_jump << endl;
    mol.evolve_jump = rp.evolve_jump;
    rp.evolve_frac = a.GetPar<double>("evolve_frac", 1.0e-4);
    cout << "evolve_frac=" << rp.evolve_frac << endl;
    mol.evolve_frac = rp.evolve_frac;
    rp.penalty = a.GetPar<string>("penalty", "pow2");
    if (rp.penalty == "pow2")
	mol.penalty = BGA::pow2;
    else if (rp.penalty == "well")
	mol.penalty = BGA::well;
    else if (rp.penalty == "fabs")
	mol.penalty = fabs;
    else
    {
	cerr << "Invalid value of penalty parameter" << endl;
	exit(EXIT_FAILURE);
    }
    cout << "penalty=" << rp.penalty << endl;
    rp.dist_trials = a.GetPar("dist_trials", 10);
    cout << "dist_trials=" << rp.dist_trials << endl;
    rp.tri_trials = a.GetPar("tri_trials", 20);
    cout << "tri_trials=" << rp.tri_trials << endl;
    rp.pyr_trials = a.GetPar("pyr_trials", 1000);
    cout << "pyr_trials=" << rp.pyr_trials << endl;
    cout << hashsep << endl << endl;
    return mol;
}

double prob_evolve(const Molecule& mol, RunPar_t& rp, double impr_rate)
{
    double pe;
    if (mol.NAtoms() == mol.max_NAtoms())
    {
	pe = 0.0;
	rp.bustprob = false;
    }
    else if (mol.NAtoms() <= 1)
	pe = 1.0;
    else if (rp.bustprob)
	pe = 1.0;
    else
	pe = impr_rate*(rp.eprob_max-rp.eprob_min)+rp.eprob_min;
    return pe;
}

void evolve_or_degenerate(Molecule& mol, RunPar_t& rp, double pe)
{
    if (pe > gsl_rng_uniform(BGA::rng))
    {
	mol.Evolve(rp.dist_trials, rp.tri_trials, rp.pyr_trials);
	cout << " E " << mol.NAtoms() << " " << mol.NormBadness();
    }
    else
    {
	int Npop = 1;
	if (mol.NormBadness() > rp.tol_bad)
	{
	    Npop = (int) ceil( mol.NAtoms()/4.0 *
		    (1.0 - rp.tol_bad/mol.NormBadness()) );
	    //		would this help?
	    Npop = 1 + gsl_rng_uniform_int(BGA::rng, Npop);
	}
	mol.Degenerate(Npop);
	cout << " D " << mol.NAtoms() << " " << mol.NormBadness();
    }
    cout << endl;
}

void save_snapshot(Molecule& mol, RunPar_t& rp)
{
//  numeric_limits<double> double_info;
    static int cnt = 0;
    static int largest = 0;
    static double bestMNB = numeric_limits<double>().max();
    if (rp.snapshot == NULL || rp.snaprate == 0 || ++cnt < rp.snaprate)
	return;
    if (mol.NAtoms() >= largest && mol.NormBadness() < bestMNB)
    {
	largest = mol.NAtoms();
	bestMNB = mol.NormBadness();
	mol.WriteAtomEye(rp.snapshot);
	cnt = 0;
    }
}

void save_frames(Molecule& mol, RunPar_t& rp, int iteration)
{
//  numeric_limits<double> double_info;
    static int cnt = 0;
    if (rp.frames == NULL || rp.framesrate == 0 || ++cnt < rp.framesrate)
	return;
    ostringstream oss;
    oss << rp.frames << "." << iteration;
    mol.WriteAtomEye(oss.str().c_str());
    cnt = 0;
}

int main(int argc, char *argv[])
{
    RunPar_t rp;
    Molecule mol = process_arguments(rp, argc, argv);
    // set bestMNBadness to a maximum double
    numeric_limits<double> double_info;
    valarray<double> bestMNBadness(double_info.max(), 1+mol.max_NAtoms());
    valarray<int> improved(1, rp.logsize);
    bool bust_now = false;

    for (int iteration = 0; ; ++iteration)
    {
	// calculate probability of evolution
	double impr_rate = 1.0*improved.sum()/rp.logsize;
	if (impr_rate >= 0.5 && rp.bustprob > gsl_rng_uniform(BGA::rng))
	    bust_now = true;
	double pe = prob_evolve(mol, rp, impr_rate);
	cout << iteration;
	evolve_or_degenerate(mol, rp, pe);
//	if (rp.show_abad)
//	    mol.PrintBadness();
	// update bestMNBadness and improved
	int ilog = iteration % rp.logsize;
	if (mol.NormBadness() < bestMNBadness[mol.NAtoms()])
	{
	    bestMNBadness[mol.NAtoms()] = mol.NormBadness();
	    improved[ilog] = 1;
	}
	else
	{
	    improved[ilog] = 0;
	    if (bestMNBadness[mol.NAtoms()] < rp.tol_bad)
		bestMNBadness[mol.NAtoms()] = rp.tol_bad;
	}
	save_snapshot(mol, rp);
	save_frames(mol, rp, iteration);
	if (mol.NAtoms() == mol.max_NAtoms() && mol.NormBadness() < rp.tol_bad)
	{
	    cout << "Solution found!!!" << endl;
	    break;
	}
    }
    // save final structure
    if (rp.outstru)
	mol.WriteAtomEye(rp.outstru);
    return EXIT_SUCCESS;
}
