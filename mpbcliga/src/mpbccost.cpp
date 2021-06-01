/*****************************************************************************
* Short Title: calculate cost value of a structure file
*
* Comments:
*
*****************************************************************************/

#include <cstdlib>
#include <memory>
#include <boost/foreach.hpp>

#include "Exceptions.hpp"
#include "RunPar_t.hpp"
#include "Liga_t.hpp"

using namespace std;

const int EXIT_INPUT_ERROR = 2;

namespace {

class RunParCost : public RunPar_t
{
    public:

        // methods
        virtual const string& getAppName() const
        {
            static string appname =  "mpbccost";
            return appname;
        }

        // data
        vector<string> strufiles;

    protected:

        virtual void process_cmdline_args()
        {
            if (this->args->args.size() > 1)
            {
                const string& parfile = this->args->args[0];
                this->args->ReadPars(parfile);
                this->args->ValidatePars(this->validpars());
                this->strufiles = this->args->args;
                this->strufiles.erase(this->strufiles.begin());
            }
            if (this->strufiles.empty())
            {
                const char* emsg = "Structure file not defined.";
                throw ParseArgsError(emsg);
            }
        }


        virtual void print_help()
        {
            const string& cmd_t = args->cmd_t;
            cout <<
                "usage: " << cmd_t << " PARFILE stru1 stru2... [par1=val1 par2=val2...]\n"
                "print cost of structure using PARFILE configuration.  Parameters can\n"
                "be set in PARFILE or on the command line, which overrides PARFILE.\n"
                "Options:\n"
                "  -p, --parfile=FILE    read parameters from FILE\n"
                "  -h, --help            display this message\n"
                "  -V, --version         show program version\n"
                "  --db-abortstop        stop process on SIGABRT, allows to attach gdb\n"
                "IO parameters:\n"
                "  distfile=FILE         target distance table\n"
                "  verbose=array         [] output flags from (all, )\n"
                "Liga parameters:\n"
                "  ndim={1,2,3}          [3] search in n-dimensional space\n"
                "  crystal=bool          [false] assume periodic crystal structure\n"
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
                "  rngseed=int           seed of random number generator\n"
                ;
        }


        virtual void print_pars()
        {
            using NS_LIGA_VERBOSE_FLAG::ALL;
            if (this->verbose[ALL])  this->RunPar_t::print_pars();
        }
};

}   // namespace

//////////////////////////////////////////////////////////////////////////////
// MAIN
//////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    RunParCost rp;
    // Catch exceptions
    try {
        rp.processArguments(argc, argv);
        BOOST_FOREACH (string stru, rp.strufiles)
        {
            rp.mol->ReadFile(stru);
            if (!rp.formula.empty())  rp.mol->setChemicalFormula(rp.formula);
            rp.mol->recalculate();
            const Molecule& m = *(rp.mol);
            cout << stru << ' ' << m.countAtoms() << ' ' << m.cost() <<
                " dc " << m.costDistance() << " oc " << m.costOverlap() <<
                endl;
            using NS_LIGA_VERBOSE_FLAG::ALL;
            if (rp.verbose[ALL])  cout << m;
        }
    }
    catch (IOError(e)) {
        cerr << e.what() << endl;
        return EXIT_INPUT_ERROR;
    }
    catch (ParseArgsError(e)) {
        cerr << e.what() << endl;
        return EXIT_INPUT_ERROR;
    }
    catch (runtime_error(e)) {
        cerr << e.what() << endl;
        return EXIT_INPUT_ERROR;
    }
    catch (invalid_argument(e)) {
        cerr << e.what() << endl;
        return EXIT_INPUT_ERROR;
    }
    return EXIT_SUCCESS;
}

// End of file
