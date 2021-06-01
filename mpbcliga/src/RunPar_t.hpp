/*****************************************************************************
* Short Title: run parameters for mpbcliga application
*
* Comments:
*
* <license text>
*****************************************************************************/

#ifndef RUNPAR_T_HPP_INCLUDED
#define RUNPAR_T_HPP_INCLUDED

#include <deque>
#include <string>
#include <memory>

#include <boost/python.hpp>

#include "ParseArgs.hpp"
#include "TraceId_t.hpp"
#include "Molecule.hpp"

class RunPar_t
{
    public:

        // constructor
        virtual ~RunPar_t()  { }
        // methods
        void processArguments(int argc, char* const argv[]);
        virtual const std::string& getAppName() const;
        bool outOfCPUTime() const;
        bool outOfWallTime() const;
        boost::python::object importMapFunction();
        boost::python::object importScoopFunction() const;
        double applyScoopFunction(Molecule* mol) const;
        void checkScoopFunction(const Molecule&) const;
        // parsed input arguments
        std::auto_ptr<ParseArgs> args;
        // Output option
        bool trace;
        // IO parameters
        std::string distfile;
        std::string inistru;
        std::string outstru;
        std::string outfmt;
        int saverate;
        bool saveall;
        std::string frames;
        int framesrate;
        std::deque<TraceId_t> framestrace;
        std::string scoopfunction;
        int scooprate;
        mutable int ncpu;
        std::vector<bool> verbose;
        // Liga parameters
        size_t ndim;
        bool crystal;
        std::vector<double> latpar;
        double rmax;
        bool distreuse;
        double tolcost;
        std::vector<double> costweights;
        int natoms;
        ChemicalFormula formula;
        AtomRadiiTable radii;
        double samepairradius;
        std::vector<int> fixed_atoms;
        double maxcputime;
        double maxwalltime;
        int rngseed;
        double promotefrac;
        bool promoterelax;
        bool demoterelax;
        int ligasize;
        double stopgame;
        int seasontrials;
        std::string trialsharing;
        // generated data
        std::auto_ptr<Molecule> mol;
        int base_level;
        // Constrains
        std::vector<double> bangle_range;
        double maxbondlength;

    protected:

        virtual void process_cmdline_args();
        std::string version_string(std::string quote="");
        const std::list<std::string>& validpars() const;
        const std::string& joined_verbose_flags() const;
        virtual void print_help();
        virtual void print_pars();

    private:

        mutable std::auto_ptr<boost::python::object> mmapfunctionobj;
        mutable std::auto_ptr<boost::python::object> mscoopfunctionobj;

};

#endif  // RUNPAR_T_HPP_INCLUDED
