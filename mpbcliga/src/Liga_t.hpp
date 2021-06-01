/*****************************************************************************
* Short Title: class Liga_t for competisions organization
*
* Comments:
*
* <license text>
*****************************************************************************/

#ifndef LIGA_T_HPP_INCLUDED
#define LIGA_T_HPP_INCLUDED

#include <memory>
#include <set>
#include <vector>
#include <string>

#include <boost/shared_ptr.hpp>

#include "Division_t.hpp"
#include "RunPar_t.hpp"
#include "Counter.hpp"
#include "LigaUtils.hpp"
#include "TrialDistributor.hpp"

namespace NS_LIGA_VERBOSE_FLAG {

enum VerboseFlag { AD, WC, BC, AV, TS, SC, ALL, VERBOSE_SIZE };

const std::string verbose_flags_array[VERBOSE_SIZE] = {
    "ad", "wc", "bc", "av", "ts", "sc", "all" };

}   // namespace NS_LIGA_VERBOSE_FLAG


class Liga_t : public std::vector<Division_t>
{
    public:

        // types
        typedef NS_LIGA_VERBOSE_FLAG::VerboseFlag VerboseFlag;

        // class data
        static const std::vector<std::string> verbose_flags;

        // class methods
        static std::vector<bool> getDefaultVerbose();
        static void setVerboseVector(std::vector<bool>& vb,
                VerboseFlag flag, bool value=true);
        static void setVerboseVector(std::vector<bool>& vb,
                std::string flag, bool value=true);

        // instance data
        int season;

        // Constructor and destructor
        Liga_t(RunPar_t* runpar);

        // Public methods
        void prepare();
        void playSeason();
        void playLevel(size_t lo_level);
        bool stopFlag() const;
        void useStopFlag(int* flag);
        bool finished() const;
        bool solutionFound() const;
        bool outOfTime() const;
        void printFramesTrace() const;
        void printSummary() const;
        void setVerbose(const std::vector<bool>& vb);
        const std::vector<bool>& getVerbose() const;

    private:

        // Types
        typedef Molecule* PMOL;
        struct epsDoubleCompare : public std::binary_function<double,double,bool>
        {
            bool operator()(const double& x0, const double& x1) const
            {
                return eps_lt(x0, x1);
            }
        };

        // Data members
        RunPar_t* rp;
        int* stopflag;
        int base_level;
        PMOL world_champ;
        std::auto_ptr<Molecule> best_champ;
        mutable bool printed_best_champ;
        mutable bool printed_scooped_structures;
        mutable bool saved_scooped_structures;
        std::auto_ptr<TrialDistributor> tdistributor;
        std::vector<bool> verbose;
        std::auto_ptr<boost::python::list> mscoop_cost_stru;

        // Private methods
        int divSize(int level);
        void shareSeasonTrials();
        PMOL updateWorldChamp();
        void updateBestChamp();
        void printWorldChamp() const;
        void printBestChamp() const;
        void printLevelAverages() const;
        void printTrialShares() const;
        void saveOutStru();
        void saveFrames();
        void recordFramesTrace(std::set<PMOL>& modified, size_t lo_level);
        void saveFramesTrace(std::set<PMOL>& modified, size_t lo_level);
        void prepareScooping();
        void updateScoopedStructures();
        void printScoopedStructures() const;
        void saveScoopedStructures() const;
        void injectBestScoop();
        void injectOverlapMinimization();
        void injectCompetitor(const Molecule*);
};

#endif  // LIGA_T_HPP_INCLUDED
