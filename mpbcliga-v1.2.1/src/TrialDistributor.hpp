/***********************************************************************
* Short Title: declaration of TrialDistributor and related classes
*
* Comments:
*
* <license text>
***********************************************************************/

#ifndef TRIALDISTRIBUTOR_HPP_INCLUDED
#define TRIALDISTRIBUTOR_HPP_INCLUDED

#include <string>
#include <sstream>
#include <stdexcept>
#include <list>
#include <map>
#include <valarray>
#include <deque>

class RunPar_t;

class TrialDistributor
{
    public:

        enum DistributorType { EQUAL, SIZE, SUCCESS };

        // class methods
        static TrialDistributor* create(RunPar_t* rp);
        static TrialDistributor* create(DistributorType tp);
        static std::list<std::string> getTypes();
        static bool isType(const std::string& tp);

        // constructor
        TrialDistributor()
        {
            tolcost = 1.0e-8;
        }
        // destructor
        virtual ~TrialDistributor() { }

        // data - calculated
        std::valarray<double> tshares;

        // methods - class registration and type info
        bool Register();
        virtual DistributorType type() = 0;
        virtual std::string typeStr() = 0;

        // methods
        void setLevelBadness(size_t lv, double bd);
        void setLevelFillRate(size_t lv, double fr);
        void resize(size_t sz);
        inline size_t size()    { return lvbadlog.size(); }
        virtual void share(int seasontrials) = 0;

    protected:

        // types
        typedef std::deque<double> BadnessHistory;

        // class data
        static const size_t histsize;

        // data members
        std::deque<BadnessHistory> lvbadlog;
        std::deque<double> fillrate;
        int base_level;
        int top_level;
        double tolcost;

        // protected methods

    private:

        // registry of derived trial distributors
        static std::map<std::string,DistributorType>& distributorsRegistry();

};  // class TrialDistributor


class TrialDistributorEqual : public TrialDistributor
{
    public:

        // constructor and destructor
        TrialDistributorEqual() : TrialDistributor() { }
        virtual ~TrialDistributorEqual() { }

        // methods
        virtual DistributorType type()  { return EQUAL; }
        virtual std::string typeStr()   { return "equal"; }
        virtual void share(int seasontrials);

};  // class TrialDistributorEqual


class TrialDistributorSize : public TrialDistributor
{
    public:

        // constructor and destructor
        TrialDistributorSize() : TrialDistributor() { }
        virtual ~TrialDistributorSize() { }

        // methods
        virtual DistributorType type()  { return SIZE; }
        virtual std::string typeStr()   { return "size"; }
        virtual void share(int seasontrials);

};  // class TrialDistributorSize


class TrialDistributorSuccess : public TrialDistributor
{
    public:

        // constructor and destructor
        TrialDistributorSuccess() : TrialDistributor() { }
        virtual ~TrialDistributorSuccess() { }

        // methods
        virtual DistributorType type()  { return SUCCESS; }
        virtual std::string typeStr()   { return "success"; }
        virtual void share(int seasontrials);

};  // class TrialDistributorSuccess


#endif  // TRIALDISTRIBUTOR_HPP_INCLUDED
