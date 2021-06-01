/***********************************************************************
* Short Title: LIGA statistics keeper
*
* Comments: class for keeping algorithm statistics
*
* <license text>
***********************************************************************/

#ifndef COUNTER_HPP_INCLUDED
#define COUNTER_HPP_INCLUDED

#include <iostream>
#include <string>
#include <map>
#include <ctime>

class Counter
{
    public:

        // friends
        friend std::ostream& operator<<(std::ostream&, const Counter&);

        // types
        typedef unsigned long long ValueType;

        // class methods
        static Counter* getCounter(std::string name);
        static double CPUTime();
        static double WallTime();
        static void printCounters();
        static void printRunStats();

        // public methods
        inline const std::string& name() const  { return _name; }
        inline ValueType value() const          { return _value; }
        inline void count(int tics=1)           { _value += tics; }
        inline void reset(ValueType cnt=0)      { _value = cnt; }

    private:

        // helper class
        class CounterStorage : public std::map<std::string,Counter*>
        {
            public:
                ~CounterStorage();
        };

        // class methods
        static CounterStorage& storage();

        // constructor
        Counter(std::string name);

        // Data Members
        const std::string _name;
        ValueType _value;
        static const time_t _start_walltime;

};

// non-member operators

std::ostream& operator<<(std::ostream& os, const Counter& cnt);

#endif  // COUNTERS_HPP_INCLUDED
