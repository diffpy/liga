#include "tests_dir.hpp"

std::string prepend_tests_dir(const std::string& f)
{
    using namespace std;
    string rv("%(tests_dir)s");
    rv = rv + '/' + f;
    return rv;
}

// vim:ft=cpp:
// End of tests_dir.cpp
