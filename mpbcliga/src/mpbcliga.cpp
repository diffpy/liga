/*****************************************************************************
* Short Title: molecule reconstruction from distance table
*
* Comments: self-tuning competition search for molecular configuration
*
*****************************************************************************/

#include <cstdlib>
#include <limits>
#include <memory>
#include <unistd.h>
#include <csignal>
#include "ParseArgs.hpp"
#include "Exceptions.hpp"
#include "Liga_t.hpp"

using namespace std;

const int EXIT_INPUT_ERROR = 2;

//////////////////////////////////////////////////////////////////////////////
// SIGHUP handling
//////////////////////////////////////////////////////////////////////////////

int SIGHUP_received = 0;

void SIGHUP_handler(int signum)
{
    // die on 2nd call
    if (SIGHUP_received)    exit(128 + signum);
    SIGHUP_received = signum;
}

//////////////////////////////////////////////////////////////////////////////
// MAIN
//////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    RunPar_t rp;
    auto_ptr<Liga_t> liga;
    // Catch exceptions
    try {
        // process arguments
        rp.processArguments(argc, argv);
        liga.reset(new Liga_t(&rp));
        // watch for HUP
        signal(SIGHUP, SIGHUP_handler);
        liga->useStopFlag(&SIGHUP_received);
        // main loop
        liga->prepare();
        while (!liga->finished())    liga->playSeason();
        liga->printSummary();
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
    // figure out exit code
    int exit_code;
    if (SIGHUP_received)            exit_code = SIGHUP + 128;
    else if (liga->solutionFound()) exit_code = EXIT_SUCCESS;
    else                            exit_code = EXIT_FAILURE;
    // all done here
    return exit_code;
}

// End of file
