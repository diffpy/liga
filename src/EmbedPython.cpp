#include <csignal>
#include "EmbedPython.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////////


void initializePython()
{
    static bool is_initialized = false;
    if (is_initialized)
    {
        return;
    }
    static int py_argc = 1;
    static wchar_t arg0[7] = L"python";
    static wchar_t* py_argv[] = {arg0};
    Py_Initialize();
    PySys_SetArgv(py_argc, py_argv);
    // Make sure Python does not eat SIGINT.
    signal(SIGINT, SIG_DFL);
    is_initialized = true;
}


// End of file
