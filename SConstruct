MY_SCONS_HELP = """\
SCons build rules for the mpbcliga program
Usage: scons [target] [var=value]

Targets:

mpbcliga            the Liga structure solver program
mpbccost            program for calculating cost of a given structure
install             install mpbcliga and mpbccost under prefix/bin
alltests            build the unit test program "alltests"
test                execute functional and unit tests for mpbcliga

Build configuration variables:
%s
Variables can be also assigned in a user-written script sconsvars.py.
SCons construction environment can be customized in sconscript.local script.
"""

import os
import platform

def subdictionary(d, keyset):
    return dict(kv for kv in d.items() if kv[0] in keyset)

# copy system environment variables related to compilation
DefaultEnvironment(ENV=subdictionary(os.environ, '''
    PATH PYTHONPATH
    CPATH CPLUS_INCLUDE_PATH LIBRARY_PATH LD_RUN_PATH
    LD_LIBRARY_PATH DYLD_LIBRARY_PATH DYLD_FALLBACK_LIBRARY_PATH
    MACOSX_DEPLOYMENT_TARGET LANG
    '''.split())
)

# Create construction environment
env = DefaultEnvironment().Clone()

# Variables definitions below work only with 0.98.1 or later.
env.EnsureSConsVersion(0, 98, 1)

# Customizable compile variables
vars = Variables('sconsvars.py')

vars.Add('tests',
    'Fixed-string patterns for selecting unit test sources.', None)
vars.Add(EnumVariable('build',
    'compiler settings', 'fast',
    allowed_values=('debug', 'fast')))
vars.Add(EnumVariable('tool',
    'C++ compiler toolkit to be used', 'default',
    allowed_values=('default', 'intelc')))
vars.Add(BoolVariable('profile',
    'build with profiling information', False))
vars.Add(PathVariable('prefix',
    'installation prefix directory', '/usr/local'))
vars.Update(env)
vars.Add(PathVariable('bindir',
    'installation directory for executable [prefix/bin]',
    env['prefix'] + '/bin'))
vars.Update(env)
env.Help(MY_SCONS_HELP % vars.GenerateHelpText(env))

btags = [env['build'], platform.machine()]
if env['profile']:  btags.append('profile')
builddir = env.Dir('build/' + '-'.join(btags))

Export('env')

def GlobSources(pattern):
    """Same as Glob but also require that source node is a valid file.
    """
    rv = [f for f in Glob(pattern) if f.srcnode().isfile()]
    return rv

Export('GlobSources')

if os.path.isfile('sconscript.local'):
    env.SConscript('sconscript.local')

env.SConscript('src/SConscript.main', variant_dir=builddir)

# vim: ft=python
