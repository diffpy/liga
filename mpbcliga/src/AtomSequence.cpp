/***********************************************************************
* Short Title: Sequence of atoms in structure
*
* Comments:
*
* <license text>
***********************************************************************/

#include "AtomSequence.hpp"
#include "Molecule.hpp"

using namespace std;


////////////////////////////////////////////////////////////////////////
// class AtomSequence
////////////////////////////////////////////////////////////////////////

// constructors

AtomSequence::AtomSequence(const Molecule* pm)
{
    Molecule* mol = const_cast<Molecule*>(pm);
    first = mol->atoms.begin();
    last = mol->atoms.end();
    rewind();
}

AtomSequence::AtomSequence(std::vector<Atom_t*>& atoms)
{
    first = atoms.begin();
    last = atoms.end();
    rewind();
}


////////////////////////////////////////////////////////////////////////
// class AtomSequenceIndex
////////////////////////////////////////////////////////////////////////

AtomSequenceIndex::AtomSequenceIndex(const Molecule* pm) :
    AtomSequence(pm), index(0)
{ }

AtomSequenceIndex::AtomSequenceIndex(vector<Atom_t*>& atoms) :
    AtomSequence(atoms), index(0)
{ }
