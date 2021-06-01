/***********************************************************************
* Short Title: Sequence of atoms in structure
*
* Comments:
*
* <license text>
***********************************************************************/

#ifndef ATOMSEQUENCE_INCLUDED
#define ATOMSEQUENCE_INCLUDED

#include <vector>

class Molecule;
class Atom_t;

class AtomSequence
{
    public:

        // constructors
        AtomSequence(const Molecule* pm);
        AtomSequence(std::vector<Atom_t*>& atoms);

        // methods
        inline Atom_t* ptr()    { return *ii; }
        inline Atom_t& ref()    { return **ii; }
        inline void rewind()    { ii = first; }
        inline void next()      { ++ii; }
        inline bool finished()  { return ii == last; }

    private:

        // data
        std::vector<Atom_t*>::iterator ii, first, last;
};

class AtomSequenceIndex : public AtomSequence
{
    public:

        // constructors
        AtomSequenceIndex(const Molecule* pm);
        AtomSequenceIndex(std::vector<Atom_t*>& atoms);

        // methods
        inline int idx()        { return index; }
        inline void next()      { AtomSequence::next(); ++index; }
        inline void rewind()    { AtomSequence::rewind(); index = 0; }

    private:

        // data
        int index;
};

#endif  // ATOMSEQUENCE_INCLUDED
