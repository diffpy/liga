/***********************************************************************
*
* Liga Algorithm    for structure determination from pair distances
*                   Pavol Juhas
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* class TestAtomOverlap
*
* Comments: unit tests for AtomOverlap class
*
***********************************************************************/

#include <cxxtest/TestSuite.h>

#include "AtomOverlap.hpp"
#include "DistanceTable.hpp"
#include "Molecule.hpp"

using namespace std;

class TestAtomOverlap : public CxxTest::TestSuite
{
    private:

        DistanceTable dst;
        Molecule mnotouch;
        Molecule moverlap;
        AtomOverlap aoc;

    public:

        TestAtomOverlap() : CxxTest::TestSuite(), aoc(&mnotouch)
        { }


        void setUp()
        {
            double dst_data[1] = { 1.0 };
            dst = DistanceTable(dst_data, 1);
            mnotouch.setDistanceTable(dst);
            mnotouch.setChemicalFormula("C2");
            mnotouch.Clear();
            mnotouch.AddAt("C", -0.5, 0.0, 0.0);
            mnotouch.AddAt("C", +0.5, 0.0, 0.0);
            moverlap = mnotouch;
            moverlap.setAtomRadiiTable("C:2.0");
        }


        void test_notouch()
        {
            const Atom_t& a1 = mnotouch.getAtom(1);
            mnotouch.Pop(1);
            aoc.resetFor(&mnotouch);
            aoc.eval(&a1);
            TS_ASSERT_EQUALS(0.0, aoc.totalCost());
            aoc.eval(&a1, AtomCost::GRADIENT);
            R3::Vector g = aoc.gradient();
            TS_ASSERT_EQUALS(0.0, R3::norm(g));
        }


        void test_overlap()
        {
            Atom_t a1 = moverlap.getAtom(1);
            moverlap.Pop(1);
            // check if radii are properly set
            TS_ASSERT_EQUALS(2.0, a1.radius);
            TS_ASSERT_EQUALS(2.0, moverlap.getAtom(0).radius);
            aoc.resetFor(&moverlap);
            aoc.eval(a1);
            TS_ASSERT_EQUALS(9.0, aoc.totalCost());
            aoc.eval(a1, AtomCost::GRADIENT);
            R3::Vector g;
            g = aoc.gradient();
            TS_ASSERT_EQUALS(-6.0, g[0]);
            TS_ASSERT_EQUALS(0.0, g[1]);
            TS_ASSERT_EQUALS(0.0, g[2]);
            a1.r = -0.5, -1, 0;
            aoc.eval(a1, AtomCost::GRADIENT);
            TS_ASSERT_EQUALS(9.0, aoc.totalCost());
            g = aoc.gradient();
            TS_ASSERT_EQUALS(0.0, g[0]);
            TS_ASSERT_EQUALS(6.0, g[1]);
            TS_ASSERT_EQUALS(0.0, g[2]);
        }

};  // class TestAtomOverlap

// End of file
