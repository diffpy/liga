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
* class TestAtomOverlapCrystal
*
* Comments: unit tests for AtomOverlapCrystal class
*
***********************************************************************/

#include <cxxtest/TestSuite.h>

#include "AtomOverlapCrystal.hpp"
#include "DistanceTable.hpp"
#include "Lattice.hpp"
#include "Crystal.hpp"
#include "LigaUtils.hpp"

using namespace std;

class TestAtomOverlapCrystal : public CxxTest::TestSuite
{
    private:

        double eps;
        DistanceTable dst;
        Lattice cubic;
        Crystal crnotouch;
        Crystal croverlap;
        AtomOverlapCrystal aoc;

    public:

        TestAtomOverlapCrystal() : CxxTest::TestSuite(),
            cubic(1.0, 1.0, 1.0, 90.0, 90.0, 90.0),
            aoc(&crnotouch)
        { }


        void setUp()
        {
            // round-off tolerance
            eps = sqrt(DOUBLE_EPS);
            double dst_data[1] = { 1.0 };
            dst = DistanceTable(dst_data, 1);
            crnotouch.setDistanceTable(dst);
            crnotouch.setLattice(cubic);
            crnotouch.setChemicalFormula("C2");
            crnotouch.Clear();
            crnotouch.AddAt("C", 0.0, 0.0, 0.0);
            crnotouch.AddAt("C", 0.5, 0.5, 0.5);
            croverlap = crnotouch;
            // nasty hack while there is no way to set atom radii
            croverlap.setAtomRadiiTable("C:0.5");
        }


        void test_notouch()
        {
            Atom_t a1 = crnotouch.getAtom(1);
            crnotouch.Pop(1);
            aoc.resetFor(&crnotouch);
            aoc.eval(a1);
            TS_ASSERT_EQUALS(0.0, aoc.totalCost());
            aoc.eval(&a1, AtomCost::GRADIENT);
            R3::Vector g = aoc.gradient();
            TS_ASSERT_EQUALS(0.0, R3::norm(g));
        }


        void test_overlap()
        {
            Atom_t a1 = croverlap.getAtom(1);
            croverlap.Pop(1);
            // check if radii are properly set
            TS_ASSERT_EQUALS(0.5, a1.radius);
            TS_ASSERT_EQUALS(0.5, croverlap.getAtom(0).radius);
            aoc.resetFor(&croverlap);
            aoc.eval(a1);
            double c;
            c = 16 * (1.75 - sqrt(3.0));
            TS_ASSERT_DELTA(c, aoc.totalCost(), eps);
            // gradient is zero due to bcc symmetry
            aoc.eval(a1, AtomCost::GRADIENT);
            R3::Vector g;
            g = aoc.gradient();
            TS_ASSERT_DELTA(0.0, R3::norm(g), eps);
            // move it so it only touches one ball
            a1.r = -0.25, 0.0, 0.0;
            aoc.eval(a1, AtomCost::GRADIENT);
            c = 2 * (0.25 * 0.25 + 0.75 * 0.75);
            TS_ASSERT_EQUALS(c, aoc.totalCost());
            g = aoc.gradient();
            TS_ASSERT_EQUALS(2.0, g[0]);
            TS_ASSERT_EQUALS(0.0, g[1]);
            TS_ASSERT_EQUALS(0.0, g[2]);
            // check sign switching
            a1.r = +0.25, 0.0, 0.0;
            aoc.eval(a1, AtomCost::GRADIENT);
            TS_ASSERT_EQUALS(c, aoc.totalCost());
            g = aoc.gradient();
            TS_ASSERT_EQUALS(-2.0, g[0]);
            TS_ASSERT_EQUALS(0.0, g[1]);
            TS_ASSERT_EQUALS(0.0, g[2]);
        }


};  // class TestAtomOverlapCrystal

// End of file
