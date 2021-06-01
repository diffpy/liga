/***********************************************************************
* Short Title: unit tests of const evaluation in Molecule
*
* Comments:
*
* <license text>
***********************************************************************/

#include <valarray>
#include <cxxtest/TestSuite.h>

#include "Molecule.hpp"

using namespace std;

class TestMolecule : public CxxTest::TestSuite
{
    private:

        double double_eps;
        DistanceTable dst_line;
        DistanceTable dst_square;

    public:

        void setUp()
        {
            double_eps = 1.0e-6;
            double line_data[1] = { 1.0 };
            dst_line = DistanceTable(line_data, 1);
            double square_data[6] = { 1.0, 1.0, 1.0, 1.0, sqrt(2.0), sqrt(2.0) };
            dst_square = DistanceTable(square_data, 6);
        }


        void test_line()
        {
            Molecule line;
            line.setDistanceTable(dst_line);
            line.AddAt("", -0.5, 0.0, 0.0);
            line.AddAt("", +0.5, 0.0, 0.0);
            TS_ASSERT_DELTA(0.0, line.cost(), double_eps);
        }


        void test_line_esd()
        {
            Molecule line;
            vector<double> esds;
            esds.push_back(0.5);
            dst_line.setESDs(esds);
            line.setDistanceTable(dst_line);
            line.AddAt("", 0.0, 0.0, 0.0);
            line.AddAt("", 0.9, 0.0, 0.0);
            TS_ASSERT_DELTA(0.01 / 0.25, line.cost(), double_eps);
            Molecule line1 = line;
            TS_ASSERT_DELTA(0.01 / 0.25, line1.cost(), double_eps);
            line1.Pop(1);
            line1.AddAt("", 1.2, 0.0, 0.0);
            TS_ASSERT_DELTA(0.04 / 0.25, line1.cost(), double_eps);
            line1.Pop(1);
            line1.AddAt("", 1.0, 0.0, 0.0);
            TS_ASSERT_DELTA(0.0, line1.cost(), double_eps);
            dst_line.clearESDs();
            Molecule line2;
            line2.setDistanceTable(dst_line);
            line2.Add(line.getAtom(0));
            line2.Add(line.getAtom(1));
            TS_ASSERT_DELTA(0.01, line2.cost(), double_eps);
        }


        void test_square()
        {
            Molecule square;
            square.setDistanceTable(dst_square);
            square.AddAt("", -0.5, -0.5, 0.0);
            square.AddAt("", +0.5, -0.5, 0.0);
            square.AddAt("", +0.5, +0.5, 0.0);
            square.AddAt("", -0.5, +0.5, 0.0);
            TS_ASSERT_DELTA(0.0, square.cost(), double_eps);
        }


        void test_bad_square()
        {
            Molecule bad_square;
            bad_square.setDistanceTable(dst_square);
            bad_square.AddAt("", -0.5, -0.5, 0.0);
            bad_square.AddAt("", +0.5, -0.5, 0.0);
            bad_square.AddAt("", +0.5, +0.5, 0.0);
            bad_square.AddAt("", +0.0, +0.0, 0.0);
            typedef struct {double dij; int i, j;} DIJ;
            TS_ASSERT_EQUALS(6, int(dst_square.size()));
            DIJ sorted_pairs[6] = {
                {sqrt(0.5), 0, 3},
                {sqrt(0.5), 1, 3},
                {sqrt(0.5), 2, 3},
                {1.0, 0, 1},
                {1.0, 1, 2},
                {sqrt(2.0), 0, 2},
            };
            valarray<double> acost(0.0, 4);
            for (size_t didx = 0; didx != 6; ++didx)
            {
                DIJ p = sorted_pairs[didx];
                double dd = dst_square[didx] - p.dij;
                double pcost = dd * dd;
                acost[p.i] += pcost / 2.0;
                acost[p.j] += pcost / 2.0;
            }
            TS_ASSERT_DELTA(acost[0],
                    bad_square.getAtom(0).Badness(), double_eps);
            TS_ASSERT_DELTA(acost[1],
                    bad_square.getAtom(1).Badness(), double_eps);
            TS_ASSERT_DELTA(acost[2],
                    bad_square.getAtom(2).Badness(), double_eps);
            TS_ASSERT_DELTA(acost[3],
                    bad_square.getAtom(3).Badness(), double_eps);
            double totalcost = acost.sum();
            TS_ASSERT_DELTA(totalcost, bad_square.Badness(), double_eps);
            // make sure recalculate returns the same non-zero cost
            bad_square.recalculate();
            TS_ASSERT_DELTA(totalcost, bad_square.Badness(), double_eps);
        }


        void test_getNearestAtom()
        {
            Molecule square;
            AtomPtr ap = square.getNearestAtom(R3::Vector(1.0, 2.0, 3.0));
            TS_ASSERT(!ap.get());
            square.setDistanceTable(dst_square);
            square.AddAt("", -0.5, -0.5, 0.0);
            square.AddAt("", +0.5, -0.5, 0.0);
            square.AddAt("", +0.5, +0.5, 0.0);
            square.AddAt("", -0.5, +0.5, 0.0);
            ap = square.getNearestAtom(R3::Vector(1.0, 2.0, 3.0));
            TS_ASSERT_EQUALS(square.getAtom(2), *ap);
        }

};  // class TestMolecule

// End of file
