/***********************************************************************
* Short Title: unit tests of cost evaluation in Crystal
*
* Comments:
*
* <license text>
***********************************************************************/

#include <cxxtest/TestSuite.h>

#include "Lattice.hpp"
#include "Crystal.hpp"
#include "LigaUtils.hpp"
#include "AtomCost.hpp"
#include "tests_dir.hpp"

using namespace std;

class TestCrystal : public CxxTest::TestSuite
{
    private:

        double double_eps;
        double gradient_eps;
        Crystal crst;
        auto_ptr<Lattice> cubic;
        auto_ptr<Lattice> rhombohedral;
        vector<double> dst_cube;
        vector<double> dst_bcc;
        vector<double> dst_fcc;

    public:

        void setUp()
        {
            double_eps = 100*DOUBLE_EPS;
            gradient_eps = 5e-6;
            crst = Crystal();
            crst.setChemicalFormula("C4");
            crst.setRmax(3.05);
            // do the rest only once
            if (!cubic.get())   cubic.reset(new Lattice(1, 1, 1, 90, 90, 90));
            if (!rhombohedral.get())    rhombohedral.reset(
                    new Lattice(sqrt(0.5), sqrt(0.5), sqrt(0.5), 60, 60, 60));
            // distance data up to first distance over 3
            double cube_data[] = { 1.0, 1.41421356237, 1.73205080757, 2.0,
                2.2360679775, 2.44948974278, 2.82842712475, 3.0, 3.16227766017 };
            double bcc_data[] = { 0.866025403784, 1.0, 1.41421356237,
                1.65831239518, 1.73205080757, 2.0, 2.17944947177, 2.2360679775,
                2.44948974278, 2.59807621135, 2.82842712475, 2.95803989155, 3.0,
                3.16227766017 };
            double fcc_data[] = { 0.707106781187, 1.0, 1.22474487139,
                1.41421356237, 1.58113883008, 1.73205080757, 1.87082869339, 2.0,
                2.12132034356, 2.2360679775, 2.34520787991, 2.44948974278,
                2.5495097568, 2.73861278753, 2.82842712475, 2.91547594742, 3.0,
                3.08220700148 };
            size_t len_cube_data = sizeof(cube_data)/sizeof(double);
            size_t len_bcc_data = sizeof(bcc_data)/sizeof(double);
            size_t len_fcc_data = sizeof(fcc_data)/sizeof(double);
            dst_cube.assign(cube_data, cube_data + len_cube_data);
            dst_bcc.assign(bcc_data, bcc_data + len_bcc_data);
            dst_fcc.assign(fcc_data, fcc_data + len_fcc_data);
        }


        void test_cube()
        {
            crst.setLattice(*cubic);
            crst.setDistanceTable(dst_cube);
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
            crst.AddAt("C", 0.0, 0.0, 0.0);
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
            crst.AddAt("C", 0.1, 0.2, 0.3);
            TS_ASSERT(crst.cost() > 0.0);
            crst.Pop(1);
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
        }


        void test_bcc()
        {
            crst.setLattice(*cubic);
            crst.setDistanceTable(dst_bcc);
            crst.setChemicalFormula("C2");
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
            crst.AddAt("C", 0.0, 0.0, 0.0);
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
            crst.AddAt("C", 0.1, 0.2, 0.3);
            TS_ASSERT(crst.cost() > 0.0);
            // make sure recalculate returns the same non-zero cost
            double cost0 = crst.cost();
            int pcnt0 = crst.countPairs();
            crst.recalculate();
            TS_ASSERT_DELTA(cost0, crst.cost(), double_eps);
            TS_ASSERT_EQUALS(pcnt0, crst.countPairs());
            crst.Pop(1);
            crst.AddAt("C", 0.5, 0.5, 0.5);
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
        }


        void test_fcc()
        {
            crst.setLattice(*cubic);
            crst.setDistanceTable(dst_fcc);
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
            crst.AddAt("C", 0.0, 0.0, 0.0);
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
            crst.AddAt("C", 0.5, 0.5, 0.5);
            TS_ASSERT(crst.cost() > 0.0);
            crst.Pop(1);
            crst.AddAt("C", 0.5, 0.5, 0.0);
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
            crst.AddAt("C", 0.5, 0.0, 0.5);
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
            crst.AddAt("C", 0.0, 0.5, 0.5);
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
            crst.setLattice(*rhombohedral);
            TS_ASSERT(crst.cost() > 0.0);
            crst.setLattice(*cubic);
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
        }


        void test_fcc_rhomb()
        {
            crst.setLattice(*rhombohedral);
            crst.setDistanceTable(dst_fcc);
            crst.setRmax(3.05);
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
            crst.AddAt("C", 0.0, 0.0, 0.0);
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
            crst.AddAt("C", 0.1, 0.2, 0.3);
            TS_ASSERT(crst.cost() > 0.0);
            crst.Pop(1);
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
        }


        R3::Vector analytical_gradient(const Atom_t& a0, AtomCost* atomcost)
        {
            atomcost->eval(a0, AtomCost::GRADIENT);
            R3::Vector grad = atomcost->gradient();
            return grad;
        }


        R3::Vector numerical_gradient(const Atom_t& a0, AtomCost* atomcost)
        {
            const double delta = 1e-8;
            double ac0 = atomcost->eval(a0);
            R3::Vector numgrad;
            for (int i = 0; i < 3; ++i)
            {
                Atom_t a1 = a0;
                a1.r[i] += delta;
                double ac1 = atomcost->eval(a1);
                numgrad[i] = (ac1 - ac0) / delta;
            }
            return numgrad;
        }


        void test_gradient_bcc()
        {
            crst.setLattice(*cubic);
            crst.setDistanceTable(dst_bcc);
            crst.setChemicalFormula("C2");
            crst.AddAt("C", 0.0, 0.0, 0.0);
            // cost gradient
            Atom_t a1("C", 0.51, 0.52, 0.53);
            AtomCost* atomcost = crst.getAtomCostCalculator();
            R3::Vector ga, gn;
            ga = analytical_gradient(a1, atomcost);
            gn = numerical_gradient(a1, atomcost);
            double gdiff = R3::distance(gn, ga);
            TS_ASSERT_DELTA(0.0, gdiff, gradient_eps);
            // overlap gradient
            AtomCost* atomoverlap = crst.getAtomOverlapCalculator();
            TS_ASSERT_EQUALS(0.0, atomoverlap->eval(a1));
            a1.radius = 0.5;
            crst.setAtomRadiiTable("C:0.5");
            TS_ASSERT(atomoverlap->eval(a1) > 0.0)
            ga = analytical_gradient(a1, atomoverlap);
            gn = numerical_gradient(a1, atomoverlap);
            gdiff = R3::distance(gn, ga);
            TS_ASSERT_DELTA(0.0, gdiff, gradient_eps);
        }


        void test_gradient_fcc()
        {
            crst.setLattice(*cubic);
            crst.setDistanceTable(dst_fcc);
            crst.AddAt("C", 0.0, 0.0, 0.0);
            crst.AddAt("C", 0.5, 0.5, 0.0);
            crst.AddAt("C", 0.5, 0.0, 0.5);
            crst.AddAt("C", 0.0, 0.5, 0.5);
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
            Atom_t a1 = crst.getAtom(1);
            crst.Pop(1);
            Atom_t arx = a1;
            R3::Vector offset(0.013, -0.07, -0.03);
            arx.r += offset;
            AtomCost* atomcost = crst.getAtomCostCalculator();
            R3::Vector ga, gn;
            ga = analytical_gradient(arx, atomcost);
            gn = numerical_gradient(arx, atomcost);
            double gdiff = R3::distance(gn, ga);
            TS_ASSERT_DELTA(0.0, gdiff, gradient_eps);
        }


        void test_gradient_triclinic()
        {
            DistanceTable dst;
            ifstream fid(prepend_tests_dir("crystals/triclinic.dst").c_str());
            fid >> dst;
            crst.setDistanceTable(dst);
            crst.ReadFile(prepend_tests_dir("crystals/triclinic.stru"));
            Atom_t a1 = crst.getAtom(1);
            crst.Pop(1);
            // cost gradient
            AtomCost* atomcost = crst.getAtomCostCalculator();
            R3::Vector ga = analytical_gradient(a1, atomcost);
            double ganorm = R3::norm(ga);
            TS_ASSERT_DELTA(0.0, ganorm, gradient_eps);
            Atom_t a2("C", 0.5, 0.5, 0.5);
            R3::Vector gn;
            ga = analytical_gradient(a2, atomcost);
            gn = numerical_gradient(a2, atomcost);
            double gdiff = R3::distance(gn, ga);
            TS_ASSERT_DELTA(0.0, gdiff, gradient_eps);
            // overlap gradient
            AtomCost* atomoverlap = crst.getAtomOverlapCalculator();
            TS_ASSERT_EQUALS(0.0, atomoverlap->eval(a2));
            AtomRadiiTable radii;
            a2.radius = radii["C"] = 0.4;
            crst.setAtomRadiiTable(radii);
            TS_ASSERT(atomoverlap->eval(a2) > 0.0);
            ga = analytical_gradient(a2, atomoverlap);
            gn = numerical_gradient(a2, atomoverlap);
            gdiff = R3::distance(gn, ga);
            TS_ASSERT_DELTA(0.0, gdiff, gradient_eps);
        }


        void test_cube_overlap()
        {
            crst.setChemicalFormula("C1");
            crst.setAtomRadiiTable("C:0.5");
            crst.setLattice(*cubic);
            crst.setDistanceTable(dst_cube);
            crst.AddAt("C", 0.0, 0.0, 0.0);
            TS_ASSERT_EQUALS(1, crst.countAtoms());
            TS_ASSERT_EQUALS(0.0, crst.Overlap());
            crst.setChemicalFormula("C1");
            crst.setAtomRadiiTable("C:0.6");
            crst.Clear();
            TS_ASSERT_EQUALS(0, crst.countAtoms());
            TS_ASSERT_EQUALS(0.0, crst.Overlap());
            crst.AddAt("C", 0.0, 0.0, 0.0);
            TS_ASSERT_EQUALS(1, crst.countAtoms());
            TS_ASSERT_EQUALS(0.6, crst.getAtom(0).radius);
            double c = 6 * 0.2 * 0.2;
            TS_ASSERT_DELTA(c, crst.Overlap(), double_eps);
            crst.recalculate();
            TS_ASSERT_DELTA(c, crst.Overlap(), double_eps);
        }


        void test_bcc_overlap()
        {
            crst.setLattice(*cubic);
            crst.setDistanceTable(dst_bcc);
            crst.setChemicalFormula("AB");
            crst.setAtomRadiiTable("A:0.3, B:0.6");
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
            crst.AddAt("A", 0.0, 0.0, 0.0);
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
            crst.AddAt("B", 0.5, 0.5, 0.5);
            double c0 = (16 * pow(0.9 - sqrt(0.75), 2) + 6 * 0.2 * 0.2) / 2;
            TS_ASSERT_DELTA(c0, crst.cost(), double_eps);
            crst.recalculate();
            TS_ASSERT_DELTA(c0, crst.cost(), double_eps);
            crst.Pop(0);
            double c1 = 6 * 0.2 * 0.2;
            TS_ASSERT_DELTA(c1, crst.cost(), double_eps);
            crst.recalculate();
            TS_ASSERT_DELTA(c1, crst.cost(), double_eps);
            crst.AddAt("A", 0.0, 0.0, 0.0);
            TS_ASSERT_DELTA(c0, crst.cost(), double_eps);
            crst.setAtomRadiiTable("A:0.1, B:0.2");
            TS_ASSERT_DELTA(0.0, crst.cost(), double_eps);
        }


        void test_getNearestAtom()
        {
            crst.setLattice(*cubic);
            crst.setDistanceTable(dst_fcc);
            crst.setChemicalFormula("ABCD");
            AtomPtr ap = crst.getNearestAtom(R3::Vector(1.0, 2.0, 3.0));
            TS_ASSERT(!ap.get());
            crst.AddAt("A", 0.0, 0.0, 0.0);
            crst.AddAt("B", 0.5, 0.5, 0.5);
            crst.AddAt("C", 0.5, 0.5, 0.0);
            crst.AddAt("D", 0.5, 0.0, 0.5);
            R3::Vector rc(1.0, 2.0, 3.0);
            ap = crst.getNearestAtom(rc);
            TS_ASSERT_DELTA(0.0, R3::norm(ap->r - rc), double_eps);
            rc = -4.39, 5.11, -3.66;
            ap = crst.getNearestAtom(rc);
            R3::Vector rcheck(-4.5, 5.0, -3.5);
            TS_ASSERT_DELTA(0.0, R3::norm(rcheck - ap->r), double_eps);
            TS_ASSERT_EQUALS(string("D"), ap->element);
        }


};  // class TestCrystal

// End of file
