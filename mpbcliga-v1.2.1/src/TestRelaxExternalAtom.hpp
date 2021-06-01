/***********************************************************************
* Short Title: unit tests of least squares atom relaxation
*
* Comments:
*
* <license text>
***********************************************************************/

#include <cxxtest/TestSuite.h>

#include "Molecule.hpp"

using namespace std;

class TestRelaxExternalAtom : public CxxTest::TestSuite
{

    private:

        double double_eps;
        DistanceTable dtgt;
        Molecule* mol_tetrahedron;
        Atom_t* vtx_tetrahedron;

    public:

        void setUp()
        {
            double_eps = 1.0e-6;
            vector<double> dst_tetrahedron(6, 1.0);
            dtgt = dst_tetrahedron;
            mol_tetrahedron = new Molecule;
            mol_tetrahedron->setDistanceTable(dtgt);
            mol_tetrahedron->setChemicalFormula(ChemicalFormula("C4"));
            mol_tetrahedron->AddAt("C", -0.5, -sqrt(0.75)*1.0/3, 0.0);
            mol_tetrahedron->AddAt("C", +0.5, -sqrt(0.75)*1.0/3, 0.0);
            mol_tetrahedron->AddAt("C", +0.0, +sqrt(0.75)*2.0/3, 0.0);
            vtx_tetrahedron = new Atom_t("C", 0.0, 0.0, sqrt(2.0/3));
        }


        void tearDown()
        {
            delete mol_tetrahedron; mol_tetrahedron = NULL;
            delete vtx_tetrahedron; vtx_tetrahedron = NULL;
        }


        void test_RelaxTetrahedron()
        {
            Atom_t vtx("C", 1.0, 2.0, 3.0);
            mol_tetrahedron->RelaxExternalAtom(&vtx);
            double dvtx = R3::distance(vtx.r, vtx_tetrahedron->r);
            TS_ASSERT_DELTA(0.0, dvtx, double_eps);
        }

};  // class TestRelaxExternalAtom

// End of file
