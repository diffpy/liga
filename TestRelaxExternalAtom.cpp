/***********************************************************************
* Short Title: unit tests of least squares atom relaxation
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#include <stdexcept>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "BGAlib.hpp"

using namespace std;

class TestRelaxExternalAtom : public CppUnit::TestFixture
{

    CPPUNIT_TEST_SUITE(TestRelaxExternalAtom);
    CPPUNIT_TEST(test_RelaxTetrahedron);
    CPPUNIT_TEST_SUITE_END();

private:

    double double_eps;
    DistanceTable dTarget;
    Molecule* mol_tetrahedron;
    Atom_t* vtx_tetrahedron;

public:

    void setUp()
    {
	double_eps = 1.0e-6;
	vector<double> dst_tetrahedron(6, 1.0);
	dTarget = dst_tetrahedron;
	mol_tetrahedron = new Molecule(dTarget);
	mol_tetrahedron->Add(-0.5, -sqrt(0.75)*1.0/3, 0.0);
	mol_tetrahedron->Add(+0.5, -sqrt(0.75)*1.0/3, 0.0);
	mol_tetrahedron->Add(+0.0, +sqrt(0.75)*2.0/3, 0.0);
	vtx_tetrahedron = new Atom_t(0.0, 0.0, sqrt(2.0/3));
    }

    void tearDown() 
    {
	delete mol_tetrahedron; mol_tetrahedron = NULL;
	delete vtx_tetrahedron; vtx_tetrahedron = NULL;
    }

    void test_RelaxTetrahedron()
    {
	Atom_t vtx(1.0, 2.0, 3.0);
	mol_tetrahedron->RelaxExternalAtom(vtx);
	double dvtx = dist(vtx, *vtx_tetrahedron);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, dvtx, double_eps);
    }

};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(TestRelaxExternalAtom);

// End of file