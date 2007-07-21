/***********************************************************************
* Short Title: unit tests for Lattice class
*
* Comments:
*
* $Id$
*
* <license text>
***********************************************************************/

#include <stdexcept>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "R3linalg.hpp"
#include "Lattice.hpp"

using namespace std;
using namespace blitz;
using R3::MatricesAlmostEqual;
using R3::VectorsAlmostEqual;

class TestLattice : public CppUnit::TestFixture
{

    CPPUNIT_TEST_SUITE(TestLattice);
    CPPUNIT_TEST(test_Lattice);
    CPPUNIT_TEST(test_setLatPar);
    CPPUNIT_TEST(test_setLatBase);
    CPPUNIT_TEST(test_dist);
    CPPUNIT_TEST(test_angle);
    CPPUNIT_TEST_SUITE_END();

private:

    Lattice* lattice;
    static const double precision = 1.0e-12;

public:

    void setUp()
    {
	lattice = new Lattice();
    }

    void tearDown()
    {
	delete lattice;
    }

    void test_Lattice()
    {
	// default lattice should be cartesian
	R3::Matrix identity;
	identity = 1, 0, 0,
		   0, 1, 0,
		   0, 0, 1;
	CPPUNIT_ASSERT(MatricesAlmostEqual(lattice->base(), identity));
	// lattice parameters constructor
	Lattice lattice1(1.0, 2.0, 3.0, 90, 90, 120);
	R3::Vector va = lattice1.va();
	R3::Vector vb = lattice1.vb();
	R3::Vector vc = lattice1.vc();
        double adotb = R3::dot(va, vb);
        double adotc = R3::dot(va, vc);
        double bdotc = R3::dot(vb, vc);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, R3::norm(va), precision);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, R3::norm(vb), precision);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, R3::norm(vc), precision);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5*1.0*2.0, adotb, precision);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, adotc, precision);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, bdotc, precision);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(sqrt(4.0/3), lattice1.ar(), precision);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(sqrt(1.0/3), lattice1.br(), precision);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/3, lattice1.cr(), precision);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(90.0, lattice1.alphar(), precision);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(90.0, lattice1.betar(), precision);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(60.0, lattice1.gammar(), precision);
	// lattice vectors constructor
	va = 1.0, 1.0, 0.0;
	vb = 0.0, 1.0, 1.0;
	vc = 1.0, 0.0, 1.0;
	Lattice lattice2(va, vb, vc);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(sqrt(2.0), lattice2.a(), precision);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(60.0, lattice2.alpha(), precision);
    }

    void test_setLatPar()
    {
        lattice->setLatPar(1.0, 2.0, 3.0, 90, 90, 120);
	R3::Matrix base_check;
	base_check = sqrt(0.75),    -0.5,      0.0,
		     0.0,           +2.0,      0.0,
		     0.0,           +0.0,      3.0;
	CPPUNIT_ASSERT(MatricesAlmostEqual(base_check,
		    lattice->base(), precision) );
	R3::Matrix recbase_check;
        recbase_check = sqrt(4.0/3),    sqrt(1.0/12),   0.0,
                        0.0,            0.5,            0.0,
                        0.0,            0.0,            1.0/3.0;
	CPPUNIT_ASSERT(MatricesAlmostEqual(recbase_check,
		    lattice->recbase(), precision) );
    }

    void test_setLatBase()
    {
        R3::Vector va, vb, vc;
        va = 1.0,  1.0,  0.0;
	vb = 0.0,  1.0,  1.0;
	vc = 1.0,  0.0,  1.0;
        lattice->setLatBase(va, vb, vc);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(sqrt(2.0), lattice->a(), precision);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(sqrt(2.0), lattice->b(), precision);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(sqrt(2.0), lattice->c(), precision);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(60.0, lattice->alpha(), precision);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(60.0, lattice->beta(), precision);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(60.0, lattice->gamma(), precision);
	// check determinant of rotation matrix
	double detRot0 = R3::determinant(lattice->baseRot());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, detRot0, precision);
	// check if rotation matrix works
	R3::Matrix base_check;
	base_check = va[0], va[1], va[2],
		     vb[0], vb[1], vb[2],
		     vc[0], vc[1], vc[2];
	CPPUNIT_ASSERT(MatricesAlmostEqual(base_check,
		    lattice->base(), precision) );
	R3::Matrix recbase_check;
	recbase_check = 0.5,   -0.5,	0.5,
			0.5,    0.5,   -0.5,
		       -0.5,	0.5,	0.5;
	CPPUNIT_ASSERT(MatricesAlmostEqual(recbase_check,
		    lattice->recbase(), precision) );
        lattice->setLatPar( lattice->a(), lattice->b(), lattice->c(),
		44.0, 66.0, 88.0 );
	CPPUNIT_ASSERT(!MatricesAlmostEqual(base_check,
		    lattice->base(), precision) );
	CPPUNIT_ASSERT(!MatricesAlmostEqual(recbase_check,
		    lattice->recbase(), precision) );
        lattice->setLatPar( lattice->a(), lattice->b(), lattice->c(),
		60.0, 60.0, 60.0 );
	CPPUNIT_ASSERT(MatricesAlmostEqual(base_check,
		    lattice->base(), precision) );
	CPPUNIT_ASSERT(MatricesAlmostEqual(recbase_check,
		    lattice->recbase(), precision) );
    }

    void test_dist()
    {
	R3::Vector va, vb;
	va = 1.0, 2.0, 2.0;
	vb = 0.0, 0.0, 0.0;
	CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, lattice->dist(va, vb), precision);
	lattice->setLatPar(2.0, 2.0, 2.0, 90.0, 90.0, 90.0);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(6.0, lattice->dist(va, vb), precision);
    }

    void test_angle()
    {
	R3::Vector va, vb;
	va = 1.0, 0.0, 0.0;
	vb = 0.0, 1.0, 0.0;
	CPPUNIT_ASSERT_DOUBLES_EQUAL(90.0,
		lattice->angledeg(va, vb), precision);
	lattice->setLatPar(2.0, 2.0, 2.0, 90.0, 90.0, 120.0);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(120.0,
		lattice->angledeg(va, vb), precision);
    }

};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(TestLattice);

// End of file
