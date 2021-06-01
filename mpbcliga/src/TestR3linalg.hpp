/***********************************************************************
* Short Title: unit tests for R3linalg
*
* Comments:
*
* <license text>
***********************************************************************/

#include <cxxtest/TestSuite.h>

#include "R3linalg.hpp"

using namespace std;
using R3::MatricesAlmostEqual;
using R3::VectorsAlmostEqual;

class TestR3linalg : public CxxTest::TestSuite
{
    private:

        double precision;

    public:

        void setUp()
        {
            precision = 1.0e-12;
        }


        void test_determinant()
        {
            // default lattice should be cartesian
            R3::Matrix A1, A2;
            A1 = 9, 1, 9,
               6, 7, 4,
               0, 2, 9;
            A2 = 9, 1, 9,
               0, 2, 9,
               6, 7, 4;
            double detA = 549;
            TS_ASSERT_EQUALS(detA, R3::determinant(A1));
            TS_ASSERT_EQUALS(-detA, R3::determinant(A2));
        }


        void test_inverse()
        {
            R3::Matrix A, invA;
            A =
                0.5359, -0.5904, 0.8670,
                -0.0053,  0.7559, 0.2692,
                -0.8926,  0.9424, 0.9692;
            invA =
                0.49063197005867, 1.42323870111089, -0.83420736316541,
                -0.24089965988852, 1.32489393619466, -0.15249839300481,
                0.68609361943181, 0.02249568627913,  0.41178393851247;
            TS_ASSERT(MatricesAlmostEqual(invA, R3::inverse(A), precision));
        }


        void test_norm()
        {
            R3::Vector v1, v2;
            v1 = 3.0, 4.0, 0.0;
            v2 = 0.67538115798129, 0.72108424545413, 0.15458914063315;
            TS_ASSERT_DELTA(5.0, R3::norm(v1), precision);
            TS_ASSERT_DELTA(1.0, R3::norm(v2), precision);
        }


        void test_dot()
        {
            R3::Vector v1, v2;
            v1 = -0.97157650177843, 0.43206192654604, 0.56318686427062;
            v2 = -0.04787719419083, 0.55895824010234, -0.34472910285751;
            double dot_v1v2 = 0.09387402846316;
            TS_ASSERT_DELTA(dot_v1v2, R3::dot(v1, v2), precision);
        }


        void test_cross()
        {
            R3::Vector v1, v2, v1xv2;
            v1 = -0.55160549932839, -0.58291452407504,  0.12378162306543;
            v2 =  0.60842511285200, -0.97946444006248, -0.02828214306095;
            v1xv2 = 0.13772577008800,  0.05971126233738,  0.89493780662849;
            TS_ASSERT( VectorsAlmostEqual(v1xv2,
                        R3::cross(v1, v2), precision) );
        }

};  // class TestR3linalg

// End of file
