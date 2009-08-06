#include "/home/juhas/arch/i686/include/dbprint.h"
/***********************************************************************
* Short Title: unit tests of least squares atom relaxation
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#include <cxxtest/TestSuite.h>

#include "AtomCost.hpp"
#include "AtomOverlapCost.hpp"
#include "LigaUtils.hpp"
#include "Crystal.hpp"

using namespace std;

class TestRelaxCrystalAtom : public CxxTest::TestSuite
{
    private:

        double eps_distance;
        auto_ptr<DistanceTable> dtgt_bcc;
        auto_ptr<DistanceTable> dtgt_fcc;
        auto_ptr<DistanceTable> dtgt_triclinic;
        auto_ptr<Crystal> crbcc;
        auto_ptr<Crystal> crfcc;
        auto_ptr<Crystal> crtriclinic;

        template <class T> T* newFromFile(const string& filename)
        {
            T instance;
            ifstream fid(filename.c_str());
            fid >> instance;
            T* rv = new T(instance); 
            return rv;
        }

    public:

        void setUp()
        {
            eps_distance = NS_LIGA::eps_distance;
            // load distance tables
            DistanceTable* ptbl;
            ptbl = this->newFromFile<DistanceTable>("crystals/bcc.dst");
            this->dtgt_bcc.reset(ptbl);
            ptbl = this->newFromFile<DistanceTable>("crystals/fcc.dst");
            this->dtgt_fcc.reset(ptbl);
            ptbl = this->newFromFile<DistanceTable>("crystals/triclinic.dst");
            this->dtgt_triclinic.reset(ptbl);
            // bcc
            this->crbcc.reset(new Crystal);
            this->crbcc->setDistanceTable(*this->dtgt_bcc);
            this->crbcc->ReadFile("crystals/bcc.stru");
            // fcc
            this->crfcc.reset(new Crystal);
            this->crfcc->setDistanceTable(*this->dtgt_fcc);
            this->crfcc->ReadFile("crystals/fcc.stru");
            // triclinic
            this->crtriclinic.reset(new Crystal);
            this->crtriclinic->setDistanceTable(*this->dtgt_triclinic);
            this->crtriclinic->ReadFile("crystals/triclinic.stru");
        }


        void tearDown() 
        { 
            this->crbcc->getAtomCostCalculator()->setScale(1.0);
            this->crbcc->getAtomOverlapCalculator()->setScale(1.0);
        }


        void test_RelaxBCC()
        {
            Atom_t a1 = this->crbcc->getAtom(1);
            this->crbcc->Pop(1);
            Atom_t arx = a1;
            R3::Vector offset(0.013, -0.07, -0.03);
            arx.r += offset;
            this->crbcc->RelaxExternalAtom(&arx);
            double drx = R3::distance(a1.r, arx.r);
            TS_ASSERT_DELTA(0.0, drx, eps_distance);
        }


        void test_RelaxFCC()
        {
            Atom_t a1 = this->crfcc->getAtom(1);
            this->crfcc->Pop(1);
            Atom_t arx = a1;
            R3::Vector offset(0.013, -0.07, -0.03);
            arx.r += offset;
            this->crfcc->RelaxExternalAtom(&arx);
            double drx = R3::distance(a1.r, arx.r);
            TS_ASSERT_DELTA(0.0, drx, eps_distance);
        }


        void test_RelaxTriclinic()
        {
            Atom_t a1 = this->crtriclinic->getAtom(1);
            this->crtriclinic->Pop(1);
            Atom_t arx = a1;
            R3::Vector offset(0.0013, -0.007, -0.003);
            arx.r += offset;
            this->crtriclinic->RelaxExternalAtom(&arx);
            double drx = R3::distance(a1.r, arx.r);
            TS_ASSERT_DELTA(0.0, drx, eps_distance);
        }


        void test_RelaxOverlapBCC()
        {
            this->crbcc->getAtomCostCalculator()->setScale(0.0);
            DBPRINT(crbcc->pairsPerAtom());
            this->crbcc->recalculate();
            TS_ASSERT_DELTA(0.0, this->crbcc->cost(), eps_distance);
            AtomRadiiTable radii;
            radii["C"] = 0.5;
            this->crbcc->fetchAtomRadii(radii);
            TS_ASSERT(this->crbcc->cost() > eps_distance);
            Atom_t a1 = this->crbcc->getAtom(1);
            DBPRINT(crbcc->countPairs());
            DBPRINT(crbcc->countAtoms());
            DBPRINT(crbcc->pairsPerAtom());
            this->crbcc->Pop(1);
            DBPRINT(crbcc->pairsPerAtomInc());
            DBPRINT(crbcc->countPairs());
            DBPRINT(crbcc->countAtoms());
            DBPRINT(crbcc->pairsPerAtom());
            Atom_t arx = a1;
            R3::Vector offset(0.013, -0.07, -0.03);
            arx.r += offset;
            this->crbcc->RelaxExternalAtom(&arx);
            double drx = R3::distance(a1.r, arx.r);
            TS_ASSERT_DELTA(0.0, drx, eps_distance);
        }

};  // class TestRelaxCrystalAtom

// End of file