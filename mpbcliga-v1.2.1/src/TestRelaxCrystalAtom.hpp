/***********************************************************************
* Short Title: unit tests of least squares atom relaxation
*
* Comments:
*
* <license text>
***********************************************************************/

#include <cxxtest/TestSuite.h>

#include "AtomCost.hpp"
#include "AtomOverlap.hpp"
#include "LigaUtils.hpp"
#include "Crystal.hpp"
#include "tests_dir.hpp"

using namespace std;

class TestRelaxCrystalAtom : public CxxTest::TestSuite
{
    private:

        double eps_distance;
        double gradient_eps;
        Crystal crbcc;
        Crystal crfcc;
        Crystal crtriclinic;

        template <class T> T* newFromFile(const string& filename)
        {
            T instance;
            string fullfile = prepend_tests_dir(filename);
            ifstream fid(fullfile.c_str());
            fid >> instance;
            T* rv = new T(instance);
            return rv;
        }

    public:

        void setUp()
        {
            gradient_eps = 5e-6;
            eps_distance = NS_LIGA::eps_distance;
            // load distance tables
            static bool dataloaded = false;
            static Crystal cached_crbcc, cached_crfcc, cached_crtriclinic;
            if (!dataloaded)
            {
                auto_ptr<DistanceTable> dtbl;
                dtbl.reset(newFromFile<DistanceTable>("crystals/bcc.dst"));
                // bcc
                cached_crbcc.setDistanceTable(*dtbl);
                cached_crbcc.ReadFile(prepend_tests_dir("crystals/bcc.stru"));
                // fcc
                dtbl.reset(newFromFile<DistanceTable>("crystals/fcc.dst"));
                cached_crfcc.setDistanceTable(*dtbl);
                cached_crfcc.ReadFile(prepend_tests_dir("crystals/fcc.stru"));
                // triclinic
                dtbl.reset(
                        newFromFile<DistanceTable>("crystals/triclinic.dst"));
                cached_crtriclinic.setDistanceTable(*dtbl);
                cached_crtriclinic.ReadFile(
                        prepend_tests_dir("crystals/triclinic.stru"));
                dataloaded = true;
            }
            this->crbcc = cached_crbcc;
            this->crfcc = cached_crfcc;
            this->crtriclinic = cached_crtriclinic;
        }


        void tearDown()
        {
            this->crbcc.getAtomCostCalculator()->setScale(1.0);
            this->crbcc.getAtomOverlapCalculator()->setScale(1.0);
        }


        void test_RelaxBCC()
        {
            Atom_t a1 = this->crbcc.getAtom(1);
            this->crbcc.Pop(1);
            Atom_t arx = a1;
            R3::Vector offset(0.013, -0.07, -0.03);
            arx.r += offset;
            this->crbcc.RelaxExternalAtom(&arx);
            double drx = R3::distance(a1.r, arx.r);
            TS_ASSERT_DELTA(0.0, drx, eps_distance);
        }


        void test_RelaxFCC()
        {
            Atom_t a1 = this->crfcc.getAtom(1);
            this->crfcc.Pop(1);
            Atom_t arx = a1;
            R3::Vector offset(0.013, -0.07, -0.03);
            arx.r += offset;
            this->crfcc.RelaxExternalAtom(&arx);
            double drx = R3::distance(a1.r, arx.r);
            TS_ASSERT_DELTA(0.0, drx, eps_distance);
        }


        void test_RelaxTriclinic()
        {
            Atom_t a1 = this->crtriclinic.getAtom(1);
            this->crtriclinic.Pop(1);
            Atom_t arx = a1;
            R3::Vector offset(0.0013, -0.007, -0.003);
            arx.r += offset;
            this->crtriclinic.RelaxExternalAtom(&arx);
            double drx = R3::distance(a1.r, arx.r);
            TS_ASSERT_DELTA(0.0, drx, eps_distance);
        }


        void test_RelaxOverlapBCC()
        {
            this->crbcc.getAtomCostCalculator()->setScale(0.0);
            this->crbcc.recalculate();
            TS_ASSERT_DELTA(0.0, this->crbcc.cost(), eps_distance);
            this->crbcc.setAtomRadiiTable("C:0.5");
            TS_ASSERT(this->crbcc.cost() > eps_distance);
            Atom_t a1 = this->crbcc.getAtom(1);
            this->crbcc.Pop(1);
            Atom_t arx = a1;
            R3::Vector offset(0.013, -0.07, -0.03);
            arx.r += offset;
            this->crbcc.RelaxExternalAtom(&arx);
            double drx = R3::distance(a1.r, arx.r);
            TS_ASSERT_DELTA(0.0, drx, eps_distance);
        }


        void test_RelaxOverlapGradientFCC()
        {
            AtomRadiiTable radii;
            radii["C"] = sqrt(2.0) / 4;
            this->crfcc.setAtomRadiiTable(radii);
            TS_ASSERT_DELTA(0.0, this->crfcc.cost(), eps_distance);
            Atom_t a3 = this->crfcc.getAtom(3);
            this->crfcc.Pop(3);
            R3::Vector offset(0.013, -0.07, -0.03);
            a3.r += offset;
            R3::Vector ga, gn;
            ga = crfcc.rxaCheckGradient(&a3);
            R3::Vector gnrxa = numgrad_rxaCheckCost(crfcc, &a3);
            double gdiff = R3::norm((gnrxa - ga) / R3::norm(ga));
            TS_ASSERT_DELTA(0.0, gdiff, gradient_eps);
            crfcc.Add(a3);
            gn = this->numgrad_internal(crfcc, 3);
            // ga is not normalized by number of pairs
            gdiff = R3::norm(gn / R3::norm(gn) - ga / R3::norm(ga));
            TS_ASSERT_DELTA(0.0, gdiff, gradient_eps);
            TS_ASSERT(R3::norm(gn) > gradient_eps*10);
            // compare numerical gradient from rxaCheckCost();
            ga = crfcc.rxaCheckGradient(&a3);
        }

    private:

        R3::Vector numgrad_internal(Molecule& mol, int idx)
        {
            const double delta = 1e-8;
            double c0 = mol.cost();
            R3::Vector numgrad;
            for (int i = 0; i < 3; ++i)
            {
                auto_ptr<Molecule> molmod(mol.copy());
                Atom_t ai = molmod->getAtom(idx);
                ai.r[i] += delta;
                molmod->Pop(idx);
                molmod->Add(ai);
                double c1 = molmod->cost();
                numgrad[i] = (c1 - c0) / delta;
            }
            return numgrad;
        }


        R3::Vector numgrad_rxaCheckCost(Molecule& mol, const Atom_t* pa)
        {
            const double delta = 1e-8;
            double c0 = mol.rxaCheckCost(pa);
            R3::Vector numgrad;
            for (int i = 0; i < 3; ++i)
            {
                Atom_t ai = *pa;
                ai.r[i] += delta;
                double c1 = mol.rxaCheckCost(&ai);
                numgrad[i] = (c1 - c0) / delta;
            }
            return numgrad;
        }


};  // class TestRelaxCrystalAtom

// End of file
