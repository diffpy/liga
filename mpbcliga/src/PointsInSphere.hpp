/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2007 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* classes PointsInSphere
*
* Constructors:
*
*     PointsInSphere(Rmin, Rmax, a, b, c, alpha, beta, gamma)
*
*     template<class Lattice> PointsInSphere(Rmin, Rmax, const Lattice&)
*
*     where class Lattice must provide methods a(), b(), c(),
*     alpha(), beta(), gamma()
*
* Examples:
*
*     PointsInSphere sph(Rmin, Rmax, a, b, c, alpha, beta, gamma)
*     for (sph.rewind(); !sph.finished(); sph.next())
*     {
*         // lattice indices are in sph.m(), sph.n(), sph.o() or sph.mno()
*         // sph.r() is distance from origin,
*         // where sph.Rmin() < sph.r() < sph.Rmax()
*     }
*
* Tip: add epsilon to Rmax to avoid roundoff issues
*
***********************************************************************/

#ifndef POINTSINSPHERE_HPP_INCLUDED
#define POINTSINSPHERE_HPP_INCLUDED

// ensure math constants get defined for MSVC
#define _USE_MATH_DEFINES
#include <cmath>


namespace NS_POINTSINSPHERE {

class LatticeParameters
{
    public:

        // data
        // input arguments
        double a, b, c, alpha, beta, gamma;
        // cosines and sines of direct lattice angles
        double ca, cb, cg, sa, sb, sg;
        // reciprocal lattice and its cosines and sines
        double ar, br, cr, alphar, betar, gammar;
        double car, cbr, cgr, sar, sbr, sgr;

        // constructor
        LatticeParameters(double _a, double _b, double _c,
                double _alpha, double _beta, double _gamma);

        // methods
        // calculate all properties from current lattice parameters
        void update();
        // return a reciprocal of this lattice
        LatticeParameters reciprocal() const;

    private:

        // methods
        static double _cosd(double x);
        static double _sind(double x);
};

}       // namespace NS_POINTSINSPHERE


class PointsInSphere
{
    public:

        // constructors
        PointsInSphere(double rmin, double rmax,
                const NS_POINTSINSPHERE::LatticeParameters& _latpar);
        PointsInSphere(double rmin, double rmax,
                double _a, double _b, double _c,
                double _alpha, double _beta, double _gamma);
        template <class L>
            PointsInSphere(double rmin, double rmax, const L&);

        // methods
        // loop control
        void rewind();
        void next();
        bool finished() const;
        // data access
        const double& Rmin() const;
        const double& Rmax() const;
        const int* mno() const;
        const int& m() const;
        const int& n() const;
        const int& o() const;
        double r() const;

    private:

        // data
        // inputs
        const double _Rmin;
        const double _Rmax;
        const NS_POINTSINSPHERE::LatticeParameters latpar;
        // output
        int _mno[3];
        int& _m;
        int& _n;
        int& _o;
        // calculated constants set by init()
        double RminSquare, RmaxSquare;
        // 2D reciprocal parameters and cosine in bc plane
        double b2r, c2r, ca2r;
        // reciprocal c
        double c1r;
        // offset of the nearest point to [0,0,0]
        double dn0dm, do0dm, do0dn;
        // loop variables
        double n0plane, o0plane, o0line;
        double mHalfSpan, nHalfSpan, oHalfSpan;
        // o indices excluded due to Rmin
        double oExclHalfSpan;
        int hi_m, hi_n, hi_o, outside_o;
        double RplaneSquare;

        // methods
        // loop advance
        void next_m();
        void next_n();
        void next_o();
        void init();
};

// template constructor

template <class L>
PointsInSphere::PointsInSphere(double rmin, double rmax, const L& lat) :
    _Rmin(rmin), _Rmax(rmax),
    latpar(lat.a(), lat.b(), lat.c(), lat.alpha(), lat.beta(), lat.gamma()),
    _m(_mno[0]), _n(_mno[1]), _o(_mno[2])
{
    init();
    rewind();
}


#endif  // POINTSINSPHERE_HPP_INCLUDED
