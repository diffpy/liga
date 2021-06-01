/***********************************************************************
* Short Title: declaration of Lattice class
*
* Comments: class for general coordinate system
*
* <license text>
***********************************************************************/

#ifndef LATTICE_HPP_INCLUDED
#define LATTICE_HPP_INCLUDED

#include "R3linalg.hpp"

class Lattice
{
    public:

        // constructors

        Lattice();
        Lattice(double a0, double b0, double c0,
                double alpha0, double beta0, double gamma0);
        // create from base vectors
        template <class V>
            Lattice(const V& va0, const V& vb0, const V& vc0);

        // methods

        // set lattice parameters
        void setLatPar(
                double a0, double b0, double c0,
                double alpha0, double beta0, double gamma0 );
        // set lattice base vectors
        void setLatBase(
                const R3::Vector& va0,
                const R3::Vector& vb0,
                const R3::Vector& vc0 );
        template <class V>
            void setLatBase(const V& va0, const V& vb0, const V& vc0);
        // direct lattice
        double a() const         { return _a; }
        double b() const         { return _b; }
        double c() const         { return _c; }
        double alpha() const     { return _alpha; }
        double beta() const      { return _beta; }
        double gamma() const     { return _gamma; }
        double cosalpha() const  { return cosa; }
        double cosbeta() const   { return cosb; }
        double cosgamma() const  { return cosg; }
        double sinalpha() const  { return sina; }
        double sinbeta() const   { return sinb; }
        double singamma() const  { return sing; }
        const R3::Vector& va() const { return _va; }
        const R3::Vector& vb() const { return _vb; }
        const R3::Vector& vc() const { return _vc; }
        // reciprocal lattice
        double ar() const        { return _ar; }
        double br() const        { return _br; }
        double cr() const        { return _cr; }
        double alphar() const    { return _alphar; }
        double betar() const     { return _betar; }
        double gammar() const    { return _gammar; }
        double cosalphar() const { return cosar; }
        double cosbetar() const  { return cosbr; }
        double cosgammar() const { return cosgr; }
        double sinalphar() const { return sinar; }
        double sinbetar() const  { return sinbr; }
        double singammar() const { return singr; }
        const R3::Vector& var() const    { return _var; }
        const R3::Vector& vbr() const    { return _vbr; }
        const R3::Vector& vcr() const    { return _vcr; }
        // vector operations using lattice coordinates
        template <class V>
            double dot(const V& u, const V& v) const;
        template <class V>
            double norm(const V& u) const;
        template <class V>
            double distance(const V& u, const V& v) const;
        // angle in degrees
        template <class V>
            double angledeg(const V& u, const V& v) const;
        // angle in radians
        template <class V>
            double anglerad(const V& u, const V& v) const;
        // conversion of coordinates and tensors
        const R3::Vector& cartesian(const R3::Vector& lv) const;
        template <class V>
            const R3::Vector& cartesian(const V& lv) const;
        const R3::Vector& fractional(const R3::Vector& cv) const;
        template <class V>
            const R3::Vector& fractional(const V& cv) const;
        const R3::Vector& ucvCartesian(const R3::Vector& cv) const;
        template <class V>
            const R3::Vector& ucvCartesian(const V& cv) const;
        const R3::Vector& ucvFractional(const R3::Vector& lv) const;
        template <class V>
            const R3::Vector& ucvFractional(const V& lv) const;
        /// Cartesian coordinates of an quivalent point nearest to zero
        const R3::Vector& nearZeroCartesian(const R3::Vector& cv) const;
        template <class V>
            const R3::Vector& nearZeroCartesian(const V& cv) const;
        const R3::Matrix& cartesianMatrix(const R3::Matrix& Ml) const;
        const R3::Matrix& fractionalMatrix(const R3::Matrix& Mc) const;
        // largest cell diagonal in fractional coordinates
        R3::Vector ucMaxDiagonal() const;
        double ucMaxDiagonalLength() const;
        // lattice related tensors
        // metrics tensor
        const R3::Matrix& metrics() const   { return _metrics; }
        // matrix of base vectors, base() = stdbase() * baserot()
        const R3::Matrix& base() const      { return _base; }
        // standard base vectors
        const R3::Matrix& stdbase() const   { return _stdbase; }
        // base rotation matrix
        const R3::Matrix& baseRot() const   { return _baserot; }
        // inverse of base matrix
        const R3::Matrix& recbase() const   { return _recbase; }

    private:

        // data - direct lattice parameters
        double _a, _b, _c;
        double _alpha, _beta, _gamma;
        double cosa, cosb, cosg;
        double sina, sinb, sing;
        R3::Vector _va, _vb, _vc;

        // data - reciprocal lattice parameters
        double _ar, _br, _cr;
        double _alphar, _betar, _gammar;
        double cosar, cosbr, cosgr;
        double sinar, sinbr, singr;
        R3::Vector _var, _vbr, _vcr;

        // data - tensors
        // _base = _stdbase * _baserot
        R3::Matrix _metrics;        // metrics tensor
        R3::Matrix _base;           // lattice base
        R3::Matrix _stdbase;        // standard unit cell base
        R3::Matrix _baserot;        // base rotation matrix
        R3::Matrix _recbase;        // inverse of base matrix
        // base multiplied by magnitudes of reciprocal vectors
        R3::Matrix _normbase;
        R3::Matrix _recnormbase;    // inverse of _normbase


};

////////////////////////////////////////////////////////////////////////
// template definitions
////////////////////////////////////////////////////////////////////////

template <class V>
Lattice::Lattice(const V& va0, const V& vb0, const V& vc0)
{
    setLatBase(va0, vb0, vc0);
}

template <class V>
void Lattice::setLatBase(const V& va0, const V& vb0, const V& vc0)
{
    R3::Vector va1(va0[0], va0[1], va0[2]);
    R3::Vector vb1(vb0[0], vb0[1], vb0[2]);
    R3::Vector vc1(vc0[0], vc0[1], vc0[2]);
    setLatBase(va1, vb1, vc1);
}

template <class V>
double Lattice::dot(const V& u, const V& v) const
{
    double dp =
        u[0]*v[0]*_metrics(0,0) +
        u[1]*v[1]*_metrics(1,1) +
        u[2]*v[2]*_metrics(2,2) +
        (u[0]*v[1] + u[1]*v[0])*_metrics(0,1) +
        (u[0]*v[2] + u[2]*v[0])*_metrics(0,2) +
        (u[2]*v[1] + u[1]*v[2])*_metrics(1,2);
    return dp;
}

template <class V>
double Lattice::norm(const V& u) const
{
    return sqrt(dot(u, u));
}

template <class V>
double Lattice::distance(const V& u, const V& v) const
{
    static R3::Vector duv;
    duv[0] = u[0] - v[0];
    duv[1] = u[1] - v[1];
    duv[2] = u[2] - v[2];
    return norm(duv);
}

template <class V>
double Lattice::angledeg(const V& u, const V& v) const
{
    return 180.0/M_PI * anglerad(u, v);
}

template <class V>
double Lattice::anglerad(const V& u, const V& v) const
{
    double ca = dot(u, v)/(norm(u) * norm(v));
    return acos(ca);
}

template <class V>
const R3::Vector& Lattice::cartesian(const V& lv) const
{
    static R3::Vector lvcopy;
    lvcopy = double(lv[0]), double(lv[1]), double(lv[2]);
    return cartesian(lvcopy);
}

template <class V>
const R3::Vector& Lattice::fractional(const V& cv) const
{
    static R3::Vector cvcopy;
    cvcopy = cv[0], cv[1], cv[2];
    return fractional(cvcopy);
}

template <class V>
const R3::Vector& Lattice::ucvCartesian(const V& cv) const
{
    static R3::Vector cvcopy;
    cvcopy = cv[0], cv[1], cv[2];
    return ucvCartesian(cvcopy);
}

template <class V>
const R3::Vector& Lattice::ucvFractional(const V& cv) const
{
    static R3::Vector cvcopy;
    cvcopy = cv[0], cv[1], cv[2];
    return ucvFractional(cvcopy);
}


template <class V>
const R3::Vector& Lattice::nearZeroCartesian(const V& cv) const
{
    static R3::Vector cvcopy;
    cvcopy = cv[0], cv[1], cv[2];
    return this->nearZeroCartesian(cvcopy);
}

#endif  // LATTICE_HPP_INCLUDED
