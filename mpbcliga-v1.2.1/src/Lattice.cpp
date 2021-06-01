/***********************************************************************
* Short Title: definition of Lattice class
*
* Comments: class for general coordinate system
*
* <license text>
***********************************************************************/

#include <list>
#include <stdexcept>
#include "Lattice.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////
// helper functions
////////////////////////////////////////////////////////////////////////

namespace NS_LATTICE {

double cosd(double x)
{
    double xp = fmod(fabs(x), 360.0);
    if (remainder(xp, 60.0) == 0.0 || remainder(xp, 90.0) == 0.0)
    {
        switch(int(round(xp)))
        {
            case 0: return 1.0;
            case 60:
            case 300: return 0.5;
            case 90:
            case 270: return 0.0;
            case 120:
            case 240: return -0.5;
            case 180: return -1.0;
        };
    }
    return cos(x/180.0*M_PI);
}

double sind(double x)
{
    return NS_LATTICE::cosd(90.0 - x);
}

double acosd(double x)
{
    if (remainder(x, 0.5) == 0.0)
    {
        switch(int(round(x/0.5)))
        {
            case 0: return 90.0;
            case 1: return 60.0;
            case -1: return 120.0;
            case 2: return 0.0;
            case -2: return 180.0;
        };
    }
    return acos(x)/M_PI*180.0;
}

} // namespace NS_LATTICE

////////////////////////////////////////////////////////////////////////
// Lattice definitions
////////////////////////////////////////////////////////////////////////

Lattice::Lattice()
{
    _baserot = 1.0, 0.0, 0.0,
               0.0, 1.0, 0.0,
               0.0, 0.0, 1.0;
    setLatPar(1.0, 1.0, 1.0, 90.0, 90.0, 90.0);
}

Lattice::Lattice( double a0, double b0, double c0,
            double alpha0, double beta0, double gamma0 )
{
    _baserot = 1.0, 0.0, 0.0,
               0.0, 1.0, 0.0,
               0.0, 0.0, 1.0;
    setLatPar(a0, b0, c0, alpha0, beta0, gamma0);
}

void Lattice::setLatPar( double a0, double b0, double c0,
            double alpha0, double beta0, double gamma0 )
{
    using namespace NS_LATTICE;
    _a = a0; _b = b0; _c = c0;
    _alpha = alpha0; _beta = beta0; _gamma = gamma0;
    cosa = NS_LATTICE::cosd(_alpha);    sina = NS_LATTICE::sind(_alpha);
    cosb = NS_LATTICE::cosd(_beta);     sinb = NS_LATTICE::sind(_beta);
    cosg = NS_LATTICE::cosd(_gamma);    sing = NS_LATTICE::sind(_gamma);
    // Vunit is a volume of unit cell with a=b=c=1
    double Vunit = sqrt(
            1.0 + 2.0*cosa*cosb*cosg - cosa*cosa - cosb*cosb - cosg*cosg);
    // reciprocal lattice
    _ar = sina/(_a*Vunit);
    _br = sinb/(_b*Vunit);
    _cr = sing/(_c*Vunit);
    cosar = (cosb*cosg - cosa)/(sinb*sing);
    cosbr = (cosa*cosg - cosb)/(sina*sing);
    cosgr = (cosa*cosb - cosg)/(sina*sinb);
    sinar = sqrt(1.0 - cosar*cosar);
    sinbr = sqrt(1.0 - cosbr*cosbr);
    singr = sqrt(1.0 - cosgr*cosgr);
    _alphar = NS_LATTICE::acosd(cosar);
    _betar = NS_LATTICE::acosd(cosbr);
    _gammar = NS_LATTICE::acosd(cosgr);
    // metric tensor
    _metrics = _a*_a,       _a*_b*cosg,   _a*_c*cosb,
               _b*_a*cosg,  _b*_b,        _b*_c*cosa,
               _c*_a*cosb,  _c*_b*cosa,   _c*_c;

    // standard cartesian coordinates of lattice vectors
    _stdbase = 1.0/_ar,     -cosgr/singr/_ar,   cosb*_a ,
               0.0,         _b*sina,            _b*cosa ,
               0.0,         0.0,                _c;
    // calculate unit cell rotation matrix _baserot,  _base = _stdbase*_baserot
    _base = R3::product(_stdbase, _baserot);
    _va = _base(0,0), _base(0,1), _base(0,2);
    _vb = _base(1,0), _base(1,1), _base(1,2);
    _vc = _base(2,0), _base(2,1), _base(2,2);
    _recbase = R3::inverse(_base);
    _var = _recbase(0,0), _recbase(1,0), _recbase(2,0);
    _vbr = _recbase(0,1), _recbase(1,1), _recbase(2,1);
    _vcr = _recbase(0,2), _recbase(1,2), _recbase(2,2);
    _normbase =
        _base(0,0)*_ar,     _base(0,1)*_ar,     _base(0,2)*_ar,
        _base(1,0)*_br,     _base(1,1)*_br,     _base(1,2)*_br,
        _base(2,0)*_cr,     _base(2,1)*_cr,     _base(2,2)*_cr;
    _recnormbase =
        _recbase(0,0)/_ar,  _recbase(0,1)/_br,  _recbase(0,2)/_cr,
        _recbase(1,0)/_ar,  _recbase(1,1)/_br,  _recbase(1,2)/_cr,
        _recbase(2,0)/_ar,  _recbase(2,1)/_br,  _recbase(2,2)/_cr;
}

void Lattice::setLatBase(const R3::Vector& va0,
        const R3::Vector& vb0,
        const R3::Vector& vc0)
{
    using namespace NS_LATTICE;
    _base = va0[0], va0[1], va0[2],
            vb0[0], vb0[1], vb0[2],
            vc0[0], vc0[1], vc0[2];
    double detbase = R3::determinant(_base);
    if (fabs(detbase) < 1.0e-8)
    {
        throw invalid_argument("base vectors are degenerate");
    }
    else if (detbase < 0.0)
    {
        throw invalid_argument("base is not right-handed");
    }
    _va = va0;
    _vb = vb0;
    _vc = vc0;
    _a = R3::norm(_va);
    _b = R3::norm(_vb);
    _c = R3::norm(_vc);
    cosa = R3::dot(_vb, _vc) / (_b*_c);
    cosb = R3::dot(_va, _vc) / (_a*_c);
    cosg = R3::dot(_va, _vb) / (_a*_b);
    sina = sqrt(1.0 - cosa*cosa);
    sinb = sqrt(1.0 - cosb*cosb);
    sing = sqrt(1.0 - cosg*cosg);
    _alpha = NS_LATTICE::acosd(cosa);
    _beta = NS_LATTICE::acosd(cosb);
    _gamma = NS_LATTICE::acosd(cosg);
    // Vunit is a volume of unit cell with a=b=c=1
    double Vunit = sqrt(
            1.0 + 2.0*cosa*cosb*cosg - cosa*cosa - cosb*cosb - cosg*cosg);
    // reciprocal lattice
    _ar = sina/(_a*Vunit);
    _br = sinb/(_b*Vunit);
    _cr = sing/(_c*Vunit);
    cosar = (cosb*cosg - cosa)/(sinb*sing);
    cosbr = (cosa*cosg - cosb)/(sina*sing);
    cosgr = (cosa*cosb - cosg)/(sina*sinb);
    sinar = sqrt(1.0 - cosar*cosar);
    sinbr = sqrt(1.0 - cosbr*cosbr);
    singr = sqrt(1.0 - cosgr*cosgr);
    _alphar = NS_LATTICE::acosd(cosar);
    _betar = NS_LATTICE::acosd(cosbr);
    _gammar = NS_LATTICE::acosd(cosgr);
    // standard orientation of lattice vectors
    _stdbase = 1.0/_ar,     -cosgr/singr/_ar,   cosb*_a,
               0.0,         _b*sina,            _b*cosa,
               0.0,         0.0,                _c;
    // calculate unit cell rotation matrix,  base = stdbase*baserot
    _baserot = R3::product(R3::inverse(_stdbase), _base);
    _recbase = R3::inverse(_base);
    _var = _recbase(0,0), _recbase(1,0), _recbase(2,0);
    _vbr = _recbase(0,1), _recbase(1,1), _recbase(2,1);
    _vcr = _recbase(0,2), _recbase(1,2), _recbase(2,2);
    // bases normalized to unit reciprocal vectors
    _normbase =
        _base(0,0)*_ar,     _base(0,1)*_ar,     _base(0,2)*_ar,
        _base(1,0)*_br,     _base(1,1)*_br,     _base(1,2)*_br,
        _base(2,0)*_cr,     _base(2,1)*_cr,     _base(2,2)*_cr;
    _recnormbase =
        _recbase(0,0)/_ar,  _recbase(0,1)/_br,  _recbase(0,2)/_cr,
        _recbase(1,0)/_ar,  _recbase(1,1)/_br,  _recbase(1,2)/_cr,
        _recbase(2,0)/_ar,  _recbase(2,1)/_br,  _recbase(2,2)/_cr;
    // update metrics tensor
    _metrics = _a*_a,     _a*_b*cosg,  _a*_c*cosb,
               _b*_a*cosg,  _b*_b,     _b*_c*cosa,
               _c*_a*cosb,  _c*_b*cosa,  _c*_c;
}

const R3::Vector& Lattice::cartesian(const R3::Vector& lv) const
{
    static R3::Vector res;
    res = R3::product(lv, _base);
    return res;
}

const R3::Vector& Lattice::fractional(const R3::Vector& cv) const
{
    static R3::Vector res;
    res = R3::product(cv, _recbase);
    return res;
}

const R3::Vector& Lattice::ucvCartesian(const R3::Vector& cv) const
{
    static R3::Vector res;
    res = cartesian(ucvFractional(fractional(cv)));
    return res;
}

const R3::Vector& Lattice::ucvFractional(const R3::Vector& lv) const
{
    static R3::Vector res;
    res = lv - floor(lv);
    return res;
}


const R3::Vector& Lattice::nearZeroCartesian(const R3::Vector& cv) const
{
    static R3::Vector rv;
    rv = this->ucvCartesian(cv);
    R3::Vector cv1 = rv;
    R3::Vector cv2 = rv;
    double rvsquare = R3::dot(rv, rv);
    for (int x0 = -1; x0 < 1; ++x0)
        for (int y0 = -1; y0 < 1; ++y0)
            for (int z0 = -1; z0 < 1; ++z0)
            {
                cv2 = cv1 + x0 * this->va() + y0 * this->vb() + z0 * this->vc();
                double cv2square = R3::dot(cv2, cv2);
                if (cv2square < rvsquare)
                {
                    rv = cv2;
                    rvsquare = cv2square;
                }
            }
    return rv;
}


const R3::Matrix& Lattice::cartesianMatrix(const R3::Matrix& Ml) const
{
    using R3::product;
    static R3::Matrix res;
    res = product(Ml, _normbase);
    res = product(R3::transpose(_normbase), res);
    return res;
}

const R3::Matrix& Lattice::fractionalMatrix(const R3::Matrix& Mc) const
{
    static R3::Matrix res0, res1;
    res0 = R3::product(Mc, _recnormbase);
    res1 = R3::product(R3::transpose(_recnormbase), res0);
    return res1;
}


R3::Vector Lattice::ucMaxDiagonal() const
{
    static list<R3::Vector> ucdiagonals;
    if (ucdiagonals.empty())
    {
        ucdiagonals.push_back(R3::Vector(+1, +1, +1));
        ucdiagonals.push_back(R3::Vector(-1, +1, +1));
        ucdiagonals.push_back(R3::Vector(+1, -1, +1));
        ucdiagonals.push_back(R3::Vector(+1, +1, -1));
    }
    double maxnorm = -1;
    list<R3::Vector>::iterator ucd;
    list<R3::Vector>::iterator maxucd = ucdiagonals.end();
    for (ucd = ucdiagonals.begin(); ucd != ucdiagonals.end(); ++ucd)
    {
        double normucd = this->norm(*ucd);
        if (normucd > maxnorm)
        {
            maxnorm = normucd;
            maxucd = ucd;
        }
    }
    assert(maxucd != ucdiagonals.end());
    return *maxucd;
}


double Lattice::ucMaxDiagonalLength() const
{
    R3::Vector ucmd = this->ucMaxDiagonal();
    double res = this->norm(ucmd);
    return res;
}
