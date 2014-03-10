
#include <complex>
#include <cfit/phasespace.hh>
#include <cfit/models/flatte.hh>



std::complex< double > Flatte::propagator( const PhaseSpace& ps, const double& mSqAB ) const
{
  const std::complex< double > I( 0., 1. );

  const double& mGamma0 = mGamma();

  const double& rho10 = rho( ps, mSq() );
  const double& rho1  = rho( ps, mSqAB );
  const double& g1    = gamma1Sq() * rho1 / rho10;

  const double& rho20 = std::sqrt( kallen( mSq(), m02aSq(), m02bSq() ) ) / mSqAB;
  const double& rho2  = std::sqrt( kallen( mSqAB, m02aSq(), m02bSq() ) ) / mSqAB;
  const double& g2    = gamma2Sq() * rho2 / rho20;

  return mGamma0 * gamma1Sq() / ( mSq() - mSqAB - I * mGamma0 * ( g1 + g2 ) * std::pow( blattWeisskopf( ps, mSqAB ), 2 ) );
}


Flatte* Flatte::copy() const
{
  return new Flatte( *this );
}

