
#include <complex>
#include <cfit/phasespace.hh>
#include <cfit/models/glass.hh>


std::complex< double > GLass::propagator( const PhaseSpace& ps, const double& mSqAB ) const
{
  const std::complex< double > I( 0., 1. );

  double qAB   = q  ( ps, mSqAB );
  double qSqAB = qSq( ps, mSqAB );

  double qCotDeltaB = 1. / lassa() + lassr() * qSqAB / 2.;
  double cotDeltaB  = qCotDeltaB / qAB;

  std::complex< double > rTerm = lassR() * std::exp( I * phiR() + 2. * I * phiB() );
  rTerm *= ( qCotDeltaB + I * qAB ) / ( qCotDeltaB - I * qAB );
  rTerm *= m() * width() / ( mSq() - mSqAB - I * m() * runningWidth( ps, mSqAB ) );

  std::complex< double > bTerm = lassB() * std::sqrt( mSqAB ) / 2. * std::exp( I * phiB() );
  bTerm *= ( std::cos( phiB() ) + std::sin( phiB() ) * cotDeltaB ) / ( qCotDeltaB - I * qAB );

  return rTerm + bTerm;
}


GLass* GLass::copy() const
{
  return new GLass( *this );
}

