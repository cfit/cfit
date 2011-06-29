
#include <complex>
#include <cfit/phasespace.hh>
#include <cfit/models/relbreitwigner.hh>

std::complex< double > RelBreitWigner::propagator( const PhaseSpace& ps, const double& mSqAB ) const
{
  const std::complex< double > I( 0., 1. );

  return 1. / ( std::pow( _mass.value(), 2 ) - mSqAB - I * _mass.value() * runningWidth( ps, mSqAB ) );
}
