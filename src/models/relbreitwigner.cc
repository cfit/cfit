
#include <complex>
#include <cfit/phasespace.hh>
#include <cfit/models/relbreitwigner.hh>

std::complex< double > RelBreitWigner::propagator( const PhaseSpace& ps, const double& mSqAB ) const
{
  std::cout << "RBW ";

  const std::complex< double > I( 0., 1. );

  return 1. / ( std::pow( mass(), 2 ) - mSqAB - I * mass() * runningWidth( ps, mSqAB ) );
}


RelBreitWigner* RelBreitWigner::copy() const
{
  return new RelBreitWigner( *this );
}

