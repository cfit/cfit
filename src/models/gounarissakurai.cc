
#include <complex>
#include <cfit/phasespace.hh>
#include <cfit/models/gounarissakurai.hh>

std::complex< double > GounarisSakurai::propagator( const PhaseSpace& ps, const double& mSqAB ) const
{
  const std::complex< double > I( 0., 1. );

  double gsd = ( 3. / M_PI ) * ( ps.mSq( _resoA ) / qSq( ps, mSq() ) );
  gsd *= log( ( m() + 2. * q( ps, mSq() ) ) / ( 2. * ps.m( _resoA ) ) );
  gsd += m() / ( 2. * M_PI * q( ps, mSq() ) ) - ( ps.mSq( _resoA ) * m() ) / ( M_PI * std::pow( q( ps, mSq() ), 3 ) );

  std::complex< double > prop = 1. + gsd * width() / mass();

  prop *= 1. / ( std::pow( mass(), 2 ) - mSqAB + gsf( ps, mSqAB ) -I * mass() * runningWidth( ps, mSqAB ) );

  return prop;
}


GounarisSakurai* GounarisSakurai::copy() const
{
  return new GounarisSakurai( *this );
}


double GounarisSakurai::gsf( const PhaseSpace& ps, const double& mSq12 ) const
{
  double factor = width() * std::pow( mass(), 2 ) / q( ps, mSq() );

  double first  = ( qSq( ps, mSq12 ) / qSq( ps, mSq() ) ) * ( gsh( ps, mSq12 ) - gsh( ps, mSq() ) );
  double second = ( mSq() - mSq12 ) * gshprime( ps, mSq() );

  return factor * ( first + second );
}


double GounarisSakurai::gsh( const PhaseSpace& ps, const double& mSq12 ) const
{
  double m12 = std::sqrt(     mSq12 );
  double p12 = q        ( ps, mSq12 );

  return ( 2. / M_PI ) * ( p12 / m12 ) * log( ( m12 + 2. * p12 ) / ( 2. * ps.m1() ) );
}

double GounarisSakurai::gshprime( const PhaseSpace& ps, const double& mSq12 ) const
{
  double first  = 1. / ( 8. * qSq( ps, mSq12 ) );
  double second = 1. / ( 2. * mSq12        );

  return ( first - second ) * gsh( ps, mSq12 ) + second / M_PI;
}

