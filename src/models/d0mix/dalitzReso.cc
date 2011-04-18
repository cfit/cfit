
#include <iostream>
#include <cstdlib>
#include <complex>

#include <cfit/models/d0mix/dalitzReso.hh>


DalitzReso::DalitzReso( char pair1, char pair2, int L, double mM, double mA, double mB, double mC )
  : _pair1( pair1 ), _pair2( pair2 ), _L( L ), _R( 0. ),
    _lassa( 1. ), _lassr( 0. ), _lassB( 0. ), _lassPhiB( 0. ), _lassR( 0. ), _lassPhiR( 0. ) // Default LASS (return 0).
{
  // Get the non resonant index from the resonant ones.
  _other = 'A' + 'B' + 'C' - _pair1 - _pair2;

  // Squared mass of the mother.
  _mSqM = std::pow( mM, 2 );

  // Squared masses of the daughters. _pair1 and _pair2 are the resonant ones.
  switch ( _pair1 )
    {
    case 'A':
      _mSq1 = std::pow( mA, 2 );
      _m1 = mA; break;
    case 'B':
      _mSq1 = std::pow( mB, 2 );
      _m1 = mB; break;
    case 'C':
      _mSq1 = std::pow( mC, 2 );
      _m1 = mC; break;
    default:
      std::cerr << "Resonant particles can only be named \'A\', \'B\' or \'C\'." << std::endl;
      exit( -1 );
    }

  switch ( _pair2 )
    {
    case 'A':
      _mSq2 = std::pow( mA, 2 ); break;
    case 'B':
      _mSq2 = std::pow( mB, 2 ); break;
    case 'C':
      _mSq2 = std::pow( mC, 2 ); break;
    default:
      std::cerr << "Resonant particles can only be named \'A\', \'B\' or \'C\'." << std::endl;
      exit( -1 );
    }

  switch ( _other )
    {
    case 'A':
      _mSq3 = std::pow( mA, 2 ); break;
    case 'B':
      _mSq3 = std::pow( mB, 2 ); break;
    case 'C':
      _mSq3 = std::pow( mC, 2 ); break;
    default:
      std::cerr << "Resonant particles can only be named \'A\', \'B\' or \'C\'." << std::endl;
      exit( -1 );
    }
}



void DalitzReso::setReso( double m0, double g0, double R )
{
  _m0 = m0;
  _g0 = g0;
  _R  = R;

  // Squared mass of the resonant particle.
  _mSq0 = std::pow( m0, 2 );

  // rho = R^2 * p^2, necessary to compute the Blatt-Weisskopf factors.
  _rho0 = rho( _mSq0 );
  _p0   = p  ( _mSq0 );
  _pSq0 = p2 ( _mSq0 );

  // Gounaris-Sakurai constant d factor.
  _gsd = ( 3. / M_PI ) * ( _mSq1 / _pSq0 ) * log( ( _m0 + 2. * _p0 ) / ( 2. * _m1 ) );
  _gsd += _m0 / ( 2. * M_PI * _p0 ) - ( _mSq1 * _m0 ) / ( M_PI * std::pow( _p0, 3 ) );
}














double DalitzReso::kallen( double x, double y, double z )
{
  double result = 0.;
  result += std::pow( x, 2 );
  result += std::pow( y, 2 );
  result += std::pow( z, 2 );
  result -= 2. * x * y;
  result -= 2. * x * z;
  result -= 2. * y * z;

  return result;
}


double DalitzReso::blattWeisskopf( const double& m2 ) const
{
  if ( _L == 0 )
    return 1.;

  if ( _L == 1 )
    return std::sqrt( ( 1. + _rho0 ) / ( 1. + rho( m2 ) ) );

  if ( _L == 2 )
    {
      double num = 9. + 3. * _rho0     + std::pow( _rho0    , 2 );
      double den = 9. + 3. * rho( m2 ) + std::pow( rho( m2 ), 2 );
      return std::sqrt( num / den );
    }

  // Maybe should throw an exception.

  return 0.;
}


double DalitzReso::angular( const double& mSq12, const double& mSq13, const double& mSq23 ) const
{
  if ( _L == 0 )
    return 1.;

  // Squared mass differences that some terms depend on.
  const double& diffSqM3 = _mSqM - _mSq3;
  const double& diffSq21 = _mSq2 - _mSq1;

  if ( _L == 1 )
    return mSq13 - mSq23 + diffSqM3 * diffSq21 / mSq12;

  if ( _L == 2 )
    {
      // Squared mass sums that some terms depend on.
      const double& sumSqM3 = _mSqM + _mSq3;
      const double& sumSq21 = _mSq2 + _mSq1;

      double first  = mSq13 - mSq23 + diffSqM3 * diffSq21 / mSq12;
      double second = mSq12 - 2. * sumSqM3 + std::pow( diffSqM3, 2 ) / mSq12;
      double third  = mSq12 - 2. * sumSq21 + std::pow( diffSq21, 2 ) / mSq12;

      return std::pow( first, 2 ) - second * third / 3.;
    }

  // Maybe should throw an exception.

  return 0.;
}


double DalitzReso::runningWidth( const double& m2 ) const
{
  double m = std::sqrt( m2 );

  double first  = std::pow( p( m2 ) / _p0, 2 * _L + 1 );
  double second = _m0 / m;
  double third  = std::pow( blattWeisskopf( m2 ), 2 );

  return _g0 * first * second * third;
}


double DalitzReso::mSquared12( const double& m2AB, const double& m2AC, const double& m2BC ) const
{
  // mSq12 is the pair of _pair1 and _pair2.
  if ( _other == 'C' )
    return m2AB;
  if ( _other == 'B' )
    return m2AC;
  if ( _other == 'A' )
    return m2BC;

  return 0.;
}


double DalitzReso::mSquared13( const double& m2AB, const double& m2AC, const double& m2BC ) const
{
  // mSq13 is the pair of _pair1 and _other.
  if ( _pair2 == 'C' )
    return m2AB;
  if ( _pair2 == 'B' )
    return m2AC;
  if ( _pair2 == 'A' )
    return m2BC;

  return 0.;
}


double DalitzReso::mSquared23( const double& m2AB, const double& m2AC, const double& m2BC ) const
{
  // mSq23 is the pair of _pair2 and _other.
  if ( _pair1 == 'C' )
    return m2AB;
  if ( _pair1 == 'B' )
    return m2AC;
  if ( _pair1 == 'A' )
    return m2BC;

  return 0.;
}



std::complex< double > DalitzReso::rbw( const double& m2AB, const double& m2AC, const double& m2BC ) const
{
  // Decide which of the arguments are mSq12, mSq13 and mSq23.
  double mSq12 = mSquared12( m2AB, m2AC, m2BC );
  double mSq13 = mSquared13( m2AB, m2AC, m2BC );
  double mSq23 = mSquared23( m2AB, m2AC, m2BC );

  std::complex< double > propagator = 1. / std::complex< double >( _mSq0 - mSq12, - _m0 * runningWidth( mSq12 ) );
  propagator *= blattWeisskopf( mSq12 ) * angular( mSq12, mSq13, mSq23 );

  return propagator;
}



std::complex< double > DalitzReso::gs( const double& m2AB, const double& m2AC, const double& m2BC ) const
{
  // Decide which of the arguments are mSq12, mSq13 and mSq23.
  double mSq12 = mSquared12( m2AB, m2AC, m2BC );
  double mSq13 = mSquared13( m2AB, m2AC, m2BC );
  double mSq23 = mSquared23( m2AB, m2AC, m2BC );

  double m12   = std::sqrt( mSq12 );

  std::complex< double > propagator = 1. / std::complex< double >( _mSq0 - mSq12 + gsf( mSq12 ), - m12 * runningWidth( mSq12 ) );
  propagator *= 1. + _gsd * _g0 / _m0;
  propagator *= blattWeisskopf( mSq12 ) * angular( mSq12, mSq13, mSq23 );

  return propagator;
}


double DalitzReso::gsf( const double& mSq12 ) const
{
  double factor = _g0 * _mSq0 / _p0;

  double first  = ( p2( mSq12 ) / _pSq0 ) * ( gsh( mSq12 ) - gsh( _mSq0 ) );
  double second = ( _mSq0 - mSq12 ) * gshprime( _mSq0 );

  return factor * ( first + second );
}


double DalitzReso::gsh( const double& mSq12 ) const
{
  double m12 = std::sqrt( mSq12 );
  double p12 = p        ( mSq12 );

  return ( 2. / M_PI ) * ( p12 / m12 ) * log( ( m12 + 2. * p12 ) / ( 2. * _m1 ) );
}

double DalitzReso::gshprime( const double& mSq12 ) const
{
  double first  = 1. / ( 8. * p2( mSq12 ) );
  double second = 1. / ( 2. * mSq12       );

  return ( first - second ) * gsh( mSq12 ) + second / M_PI;
}


std::complex< double > DalitzReso::lass( const double& m2AB, const double& m2AC, const double& m2BC ) const
{
  // Decide which of the arguments is mSq12.
  double mSq12 = mSquared12( m2AB, m2AC, m2BC );

  double p12 = p( mSq12 );

  double cotDeltaB = 1. / ( _lassa * p12 ) + _lassr * p12 / 2.;
  double deltaB    = atan ( 1. / cotDeltaB );
  double deltaR    = atan( ( _m0 * runningWidth( mSq12 ) ) / ( _mSq0 - mSq12 ) );

  double phaseB = deltaB + _lassPhiB;
  double phaseR = deltaR + _lassPhiR;

  std::complex< double > I( 0., 1. );

  std::complex< double > phaseTermB = exp( I * phaseB );
  std::complex< double > phaseTermR = exp( I * phaseR );

  std::complex< double > first  = _lassB * std::sin( phaseB );
  std::complex< double > second = _lassR * std::sin( deltaR ) * phaseTermR * phaseTermB;

  return phaseTermB * ( first + second );
}
