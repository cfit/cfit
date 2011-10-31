
#include <cmath>

#include <cfit/phasespace.hh>

const double PhaseSpace::kallen( const double& x, const double& y, const double& z )
{
  double result = 0.;
  result += pow( x, 2 );
  result += pow( y, 2 );
  result += pow( z, 2 );
  result -= 2. * x * y;
  result -= 2. * x * z;
  result -= 2. * y * z;

  return result;
}


const double PhaseSpace::mSq13min( const double& mSq12 ) const
{
  double first  = pow( pow( _mMother, 2 ) + pow( _m1, 2 ) - pow( _m2, 2 ) - pow( _m3, 2 ), 2 );
  double second = sqrt( kallen( mSq12, pow( _m1     , 2 ), pow( _m2, 2 ) ) );
  double third  = sqrt( kallen( mSq12, pow( _mMother, 2 ), pow( _m3, 2 ) ) );

  return ( first - pow( second + third, 2 ) ) / ( 4. * mSq12 );
}


const double PhaseSpace::mSq13max( const double& mSq12 ) const
{
  double first  = pow( pow( _mMother, 2 ) + pow( _m1, 2 ) - pow( _m2, 2 ) - pow( _m3, 2 ), 2 );
  double second = sqrt( kallen( mSq12, pow( _m1     , 2 ), pow( _m2, 2 ) ) );
  double third  = sqrt( kallen( mSq12, pow( _mMother, 2 ), pow( _m3, 2 ) ) );

  return ( first - pow( second - third, 2 ) ) / ( 4. * mSq12 );
}


const double PhaseSpace::mSq23min( const double& mSq12 ) const
{
  double first  = pow( pow( _mMother, 2 ) - pow( _m1, 2 ) + pow( _m2, 2 ) - pow( _m3, 2 ), 2 );
  double second = sqrt( kallen( mSq12, pow( _m1     , 2 ), pow( _m2, 2 ) ) );
  double third  = sqrt( kallen( mSq12, pow( _mMother, 2 ), pow( _m3, 2 ) ) );

  return ( first - pow( second + third, 2 ) ) / ( 4. * mSq12 );
}


const double PhaseSpace::mSq23max( const double& mSq12 ) const
{
  double first  = pow( pow( _mMother, 2 ) - pow( _m1, 2 ) + pow( _m2, 2 ) - pow( _m3, 2 ), 2 );
  double second = sqrt( kallen( mSq12, pow( _m1     , 2 ), pow( _m2, 2 ) ) );
  double third  = sqrt( kallen( mSq12, pow( _mMother, 2 ), pow( _m3, 2 ) ) );

  return ( first - pow( second - third, 2 ) ) / ( 4. * mSq12 );
}


const double PhaseSpace::mSq( unsigned index ) const
{
  if ( index == 0 )
    return _mSqMother;
  if ( index == 1 )
    return _mSq1;
  if ( index == 2 )
    return _mSq2;
  if ( index == 3 )
    return _mSq3;

  return 0.;
}


// Check if the kinematically allowed region contains a given point.
bool PhaseSpace::contains( const double& mSq12, const double& mSq13, const double& mSq23 ) const
{
  return ( mSq13 > mSq13min( mSq12 ) ) && ( mSq13 < mSq13max( mSq12 ) ) &&
         ( mSq23 > mSq23min( mSq12 ) ) && ( mSq23 < mSq23max( mSq12 ) );
}

