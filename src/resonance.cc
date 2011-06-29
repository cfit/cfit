
#include <cfit/resonance.hh>
#include <cfit/phasespace.hh>

double Resonance::kallen( const double& x, const double& y, const double& z )
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

double Resonance::m2AB( const double& mSq12, const double& mSq13, const double& mSq23 ) const
{
  if ( _noRes == 3 ) return mSq12;
  if ( _noRes == 2 ) return mSq13;
  if ( _noRes == 1 ) return mSq23;

  return 0.;
}

double Resonance::m2AC( const double& mSq12, const double& mSq13, const double& mSq23 ) const
{
  if ( _resoB == 3 ) return mSq12;
  if ( _resoB == 2 ) return mSq13;
  if ( _resoB == 1 ) return mSq23;

  return 0.;
}

double Resonance::m2BC( const double& mSq12, const double& mSq13, const double& mSq23 ) const
{
  if ( _resoA == 3 ) return mSq12;
  if ( _resoA == 2 ) return mSq13;
  if ( _resoA == 1 ) return mSq23;

  return 0.;
}


double Resonance::q( const PhaseSpace& ps, const double& mSqAB ) const
{
  return std::sqrt( kallen( mSqAB, ps.mSq( _resoA ), ps.mSq( _resoB ) ) ) / ( 2. * std::sqrt( mSqAB ) );
}


// Phase space factor, 2 q / m, where q is the momentum of a resonant particle in the
//    rest frame of the resonant pair, and m is the invariant mass of the resonant pair.
double Resonance::rho( const PhaseSpace& ps, const double& mSqAB ) const
{
  return std::sqrt( kallen( mSqAB, ps.mSq( _resoA ), ps.mSq( _resoB ) ) ) / mSqAB;
}



double Resonance::runningWidth( const PhaseSpace& ps, const double& mSqAB ) const
{
  const double& rho0 = rho( ps, std::pow( _mass.value(), 2 ) );
  return _width.value() * ( rho( ps, mSqAB ) / rho0 ) * std::pow( blattWeisskopf( ps, mSqAB ), 2 );
}



double Resonance::zemach( const PhaseSpace& ps, const double& mSqAB, const double& mSqAC, const double& mSqBC ) const
{
//   std::cout << "DEBUG zemach: " << _l << " " << _mMotherSq << " " << _mNonResoSq << " " << _mResoSqA << " " << _mResoSqB << std::endl;

  if ( _l == 0 )
    return 1.;

  // Squared mass differences that some terms depend on.
  const double& diffSqMC = ps.mSqMother()   - ps.mSq( _noRes );
  const double& diffSqAB = ps.mSq( _resoA ) - ps.mSq( _resoB );

  // Zemach tensor for l = 1.
  const double& zemach1  = mSqAC - mSqBC - diffSqMC * diffSqAB / mSqAB;

  if ( _l == 1 )
    return zemach1;

  if ( _l == 2 )
    {
      // Squared mass sums that some terms depend on.
      const double& sumSqMC = ps.mSqMother()   + ps.mSq( _noRes );
      const double& sumSqAB = ps.mSq( _resoA ) + ps.mSq( _resoB );

      double first  = mSqAB - 2. * sumSqMC + std::pow( diffSqMC, 2 ) / mSqAB;
      double second = mSqAB - 2. * sumSqAB + std::pow( diffSqAB, 2 ) / mSqAB;

      return std::pow( zemach1, 2 ) - first * second / 3.;
    }

  // Maybe should throw an exception.

  return 0.;
}



double Resonance::blattWeisskopfPrime( const PhaseSpace& ps, const double& mSqAB ) const
{
  if ( _l == 0 )
    return 1.;

  const double& q0    = q( ps, std::pow( _mass.value(), 2 ) );
  const double& qm    = q( ps, mSqAB );
  const double& rqSq0 = std::pow( _R.value() * q0, 2 );
  const double& rqSq  = std::pow( _R.value() * qm, 2 );

  if ( _l == 1 )
    return std::sqrt( ( 1. + rqSq0 ) / ( 1. + rqSq ) );

  if ( _l == 2 )
    {
      double num = 9. + 3. * rqSq0 + std::pow( rqSq0, 2 );
      double den = 9. + 3. * rqSq  + std::pow( rqSq , 2 );
      return std::sqrt( num / den );
    }

  // Maybe should throw an exception.

  return 0.;
}


// Blatt-Weisskopf centrifugal barrier factors.
// Maybe should throw an exception.
double Resonance::blattWeisskopf( const PhaseSpace& ps, const double& mSqAB ) const
{
  if ( _l == 0 )
    return 1.;

  const double& q0 = q( ps, std::pow( _mass.value(), 2 ) );
  const double& qm = q( ps, mSqAB );

  return std::pow( qm / q0, _l ) * blattWeisskopfPrime( ps, mSqAB );
}


std::complex< double > Resonance::evaluate( const PhaseSpace& ps, 
					    const double&     mSq12,
					    const double&     mSq13,
					    const double&     mSq23 ) const
{
  // Determine the resonant pair.
  const double& mSqAB = m2AB( mSq12, mSq13, mSq23 );

  // Determine the other pairs.
  const double& mSqAC = m2AC( mSq12, mSq13, mSq23 );
  const double& mSqBC = m2BC( mSq12, mSq13, mSq23 );

  return propagator( ps, mSqAB ) * zemach( ps, mSqAB, mSqAC, mSqBC ) * blattWeisskopf( ps, mSqAB );
}
