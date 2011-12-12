
#include <cfit/fvector.hh>
#include <cfit/phasespace.hh>


const double& Fvector::_mPi = 0.139570;



void Fvector::pushBeta( const std::vector< Coef >& beta )
{
  typedef std::vector< Coef >::const_iterator kIter;

  for ( kIter c = beta.begin(); c != beta.end(); ++c )
  {
    _parMap[ c->real().name() ] = c->real();
    _parMap[ c->imag().name() ] = c->imag();
    _beta.push_back( std::make_pair( c->real().name(), c->imag().name() ) );
  }
}


void Fvector::pushfPr( const std::vector< Coef >& fPr  )
{
  typedef std::vector< Coef >::const_iterator kIter;

  for ( kIter p = fPr.begin(); p != fPr.end(); ++p )
  {
    _parMap[ p->real().name() ] = p->real();
    _parMap[ p->imag().name() ] = p->imag();
    _fPr.push_back( std::make_pair( p->real().name(), p->imag().name() ) );
  }
}



void Fvector::pushS0pr( const Parameter& par )
{
  _parMap[ par.name() ] = par;
  _s0pr = par.name();
}


// void Fvector::push( const Coef& coef )
// {
//   _parMap[ coef.real().name() ] = coef.real();
//   _parMap[ coef.imag().name() ] = coef.imag();

//   _parOrder.push_back( coef.real().name() );
//   _parOrder.push_back( coef.imag().name() );
// }


// // For resonances with larger number of parameters, be able to get them by index.
// //    Important: the zeroth extra parameter is the 3rd element in the vector.
// double Fvector::getPar( const unsigned index ) const throw( PdfException )
// {
//   if ( _parOrder.size() > index + 3 )
//     return _parMap.find( _parOrder[ index + 3 ] )->second.value();

//   throw PdfException( "Trying to access unexisting parameter." );
// }


void Fvector::setPars( const std::map< std::string, Parameter >& pars )
{
  typedef std::map< const std::string, Parameter >::iterator pIter;
  for ( pIter par = _parMap.begin(); par != _parMap.end(); ++par )
    par->second.setValue( pars.find( par->first )->second.value() );
}

// Kallen function lambda( x, y, z ) = x^2 + y^2 + z^2 - 2xy - 2xz - 2yz.
double Fvector::kallen( const double& x, const double& y, const double& z )
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

// Invariant mass of the resonant pair.
double Fvector::m2AB( const double& mSq12, const double& mSq13, const double& mSq23 ) const
{
  if ( _noRes == 3 ) return mSq12;
  if ( _noRes == 2 ) return mSq13;
  if ( _noRes == 1 ) return mSq23;

  return 0.;
}

// Invariant mass of the first non-resonant pair.
double Fvector::m2AC( const double& mSq12, const double& mSq13, const double& mSq23 ) const
{
  if ( _resoB == 3 ) return mSq12;
  if ( _resoB == 2 ) return mSq13;
  if ( _resoB == 1 ) return mSq23;

  return 0.;
}

// Invariant mass of the second non-resonant pair.
double Fvector::m2BC( const double& mSq12, const double& mSq13, const double& mSq23 ) const
{
  if ( _resoA == 3 ) return mSq12;
  if ( _resoA == 2 ) return mSq13;
  if ( _resoA == 1 ) return mSq23;

  return 0.;
}

// Momentum of a resonant particle in the rest frame of the resonant pair.
double Fvector::q( const PhaseSpace& ps, const double& mSqAB ) const
{
  return std::sqrt( kallen( mSqAB, ps.mSq( _resoA ), ps.mSq( _resoB ) ) ) / ( 2. * std::sqrt( mSqAB ) );
}

// Squared momentum of a resonant particle in the rest frame of the resonant pair.
double Fvector::qSq( const PhaseSpace& ps, const double& mSqAB ) const
{
  return kallen( mSqAB, ps.mSq( _resoA ), ps.mSq( _resoB ) ) / ( 4. * mSqAB );
}


// Phase space factor, 2 q / m, where q is the momentum of a resonant particle in the
//    rest frame of the resonant pair, and m is the invariant mass of the resonant pair.
double Fvector::rho( const PhaseSpace& ps, const double& mSqAB ) const
{
  return std::sqrt( kallen( mSqAB, ps.mSq( _resoA ), ps.mSq( _resoB ) ) ) / mSqAB;
}


std::complex< double > Fvector::rho( const double& mCh, const double& mSqAB ) const
{
  double rhoSq = 1. - std::pow( mCh, 2 ) / mSqAB;

  if ( rhoSq >= 0. )
    return std::sqrt( rhoSq );

  const std::complex< double > I( 0., 1. );

  return I * std::sqrt( - rhoSq );
}


std::complex< double > Fvector::rho4pi( const double& mSqAB ) const
{
  const double& mPi = 0.139570;

  if ( mSqAB > 1. )
    return rho( 4. * mPi, mSqAB );

  const double& m4AB = std::pow( mSqAB, 2 );
  const double& m6AB = std::pow( mSqAB, 3 );
  const double& m8AB = std::pow( mSqAB, 4 );

  double term = 0.;
  term += 0.00370909 / m4AB;
  term -= 0.111203   / mSqAB;
  term += 1.2274;
  term -= 6.39017    * mSqAB;
  term += 16.8358    * m4AB;
  term -= 21.8845    * m6AB;
  term += 11.3153    * m8AB;

  return rho( 4. * mPi, 1. ) * term;
}


std::complex< double > Fvector::rho( const int& index, const double& mSqAB ) const
{
  const double& mPi   = 0.139570;
  const double& mK    = 0.49368; // Charged K mass.
  const double& mEta  = 0.54730;
  const double& mEtaP = 0.95777;

  if ( index == 0 ) return rho   ( 2. * mPi    , mSqAB );
  if ( index == 1 ) return rho   ( 2. * mK     , mSqAB );
  if ( index == 2 ) return rho4pi(               mSqAB );
  if ( index == 3 ) return rho   ( 2. * mEta   , mSqAB );
  if ( index == 4 ) return rho   ( mEta + mEtaP, mSqAB );

  return std::complex< double >( 1., 0. );
}


std::complex< double > Fvector::propagator( const PhaseSpace& ps, const double& mSqAB ) const
{
  Matrix< double > K( 5 );

  // Initialize the K matrix to zero.
  for ( int row = 0; row < 5; ++row )
    for ( int col = 0; col < 5; ++col )
      K( row, col ) = 0.;

  // Resonant contribution.
  for ( int row = 0; row < 5; ++row )
    for ( int col = 0; col < 5; ++col )
      for ( int pole = 0; pole < 5; ++pole )
        K( row, col ) += _g0( pole, row ) * _g0( pole, col ) / ( std::pow( _m0[ pole ], 2 ) - mSqAB );

  // Non-resonant contribution.
  K += _fSc * ( 1. - _s0sc ) / ( mSqAB - _s0sc );

  // Adler term.
  K *= ( 1. - _s0A ) / ( mSqAB - _s0A ) * ( mSqAB - _sA * std::pow( _mPi, 2 ) / 2. );

  Matrix< std::complex< double > > M( 5 );
  const std::complex< double > I( 0., 1. );

  // Build M = ( 1 - i K rho ).
  for ( int row = 0; row < 5; ++row )
    for ( int col = 0; col < 5; ++col )
      M( row, col ) = 1. * ( row == col ) - I * K( row, col ) * rho( col, mSqAB );

  // Invert M: ( 1 - i K rho )^{ -1 }
  Matrix< std::complex< double > > invM;
  invM = M.inverse();

  // Build the P vector.
  double svp = ( 1. - getPar( _s0pr ) ) / ( mSqAB - getPar( _s0pr ) ); // Slowly varying part.
  std::vector< std::complex< double > > P;
  for ( int row = 0; row < 5; ++row )
  {
    std::complex< double > poles = 0.;
    for ( int pole = 0; pole < 5; ++pole )
      poles += getCoef( _beta[ pole ] ) * _g0( pole, row ) / ( std::pow( _m0[ pole ], 2 ) - mSqAB );

    P.push_back( poles + getCoef( _fPr[ row ] ) * svp );
  }

  std::complex< double > result = 0.;

  // Compute the first component of the F vector: F_0 = ( 1 - i K rho )_{0j}^{ -1 } P_j.
  for ( int line = 0; line < 5; ++line )
    result += invM( 0, line ) * P[ line ];

  return result;
}


std::complex< double > Fvector::evaluate( const PhaseSpace& ps,
                                          const double&     mSq12,
                                          const double&     mSq13,
                                          const double&     mSq23 ) const
{
  // Determine the resonant pair.
  const double& mSqAB = m2AB( mSq12, mSq13, mSq23 );

  return propagator( ps, mSqAB );
}
