
#include <iostream>
#include <algorithm>
#include <functional>

#include <cfit/models/d0mix.hh>
#include <cfit/models/d0mix/dalitzReso.hh>


// Values implemented in BaBar.
const double& DalitzD0mix::_mPi = 0.139570; //.13957018;
const double& DalitzD0mix::_mKs = 0.49767;
const double& DalitzD0mix::_mD0 = 1.8645;


// Constructor.
DalitzD0mix::DalitzD0mix( const std::vector< Variable >& vars, const std::vector< Parameter >& pars )
{
  const int nReso = 15;

  _qp = std::complex< double >( 1., 0. );

  for( std::vector< Variable >::const_iterator var = vars.begin(); var != vars.end(); ++var )
    push( *var );

  for( std::vector< Parameter >::const_iterator par = pars.begin(); par != pars.end(); ++par )
    push( *par );

  _reso = new DalitzReso[ nReso ];

  // Defining resonances.
  _reso[  0 ] = DalitzReso( 'A', 'C', 1 ); // KStarm
  _reso[  1 ] = DalitzReso( 'A', 'B', 1 ); // KStarp
  _reso[  2 ] = DalitzReso( 'B', 'C', 1 ); // rho0
  _reso[  3 ] = DalitzReso( 'B', 'C', 1 ); // omega
  _reso[  4 ] = DalitzReso( 'B', 'C', 0 ); // f0_980
  _reso[  5 ] = DalitzReso( 'B', 'C', 0 ); // f0_1370
  _reso[  6 ] = DalitzReso( 'B', 'C', 2 ); // f2(1270)
  _reso[  7 ] = DalitzReso( 'A', 'C', 0 ); // K0Starm(1430)
  _reso[  8 ] = DalitzReso( 'A', 'B', 0 ); // K0Starp(1430)
  _reso[  9 ] = DalitzReso( 'A', 'C', 2 ); // K2Starm(1430)
  _reso[ 10 ] = DalitzReso( 'A', 'B', 2 ); // K2Starp(1430)
  _reso[ 11 ] = DalitzReso( 'B', 'C', 0 ); // sigma
  _reso[ 12 ] = DalitzReso( 'B', 'C', 0 ); // sigma2
  _reso[ 13 ] = DalitzReso( 'A', 'C', 1 ); // KStarm(1680)

  // Set the parameters of the resonances.
  _reso[  0 ].setReso( mKstarc      (), wKstarc      () );
  _reso[  1 ].setReso( mKstarc      (), wKstarc      () );
  _reso[  2 ].setReso( mRho         (), wRho         () );
  _reso[  3 ].setReso( mOmega       (), wOmega       () );
  _reso[  4 ].setReso( mF0_980      (), wF0_980      () );
  _reso[  5 ].setReso( mF0_1370     (), wF0_1370     () );
  _reso[  6 ].setReso( mF2_1270     (), wF2_1270     () );
  _reso[  7 ].setReso( mK0starc_1430(), wK0starc_1430() );
  _reso[  8 ].setReso( mK0starc_1430(), wK0starc_1430() );
  _reso[  9 ].setReso( mK2starc_1430(), wK2starc_1430() );
  _reso[ 10 ].setReso( mK2starc_1430(), wK2starc_1430() );
  _reso[ 11 ].setReso( mSigma       (), wSigma       () );
  _reso[ 12 ].setReso( mSigma2      (), wSigma2      () );
  _reso[ 13 ].setReso( mKstarm_1680 (), wKstarm_1680 () );

  // Check if all the propagators are fixed.
  _propagatorsFixed = true;
  for ( int par = 3; ( par < 25 ) && _propagatorsFixed; par++ )
    _propagatorsFixed &= getPar( par ).isFixed();

  //_propagatorsFixed = false; // DEBUG only.

  // Check if all amplitudes and propagators are fixed.
  _integralsFixed = _propagatorsFixed;
  for ( int par = 25; ( par < 55 ) && _integralsFixed; par++ )
    _integralsFixed &= getPar( par ).isFixed();

  // Decide which amplitudes will need to be recomputed at each minuit iteration.
  _unfixed = new bool[ nReso ];
  _unfixed[  0 ] = getPar(  3 ).isReleased() || getPar(  4 ).isReleased();   // KStarm
  _unfixed[  1 ] = _unfixed[ 0 ];                                            // KStarp
  _unfixed[  2 ] = getPar(  5 ).isReleased() || getPar(  6 ).isReleased();   // Rho
  _unfixed[  3 ] = getPar(  7 ).isReleased() || getPar(  8 ).isReleased();   // Omega
  _unfixed[  4 ] = getPar(  9 ).isReleased() || getPar( 10 ).isReleased();   // F0_980
  _unfixed[  5 ] = getPar( 11 ).isReleased() || getPar( 12 ).isReleased();   // F0_1370
  _unfixed[  6 ] = getPar( 13 ).isReleased() || getPar( 14 ).isReleased();   // F2_1270
  _unfixed[  7 ] = getPar( 15 ).isReleased() || getPar( 16 ).isReleased();   // K0Starm
  _unfixed[  8 ] = _unfixed[ 7 ];                                            // K0Starp
  _unfixed[  9 ] = getPar( 17 ).isReleased() || getPar( 18 ).isReleased();   // K2Starm
  _unfixed[ 10 ] = _unfixed[ 9 ];                                            // K2Starp
  _unfixed[ 11 ] = getPar( 19 ).isReleased() || getPar( 20 ).isReleased();   // sigma
  _unfixed[ 12 ] = getPar( 21 ).isReleased() || getPar( 22 ).isReleased();   // sigma2
  _unfixed[ 13 ] = getPar( 23 ).isReleased() || getPar( 24 ).isReleased();   // KStarm1680
  _unfixed[ 14 ] = false;                                                    // Non-resonant

  // Decide which amplitudes are the conjugate of their previous ones to avoid recomputing.
  _conjOfPrevious = new bool[ nReso ];
  _conjOfPrevious[  0 ] = false; // KStarm
  _conjOfPrevious[  1 ] = true;  // KStarp
  _conjOfPrevious[  2 ] = false; // Rho
  _conjOfPrevious[  3 ] = false; // Omega
  _conjOfPrevious[  4 ] = false; // f0_980
  _conjOfPrevious[  5 ] = false; // f0_1370
  _conjOfPrevious[  6 ] = false; // f2_1270
  _conjOfPrevious[  7 ] = false; // K0Starm
  _conjOfPrevious[  8 ] = true;  // K0Starp
  _conjOfPrevious[  9 ] = false; // K2Starm
  _conjOfPrevious[ 10 ] = true;  // K2Starp
  _conjOfPrevious[ 11 ] = false; // sigma
  _conjOfPrevious[ 12 ] = false; // sigma2
  _conjOfPrevious[ 13 ] = false; // KStarm1680
  _conjOfPrevious[ 14 ] = false; // Non-resonant

  /*
  for ( int i = 0; i < 10; i++ ) // DEBUG only.
    {
      _unfixed       [ i ] = true;
      _conjOfPrevious[ i ] = false;
    }
  */

  // Allocate memory for the direct and crossed matrices.
  _Id = new std::complex< double >*[ nReso ];
  _Ix = new std::complex< double >*[ nReso ];
  for ( int reso = 0; reso < nReso; reso++ )
    {
      _Id[ reso ] = new std::complex< double >[ nReso ];
      _Ix[ reso ] = new std::complex< double >[ nReso ];
    }

  // If all the resonances are fixed, compute the matrices of integrals just
  //    once here in the constructor.
  if ( _propagatorsFixed )
    cacheIntegralsMatrix( true );

  // If all the resonances and amplitudes are fixed, compute the I_1 and I_\chi
  //    integrals just once here in the constructor.
  if ( _integralsFixed )
    cacheDalitzIntegrals();
}


void DalitzD0mix::setPars( const std::vector< double >& pars ) throw( PdfException )
{
  PdfModel::setPars( pars );

  _reso[  0 ].setReso( mKstarc      (), wKstarc      () );
  _reso[  1 ].setReso( mKstarc      (), wKstarc      () );
  _reso[  2 ].setReso( mRho         (), wRho         () );
  _reso[  3 ].setReso( mOmega       (), wOmega       () );
  _reso[  4 ].setReso( mF0_980      (), wF0_980      () );
  _reso[  5 ].setReso( mF0_1370     (), wF0_1370     () );
  _reso[  6 ].setReso( mF2_1270     (), wF2_1270     () );
  _reso[  7 ].setReso( mK0starc_1430(), wK0starc_1430() );
  _reso[  8 ].setReso( mK0starc_1430(), wK0starc_1430() );
  _reso[  9 ].setReso( mK2starc_1430(), wK2starc_1430() );
  _reso[ 10 ].setReso( mK2starc_1430(), wK2starc_1430() );
  _reso[ 11 ].setReso( mSigma       (), wSigma       () );
  _reso[ 12 ].setReso( mSigma2      (), wSigma2      () );
  _reso[ 13 ].setReso( mKstarm_1680 (), wKstarm_1680 () );

  return;
}



double DalitzD0mix::evaluate() const throw( PdfException )
{
//   if ( ( m2AC() < m2ACmin( m2AB() ) ) || ( m2AC() > m2ACmax( m2AB() ) ) )
//     return 0.;

  // Dalitz amplitudes of the decay of the particle and that of the antiparticle.
  std::complex< double > ampDalitz     = dalitzKsPiPi( m2AB(), m2AC(), m2BC() );
  std::complex< double > ampAntiDalitz = dalitzKsPiPi( m2AC(), m2AB(), m2BC() );

  // Assume there's no direct CP violation.
  std::complex< double > barAOverA = ampAntiDalitz / ampDalitz;

  // CP violation in the interference. _qp implements CP violation in the mixing.
  std::complex< double > chi = _qp * barAOverA;

  // Gamma t. Notice that Gamma is not 1 / tau, since it also depends on x and y,
  //    as well as on integrals of products of amplitudes over the Dalitz plot.
  double gt = t() / invGamma();

  // Compute time dependent amplitude.
  std::complex< double > amp = .5 * ampDalitz * ( ( 1. + chi ) * e1( gt ) + ( 1. - chi ) * e2( gt ) );

  // std::norm returns the squared modulus of the complex number, not its norm.
  return std::norm( amp ) / _norm;
}



std::complex< double > DalitzD0mix::evaluateReso( int index, double m2AB, double m2AC, double m2BC ) const
{
  if ( index ==  0 ) return _reso[  0 ].rbw ( m2AB, m2AC, m2BC );
  if ( index ==  1 ) return _reso[  1 ].rbw ( m2AB, m2AC, m2BC );
  if ( index ==  2 ) return _reso[  2 ].gs  ( m2AB, m2AC, m2BC );
  if ( index ==  3 ) return _reso[  3 ].rbw ( m2AB, m2AC, m2BC );
  if ( index ==  4 ) return _reso[  4 ].rbw ( m2AB, m2AC, m2BC );
  if ( index ==  5 ) return _reso[  5 ].rbw ( m2AB, m2AC, m2BC );
  if ( index ==  6 ) return _reso[  6 ].rbw ( m2AB, m2AC, m2BC );
  if ( index ==  7 ) return _reso[  7 ].rbw ( m2AB, m2AC, m2BC );
  if ( index ==  8 ) return _reso[  8 ].rbw ( m2AB, m2AC, m2BC );
  if ( index ==  9 ) return _reso[  9 ].rbw ( m2AB, m2AC, m2BC );
  if ( index == 10 ) return _reso[ 10 ].rbw ( m2AB, m2AC, m2BC );
  if ( index == 11 ) return _reso[ 11 ].rbw ( m2AB, m2AC, m2BC );
  if ( index == 12 ) return _reso[ 12 ].rbw ( m2AB, m2AC, m2BC );
  if ( index == 13 ) return _reso[ 13 ].rbw ( m2AB, m2AC, m2BC );
  if ( index == 14 ) return std::complex< double >( 1., 0. ); // Non-resonant term.

  return std::complex< double >();
}


std::complex< double > DalitzD0mix::coef( int index ) const
{
  if ( index ==  0 ) return std::complex< double >( rKstarm      (), iKstarm      () );
  if ( index ==  1 ) return std::complex< double >( rKstarp      (), iKstarp      () );
  if ( index ==  2 ) return std::complex< double >( rRho         (), iRho         () );
  if ( index ==  3 ) return std::complex< double >( rOmega       (), iOmega       () );
  if ( index ==  4 ) return std::complex< double >( rF0_980      (), iF0_980      () );
  if ( index ==  5 ) return std::complex< double >( rF0_1370     (), iF0_1370     () );
  if ( index ==  6 ) return std::complex< double >( rF2_1270     (), iF2_1270     () );
  if ( index ==  7 ) return std::complex< double >( rK0starm_1430(), iK0starm_1430() );
  if ( index ==  8 ) return std::complex< double >( rK0starp_1430(), iK0starp_1430() );
  if ( index ==  9 ) return std::complex< double >( rK2starm_1430(), iK2starm_1430() );
  if ( index == 10 ) return std::complex< double >( rK2starp_1430(), iK2starp_1430() );
  if ( index == 11 ) return std::complex< double >( rSigma       (), iSigma       () );
  if ( index == 12 ) return std::complex< double >( rSigma2      (), iSigma2      () );
  if ( index == 13 ) return std::complex< double >( rKstarm_1680 (), iKstarm_1680 () );
  if ( index == 14 ) return std::complex< double >( rNonReson    (), iNonReson    () );

  return std::complex< double >();
}



double DalitzD0mix::norm() const
{
  // Add time dependent terms to the norm.
  double x2 = std::pow( x(), 2 );
  double y2 = std::pow( y(), 2 );

  double first  = ( 1. + std::norm( _qp ) ) / 2. - y() * real( _iChi );
  double second = ( 1. - std::norm( _qp ) ) / 2. + x() * imag( _iChi );

  return real( _i1 ) * invGamma() * ( first / ( 1. - y2 ) + second / ( 1. + x2 ) );
}



// If requested, compute all the terms, regardless of being changed or not.
void DalitzD0mix::cacheIntegralsMatrix( bool all ) const
{
  // Define the amplitudes that are gonna be computed at each point.
  std::complex< double > A   [ 15 ];
  std::complex< double > ABar[ 15 ];

  // Define the properties of the integration method.
  int    nBins = 400;
  double min   = pow( mKs() + mPi(), 2 );
  double max   = pow( mD0() - mPi(), 2 );
  double step  = ( max - min ) / double( nBins );

  // Define the variables at each bin.
  double m2AB;
  double m2AC;
  double m2BC;

  // Initialize the matrix elements that need to be recomputed.
  for ( int i = 0; i < 15; i++ )
    for ( int j = i; j < 15; j++ )
      if ( all || _unfixed[ i ] || _unfixed[ j ] )
	{
	  _Id[ i ][ j ] = std::complex< double >();
	  _Ix[ i ][ j ] = std::complex< double >();
	}

  // Compute the integral on the grid.
  for ( int binX = 0; binX < nBins; binX++ )
    for ( int binY = binX; binY < nBins; binY++ )
      {
	m2AB = min + step * ( binX + .5 );
	m2AC = min + step * ( binY + .5 );
	m2BC = pow( mD0(), 2 ) + pow( mKs(), 2 ) + 2. * pow( mPi(), 2 ) - m2AB - m2AC;

	// Proceed only if the point lies inside the kinematically allowed Dalitz region.
	if ( ( m2AC > m2ACmin( m2AB ) ) && ( m2AC < m2ACmax( m2AB ) ) )
	  {
	    // Compute the resonances that may have varied.
	    for ( int reso = 0; reso < 15; reso++ )
	      if ( _conjOfPrevious[ reso ] )
		{
		  A   [ reso ] = ABar[ reso - 1 ];
		  ABar[ reso ] = A   [ reso - 1 ];
		}
	      else
		{
		  A   [ reso ] = evaluateReso( reso, m2AB, m2AC, m2BC );
		  ABar[ reso ] = evaluateReso( reso, m2AC, m2AB, m2BC );
		}

	    // Add the terms to the matrices.
	    for ( int i = 0; i < 15; i++ )
	      for ( int j = i; j < 15; j++ )
		if ( all || _unfixed[ i ] || _unfixed[ j ] )
		  {
		    if ( binX == binY )
		      {
			_Id[ i ][ j ] += ( conj( A[ i ] ) * A   [ j ] + conj( ABar[ i ] ) * ABar[ j ] ) / 2.;
			_Ix[ i ][ j ] += ( conj( A[ i ] ) * ABar[ j ] + conj( ABar[ i ] ) * A   [ j ] ) / 2.;
		      }
		    else
		      {
			_Id[ i ][ j ] +=   conj( A[ i ] ) * A   [ j ] + conj( ABar[ i ] ) * ABar[ j ];
			_Ix[ i ][ j ] +=   conj( A[ i ] ) * ABar[ j ] + conj( ABar[ i ] ) * A   [ j ];
		      }
		  }
	  }
      }

  // Finish the computation of the integrals and compute the lower triangle elements.
  for ( int i = 0; i < 15; i++ )
    for ( int j = i; j < 15; j++ )
      {
	_Id[ i ][ j ] *= pow( step, 2 );
	_Ix[ i ][ j ] *= pow( step, 2 );

	if ( i != j )
	  {
	    _Id[ j ][ i ] = conj( _Id[ i ][ j ] );
	    _Ix[ j ][ i ] = conj( _Ix[ i ][ j ] );
	  }
      }

  return;
}


// Compute _i1 and _iChi after the matrix elements have been cached.
void DalitzD0mix::cacheDalitzIntegrals()
{
  _i1   = std::complex< double >();
  _iChi = std::complex< double >();

  for ( int i = 0; i < 15; i++ )
    for ( int j = 0; j < 15; j++ )
      {
	_i1   += conj( coef( i ) ) * _Id[ i ][ j ] * coef( j );
	_iChi += conj( coef( i ) ) * _Ix[ i ][ j ] * coef( j );
      }

  _iChi *= _qp / _i1;

  return;
}



void DalitzD0mix::cache()
{
  // Compute the terms of the matrices of integrals (Id and Ix) that need recalculation.
  if ( ! _propagatorsFixed )
    cacheIntegralsMatrix();

  // Compute the integrals I_1 and I_\chi if the amplitudes are constant.
  if ( ! _integralsFixed )
    cacheDalitzIntegrals();

  _norm = norm();
}



double DalitzD0mix::kallen( double x, double y, double z )
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


double DalitzD0mix::m2ACmin( double m2AB )
{
  double first  = pow( pow( _mD0, 2 ) + pow( _mKs, 2 ) - 2. * pow( _mPi, 2 ), 2 );
  double second = sqrt( kallen( m2AB, pow( _mKs, 2 ), pow( _mPi, 2 ) ) );
  double third  = sqrt( kallen( m2AB, pow( _mD0, 2 ), pow( _mPi, 2 ) ) );

  return ( first - pow( second + third, 2 ) ) / ( 4. * m2AB );
}


double DalitzD0mix::m2ACmax( double m2AB )
{
  double first  = pow( pow( _mD0, 2 ) + pow( _mKs, 2 ) - 2. * pow( _mPi, 2 ), 2 );
  double second = sqrt( kallen( m2AB, pow( _mKs, 2 ), pow( _mPi, 2 ) ) );
  double third  = sqrt( kallen( m2AB, pow( _mD0, 2 ), pow( _mPi, 2 ) ) );

  return ( first - pow( second - third, 2 ) ) / ( 4. * m2AB );
}


std::complex< double > DalitzD0mix::dalitzKsPiPi( double m2AB, double m2AC, double m2BC ) const
{
  // Adding terms to the amplitude with their corresponding amplitude and phase terms.
  std::complex< double > amp( 0., 0. );
  for ( int reso = 0; reso < 15; reso++ )
    amp += coef( reso ) * evaluateReso( reso, m2AB, m2AC, m2BC );

  return amp;
}

// < f | H | D^0 (t) > = 1/2 * [ ( 1 + \chi_f ) * A_f * e_1(gt) + ( 1 - \chi_f ) * A_f * e_2(gt) ]
std::complex< double > DalitzD0mix::e1( const double& gt ) const
{
  double real = exp ( - ( 1. + y() ) * gt / 2. ) * cos(   x() * gt / 2. );
  double imag = exp ( - ( 1. + y() ) * gt / 2. ) * sin( - x() * gt / 2. );

  return std::complex< double >( real, imag );
}

std::complex< double > DalitzD0mix::e2( const double& gt ) const
{
  double real = exp ( - ( 1. - y() ) * gt / 2. ) * cos( x() * gt / 2. );
  double imag = exp ( - ( 1. - y() ) * gt / 2. ) * sin( x() * gt / 2. );

  return std::complex< double >( real, imag );
}

double DalitzD0mix::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  throw PdfException( "Do not use this function yet." );
  return 0.;
}
