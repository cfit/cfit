
#include <algorithm>
#include <functional>

#include <cfit/models/d0mix.hh>
#include <cfit/models/d0mix/dalitzReso.hh>

// Initialize the static const values.
const double& DalitzD0mix::_mPi = .139570; //.13957018;
const double& DalitzD0mix::_mKs = .49767;  //.497648;
const double& DalitzD0mix::_mD0 = 1.8645;




// Constructor.
DalitzD0mix::DalitzD0mix( const std::vector< Variable >& vars, const std::vector< Parameter >& pars )
{
  const int nReso = 10;

  _qp = std::complex< double >( 1., 0. );

  for( std::vector< Variable >::const_iterator var = vars.begin(); var != vars.end(); ++var )
    push( *var );

  for( std::vector< Parameter >::const_iterator par = pars.begin(); par != pars.end(); ++par )
    push( *par );

  _reso = new DalitzReso[ nReso ];

  // Defining resonances.
  _reso[ 0 ] = DalitzReso( 'A', 'C', 1 ); // KStarm
  _reso[ 1 ] = DalitzReso( 'A', 'B', 1 ); // KStarp
  _reso[ 2 ] = DalitzReso( 'B', 'C', 1 ); // rho0
  _reso[ 3 ] = DalitzReso( 'B', 'C', 1 ); // omega
  _reso[ 4 ] = DalitzReso( 'B', 'C', 2 ); // f2(1270)
  _reso[ 5 ] = DalitzReso( 'A', 'C', 0 ); // K0Starm(1430)
  _reso[ 6 ] = DalitzReso( 'A', 'B', 0 ); // K0Starp(1430)
  _reso[ 7 ] = DalitzReso( 'A', 'C', 2 ); // K2Starm(1430)
  _reso[ 8 ] = DalitzReso( 'A', 'B', 2 ); // K2Starp(1430)
  _reso[ 9 ] = DalitzReso( 'A', 'C', 1 ); // KStarm(1680)

  // Set the parameters of the resonances.
  _reso[ 0 ].setReso( mKStarc()    , wKStarc()     );
  _reso[ 1 ].setReso( mKStarc()    , wKStarc()     );
  _reso[ 2 ].setReso( mRho()       , wRho()        );
  _reso[ 3 ].setReso( mOmega()     , wOmega()      );
  _reso[ 4 ].setReso( mF2()        , wF2()         );
  _reso[ 5 ].setReso( mK0Starc()   , wK0Starc()    );
  _reso[ 6 ].setReso( mK0Starc()   , wK0Starc()    );
  _reso[ 7 ].setReso( mK2Starc()   , wK2Starc()    );
  _reso[ 8 ].setReso( mK2Starc()   , wK2Starc()    );
  _reso[ 9 ].setReso( mKStarm1680(), wKStarm1680() );

  _reso[ 5 ].setLassa   ( aLass()  );
  _reso[ 5 ].setLassr   ( rLass()  );
  _reso[ 5 ].setLassB   ( bLass()  );
  _reso[ 5 ].setLassPhiB( pBLass() );
  _reso[ 5 ].setLassR   ( RLass()  );
  _reso[ 5 ].setLassPhiR( prLass() );

  _reso[ 6 ].setLassa   ( aLass()  );
  _reso[ 6 ].setLassr   ( rLass()  );
  _reso[ 6 ].setLassB   ( bLass()  );
  _reso[ 6 ].setLassPhiB( pBLass() );
  _reso[ 6 ].setLassR   ( RLass()  );
  _reso[ 6 ].setLassPhiR( prLass() );

  // Check if all the propagators are fixed.
  _propagatorsFixed = true;
  for ( int par = 3; ( par < 23 ) && _propagatorsFixed; par++ )
    _propagatorsFixed &= getPar( par ).isFixed();

  //_propagatorsFixed = false; // DEBUG only.

  // Check if all amplitudes and propagators are fixed.
  _integralsFixed = _propagatorsFixed;
  for ( int par = 23; ( par < 53 ) && _integralsFixed; par++ )
    _integralsFixed &= getPar( par ).isFixed();

  // Decide which amplitudes will need to be recomputed at each minuit iteration.
  _unfixed = new bool[ nReso ];
  _unfixed[ 0 ] = getPar(  3 ).isReleased() || getPar(  4 ).isReleased();   // KStarm
  _unfixed[ 1 ] = _unfixed[ 0 ];                                            // KStarp
  _unfixed[ 2 ] = getPar(  5 ).isReleased() || getPar(  6 ).isReleased();   // Rho
  _unfixed[ 3 ] = getPar(  7 ).isReleased() || getPar(  8 ).isReleased();   // Omega
  _unfixed[ 4 ] = getPar(  9 ).isReleased() || getPar( 10 ).isReleased();   // F2
  _unfixed[ 5 ] = getPar( 11 ).isReleased() || getPar( 12 ).isReleased() || // K0Starm (LASS)
                  getPar( 13 ).isReleased() || getPar( 14 ).isReleased() ||
                  getPar( 15 ).isReleased() || getPar( 16 ).isReleased() ||
                  getPar( 17 ).isReleased() || getPar( 18 ).isReleased();
  _unfixed[ 6 ] = _unfixed[ 5 ];                                            // K0Starp
  _unfixed[ 7 ] = getPar( 19 ).isReleased() || getPar( 20 ).isReleased();   // K2Starm
  _unfixed[ 8 ] = _unfixed[ 7 ];                                            // K2Starp
  _unfixed[ 9 ] = getPar( 21 ).isReleased() || getPar( 22 ).isReleased();   // KStarm1680

  // Decide which amplitudes are the conjugate of their previous ones to avoid recomputing.
  _conjOfPrevious = new bool[ nReso ];
  _conjOfPrevious[ 0 ] = false; // KStarm
  _conjOfPrevious[ 1 ] = true;  // KStarp
  _conjOfPrevious[ 2 ] = false; // Rho
  _conjOfPrevious[ 3 ] = false; // Omega
  _conjOfPrevious[ 4 ] = false; // F2
  _conjOfPrevious[ 5 ] = false; // K0Starm
  _conjOfPrevious[ 6 ] = true;  // K0Starp
  _conjOfPrevious[ 7 ] = false; // K2Starm
  _conjOfPrevious[ 8 ] = true;  // K2Starp
  _conjOfPrevious[ 9 ] = false; // KStarm1680

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

  _reso[ 0 ].setReso( mKStarc()    , wKStarc()     );
  _reso[ 1 ].setReso( mKStarc()    , wKStarc()     );
  _reso[ 2 ].setReso( mRho()       , wRho()        );
  _reso[ 3 ].setReso( mOmega()     , wOmega()      );
  _reso[ 4 ].setReso( mF2()        , wF2()         );
  _reso[ 5 ].setReso( mK0Starc()   , wK0Starc()    );
  _reso[ 6 ].setReso( mK0Starc()   , wK0Starc()    );
  _reso[ 7 ].setReso( mK2Starc()   , wK2Starc()    );
  _reso[ 8 ].setReso( mK2Starc()   , wK2Starc()    );
  _reso[ 9 ].setReso( mKStarm1680(), wKStarm1680() );


  _reso[ 5 ].setLassa   ( aLass()  );
  _reso[ 5 ].setLassr   ( rLass()  );
  _reso[ 5 ].setLassB   ( bLass()  );
  _reso[ 5 ].setLassPhiB( pBLass() );
  _reso[ 5 ].setLassR   ( RLass()  );
  _reso[ 5 ].setLassPhiR( prLass() );

  _reso[ 6 ].setLassa   ( aLass()  );
  _reso[ 6 ].setLassr   ( rLass()  );
  _reso[ 6 ].setLassB   ( bLass()  );
  _reso[ 6 ].setLassPhiB( pBLass() );
  _reso[ 6 ].setLassR   ( RLass()  );
  _reso[ 6 ].setLassPhiR( prLass() );

  return;
}



double DalitzD0mix::evaluate() const throw( PdfException )
{
  if ( ( m2AC() < m2ACmin( m2AB() ) ) || ( m2AC() > m2ACmax( m2AB() ) ) )
    return 0.;

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

  return std::norm( amp ) / _norm;
}



std::complex< double > DalitzD0mix::evaluateReso( int index, double m2AB, double m2AC, double m2BC ) const
{
  if ( index == 0 ) return _reso[ 0 ].rbw ( m2AB, m2AC, m2BC );
  if ( index == 1 ) return _reso[ 1 ].rbw ( m2AB, m2AC, m2BC );
  if ( index == 2 ) return _reso[ 2 ].gs  ( m2AB, m2AC, m2BC );
  if ( index == 3 ) return _reso[ 3 ].rbw ( m2AB, m2AC, m2BC );
  if ( index == 4 ) return _reso[ 4 ].rbw ( m2AB, m2AC, m2BC );
  if ( index == 5 ) return _reso[ 5 ].lass( m2AB, m2AC, m2BC );
  if ( index == 6 ) return _reso[ 6 ].lass( m2AB, m2AC, m2BC );
  if ( index == 7 ) return _reso[ 7 ].rbw ( m2AB, m2AC, m2BC );
  if ( index == 8 ) return _reso[ 8 ].rbw ( m2AB, m2AC, m2BC );
  if ( index == 9 ) return _reso[ 9 ].rbw ( m2AB, m2AC, m2BC );

  return std::complex< double >();
}


std::complex< double > DalitzD0mix::coef( int index ) const
{
  if ( index == 0 ) return std::complex< double >( rKStarm()    , iKStarm()     );
  if ( index == 1 ) return std::complex< double >( rKStarp()    , iKStarp()     );
  if ( index == 2 ) return std::complex< double >( rRho()       , iRho()        );
  if ( index == 3 ) return std::complex< double >( rOmega()     , iOmega()      );
  if ( index == 4 ) return std::complex< double >( rF2()        , iF2()         );
  if ( index == 5 ) return std::complex< double >( rK0Starm()   , iK0Starm()    );
  if ( index == 6 ) return std::complex< double >( rK0Starp()   , iK0Starp()    );
  if ( index == 7 ) return std::complex< double >( rK2Starm()   , iK2Starm()    );
  if ( index == 8 ) return std::complex< double >( rK2Starp()   , iK2Starp()    );
  if ( index == 9 ) return std::complex< double >( rKStarm1680(), iKStarm1680() );

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
  std::complex< double > A   [ 10 ];
  std::complex< double > ABar[ 10 ];

  // Define the properties of the integration method.
  int    nBins = 300;
  double min   = pow( mKs() + mPi(), 2 );
  double max   = pow( mD0() - mPi(), 2 );
  double step  = ( max - min ) / double( nBins );

  // Define the variables at each bin.
  double m2AB;
  double m2AC;
  double m2BC;

  // Initialize the matrix elements that need to be recomputed.
  for ( int i = 0; i < 10; i++ )
    for ( int j = i; j < 10; j++ )
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
	    for ( int reso = 0; reso < 10; reso++ )
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
	    for ( int i = 0; i < 10; i++ )
	      for ( int j = i; j < 10; j++ )
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
  for ( int i = 0; i < 10; i++ )
    for ( int j = i; j < 10; j++ )
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

  for ( int i = 0; i < 10; i++ )
    for ( int j = 0; j < 10; j++ )
      {
	_i1   += conj( coef( i ) ) * _Id[ i ][ j ] * coef( j );
	_iChi += conj( coef( i ) ) * _Ix[ i ][ j ] * coef( j );
      }

  _iChi *= _qp / _i1;

  return;
}



double DalitzD0mix::kallen( double m1, double m2, double m3 )
{
  double result = 0.;
  result += pow( m1, 2 );
  result += pow( m2, 2 );
  result += pow( m3, 2 );
  result -= 2. * m1 * m2;
  result -= 2. * m1 * m3;
  result -= 2. * m2 * m3;

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
  /*
  static const EvtDalitzPlot plot( _mKs, _mPi, _mPi, _mD0 );

  EvtDalitzPoint point( _mKs, _mPi, _mPi, m2AB, m2BC, m2AC );

  // Defining resonances.
  EvtDalitzReso KStarmOld  ( plot, _BC, _AC, _VECTOR, mKStarc()    , wKStarc()    , _RBW );
  EvtDalitzReso KStarpOld  ( plot, _BC, _AB, _VECTOR, mKStarc()    , wKStarc()    , _RBW );
  EvtDalitzReso rho0Old    ( plot, _AC, _BC, _VECTOR, mRho()       , wRho()       , _GS  );
  EvtDalitzReso omegaOld   ( plot, _AC, _BC, _VECTOR, mOmega()     , wOmega()     , _RBW );
  EvtDalitzReso f2_1270Old ( plot, _AC, _BC, _TENSOR, mF2()        , wF2()        , _RBW );

  // LASS
  EvtDalitzReso K0Starm1430Old( plot, _AC, mK0Starc(), wK0Starc(), aLass(), rLass(), bLass(), pBLass(), RLass(), prLass() );
  EvtDalitzReso K0Starp1430Old( plot, _AB, mK0Starc(), wK0Starc(), aLass(), rLass(), bLass(), pBLass(), RLass(), prLass() );

  EvtDalitzReso K2Starm1430Old( plot, _BC, _AC, _TENSOR, mK2Starc()   , wK2Starc()   , _RBW );
  EvtDalitzReso K2Starp1430Old( plot, _BC, _AB, _TENSOR, mK2Starc()   , wK2Starc()   , _RBW );
  EvtDalitzReso KStarm1680Old ( plot, _BC, _AC, _VECTOR, mKStarm1680(), wKStarm1680(), _RBW );
  */

  // Defining K-matrix.
  /*
  EvtComplex fr12Old( 1.87981, -.628378 );
  EvtComplex fr13Old( 4.3242 , 2.75019  );
  EvtComplex fr14Old( 3.22336,  .271048 );
  EvtComplex fr15Old(  .0    ,  .0      );
  EvtDalitzReso Pole1Old  ( plot, _BC, "Pole1"  , _KM, fr12Old, fr13Old, fr14Old, fr15Old, -.0694725 );
  EvtDalitzReso Pole2Old  ( plot, _BC, "Pole2"  , _KM, fr12Old, fr13Old, fr14Old, fr15Old, -.0694725 );
  EvtDalitzReso Pole3Old  ( plot, _BC, "Pole3"  , _KM, fr12Old, fr13Old, fr14Old, fr15Old, -.0694725 );
  EvtDalitzReso Pole4Old  ( plot, _BC, "Pole4"  , _KM, fr12Old, fr13Old, fr14Old, fr15Old, -.0694725 );
  EvtDalitzReso Pole0Old  ( plot, _BC, "f11prod", _KM, fr12Old, fr13Old, fr14Old, fr15Old, -.0694725 );
  */

  /*
  Complex ampOld( 0., 0. );
  ampOld += Complex( rKStarm()    , iKStarm()      ) * KStarmOld     .evaluate( point );
  ampOld += Complex( rKStarp()    , iKStarp()      ) * KStarpOld     .evaluate( point );
  ampOld += Complex( rRho()       , iRho()         ) * rho0Old       .evaluate( point );
  ampOld += Complex( rOmega()     , iOmega()       ) * omegaOld      .evaluate( point );
  ampOld += Complex( rF2()        , iF2()          ) * f2_1270Old    .evaluate( point );
  ampOld += Complex( rK0Starm()   , iK0Starm()     ) * K0Starm1430Old.evaluate( point );
  ampOld += Complex( rK0Starp()   , iK0Starp()     ) * K0Starp1430Old.evaluate( point );
  ampOld += Complex( rK2Starm()   , iK2Starm()     ) * K2Starm1430Old.evaluate( point );
  ampOld += Complex( rK2Starp()   , iK2Starp()     ) * K2Starp1430Old.evaluate( point );
  ampOld += Complex( rKStarm1680(), iKStarm1680()  ) * KStarm1680Old .evaluate( point );
  ampOld += Complex( rPole0()     , iPole0()       ) * Pole0      .evaluate( point );
  ampOld += Complex( rPole1()     , iPole1()       ) * Pole1      .evaluate( point );
  ampOld += Complex( rPole2()     , iPole2()       ) * Pole2      .evaluate( point );
  ampOld += Complex( rPole3()     , iPole3()       ) * Pole3      .evaluate( point );
  ampOld += Complex( rPole4()     , iPole4()       ) * Pole4      .evaluate( point );
  */

  // Adding terms to the amplitude with their corresponding amplitude and phase terms.
  std::complex< double > amp( 0., 0. );
  for ( int reso = 0; reso < 10; reso++ )
    amp += coef( reso ) * evaluateReso( reso, m2AB, m2AC, m2BC );

  //amp += Complex( rPole1()      , iPole1()       ) * Pole1      .evaluate( point );
  //amp += Complex( rPole2()      , iPole2()       ) * Pole2      .evaluate( point );
  //amp += Complex( rPole3()      , iPole3()       ) * Pole3      .evaluate( point );
  //amp += Complex( rPole4()      , iPole4()       ) * Pole4      .evaluate( point );
  //amp += Complex( rPole0()      , iPole0()       ) * Pole0      .evaluate( point );

  //std::cout << "Debug diff total: " << ( ( amp - ampOld ).mod() > 1.e-12 ) << std::endl;

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
  return 0.;
}
