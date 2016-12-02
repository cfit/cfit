
#include <complex>

#include <cfit/dataset.hh>
#include <cfit/function.hh>
#include <cfit/random.hh>

#include <cfit/models/decay3bodybin.hh>


Decay3BodyBin::Decay3BodyBin( const Variable&        mSq12  ,
                              const Variable&        mSq13  ,
                              const Variable&        mSq23  ,
                              const BinnedAmplitude& amp    ,
                              const Binning&         binning,
                              const CoefExpr&        phi    ,
                              const PhaseSpace&      ps     ,
                              bool                   docache  )
  : DecayModel( mSq12, mSq13, mSq23, amp, ps ),
    _binning( binning ),
    _hasKappa( false ), _phi( phi ),
    _nDir( 0.0 ), _nXed( 0.0 ), _norm( 1.0 ), _fixedAmp( false ),
    _binIndex( 0 )
{
  push( phi );

  // Do calculations common to all values of variables
  //    (usually compute norm).
  if ( docache )
    cache();
}



Decay3BodyBin::Decay3BodyBin( const Variable&        mSq12  ,
                              const Variable&        mSq13  ,
                              const Variable&        mSq23  ,
                              const BinnedAmplitude& amp    ,
                              const Binning&         binning,
                              const CoefExpr&        phi    ,
                              const Parameter&       kappa  ,
                              const PhaseSpace&      ps     ,
                              bool                   docache  )
  : DecayModel( mSq12, mSq13, mSq23, amp, ps ),
    _binning( binning ),
    _hasKappa( true ), _kappa( kappa ), _phi( phi ),
    _nDir( 0.0 ), _nXed( 0.0 ), _norm( 1.0 ), _fixedAmp( false ),
    _binIndex( 0 )
{
  push( phi   );
  push( kappa );

  // Do calculations common to all values of variables
  //    (usually compute norm).
  if ( docache )
    cache();
}




Decay3BodyBin::Decay3BodyBin( const Variable&        mSq12  ,
                              const Variable&        mSq13  ,
                              const Variable&        mSq23  ,
                              const BinnedAmplitude& amp    ,
                              const Binning&         binning,
                              const CoefExpr&        phi    ,
                              const ParameterExpr&   kappa  ,
                              const PhaseSpace&      ps     ,
                              bool                   docache  )
  : DecayModel( mSq12, mSq13, mSq23, amp, ps ),
    _binning( binning ),
    _hasKappa( true ), _kappa( kappa ), _phi( phi ),
    _nDir( 0.0 ), _nXed( 0.0 ), _norm( 1.0 ), _fixedAmp( false ),
    _binIndex( 0 )
{
  push( phi   );
  push( kappa );

  // Do calculations common to all values of variables
  //    (usually compute norm).
  if ( docache )
    cache();
}





Decay3BodyBin* Decay3BodyBin::copy() const
{
  return new Decay3BodyBin( *this );
}



void Decay3BodyBin::setParExpr()
{
  _phi.setPars( _parMap );

  // Propagate the kappa parameter value if necessary.
  if ( _hasKappa )
    _kappa.setPars( _parMap );
}



void Decay3BodyBin::cacheNormComponents()
{
  // If the amplitude is fixed and the components have already
  //    been computed, just return without recomputing anything.
  if ( _fixedAmp )
    return;

  // Define the properties of the integration method.
  const unsigned& nBins = 100;
  const double min   = _ps.mSq12min();
  const double max   = _ps.mSq12max();
  const double step  = ( max - min ) / double( nBins );

  // Determine whether the bins in the 100x100 bins should be cached.
  //    Cache them if it has not yet been cached.
  bool needToCache = _binCache.empty();
  if ( needToCache )
    _binCache.resize( std::pow( nBins, 2 ) );

  double mSq12;
  double mSq13;

  // Initialize the integrated efficiencies.
  std::vector< double > effs;
  effs.resize( _amp.nBins() );
  for ( unsigned bin = 0; bin < _amp.nBins(); ++bin )
    effs[ bin ] = 0.0;

  unsigned bin;
  // Integrate only over the region mSq13 < mSq12 (positive bins).
  for ( unsigned binX = 0; binX < nBins; ++binX )
    for ( unsigned binY = 0; binY < binX; ++binY )
    {
      mSq12 = min + step * ( binX + 0.5 );
      mSq13 = min + step * ( binY + 0.5 );

      if ( _ps.contains( mSq12, mSq13 ) )
      {
        if ( needToCache )
        {
          bin = std::abs( _binning.bin( mSq12, mSq13 ) );

          if ( needToCache )
            _binCache[ nBins * binX + binY ] = bin;
        }
        else
        {
          bin = _binCache[ nBins * binX + binY ];
        }

        // std::cout << "DEBUG Decay3BodyBin: " << mSq12 << " " << mSq13 << " " << _binning.bin( mSq12, mSq13 ) << std::endl;
        // bin = std::abs( _binning.bin( mSq12, mSq13 ) );
        effs[ bin - 1 ] += evaluateFuncs( mSq12, mSq13 );
      }
    }

  for ( unsigned bin = 0; bin < _amp.nBins(); ++bin )
    effs[ bin ] *= std::pow( step, 2 );

  // Initialize the value of the norm components and the norm.
  _nDir = 0.0;
  _nXed = 0.0;
  _norm = 0.0;

  std::vector< Parameter >::const_iterator tpb = _amp.tpb().begin();
  std::vector< Parameter >::const_iterator tmb = _amp.tmb().begin();
  std::vector< CoefExpr  >::const_iterator xb  = _amp.xb() .begin();
  std::vector< double    >::const_iterator eff = effs      .begin();

  double                 ef;
  double                 tp;
  double                 tm;
  std::complex< double > tx;

  // Calculate the norm components.
  while ( tpb != _amp.tpb().end() )
  {
    ef = *eff++;

    tp = (tpb++)->value();
    tm = (tmb++)->value();
    tx = std::sqrt( tp * tm ) * (xb++)->evaluate();

    _nDir += ef * ( tp + tm );
    _nXed += ef * 2.0 * std::real( tx );
  }

  return;
}


void Decay3BodyBin::cache()
{
  // Compute the norm components, only if the amplitude
  //    is not fixed or they have not yet been computed.
  cacheNormComponents();

  // _fixedAmp = true; // _amp.isFixed();

  // Calculate the norm.
  const std::complex< double >& vz     = z();
  const double&                 vKappa = kappa();
  _norm = _nDir * ( 1.0 + std::norm( vz ) ) + 2.0 * vKappa * std::real( vz * _nXed );

  return;
}





const std::map< unsigned, std::vector< double > > Decay3BodyBin::cacheReal( const Dataset& data )
{
  std::map< unsigned, std::vector< double > > cached;

  // Get an index for the cached bin.
  _binIndex = _cacheIdxReal++;

  const std::string& mSq12name = getVar( 0 ).name();
  const std::string& mSq13name = getVar( 1 ).name();

  double mSq12;
  double mSq13;

  const std::size_t& size = data.size();
  for ( std::size_t entry = 0; entry < size; ++entry )
  {
    mSq12 = data.value( mSq12name, entry );
    mSq13 = data.value( mSq13name, entry );

    cached[ _binIndex ].push_back( _binning.bin( mSq12, mSq13 ) );
  }

  return cached;
}






// Unnormalized evaluation.
const double Decay3BodyBin::evaluateUnnorm( const double& mSq12, const double& mSq13 ) const throw( PdfException )
{
  // Calculate the bin number from the binning.
  int bin = _binning.bin( mSq12, mSq13 );

  const std::tuple< double, double, std::complex< double > >&& tx = _amp.evaluate( bin );

  const std::complex< double >&& vz     = z();
  const double&&                 vKappa = kappa();

  const double&&                 funcs  = evaluateFuncs( mSq12, mSq13 );

  return ( std::get< 0 >( tx )                   +
           std::get< 1 >( tx ) * std::norm( vz ) +
           2.0 * vKappa * std::real( vz * std::get< 2 >( tx ) ) ) * funcs;
}


const double Decay3BodyBin::evaluate( const double& mSq12, const double& mSq13, const double& ) const throw( PdfException )
{
  // Ignore any eventual 3rd argument.
  return evaluate( mSq12, mSq13 );
}



const double Decay3BodyBin::evaluate( const double& mSq12, const double& mSq13 ) const throw( PdfException )
{
  return evaluateUnnorm( mSq12, mSq13 ) / _norm;
}



const double Decay3BodyBin::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  if ( vars.size() != 2 && vars.size() != 3 )
    throw PdfException( "Decay3BodyBin can only take 2 or 3 arguments." );

  // In any case, ignore the 3rd argument, mSq23.
  return evaluate( vars[ 0 ], vars[ 1 ] );
}



const double Decay3BodyBin::evaluate( const std::vector< double >&                 vars  ,
                                      const std::vector< double >&                 cacheR,
                                      const std::vector< std::complex< double > >& cacheC ) const throw( PdfException )
{
  if ( cacheR.empty() )
    return evaluate( vars );

  const std::size_t& size = vars.size();

  if ( ( size != 2 ) && ( size != 3 ) )
    throw PdfException( "Decay3BodyBin can only take either 2 or 3 arguments." );

  int bin = cacheR[ _binIndex ];

  const std::complex< double >&& vz     = z();
  const double&&                 vKappa = kappa();

  // Evaluate the functions that describe the efficiency.
  double funcs = 0.0;
  if ( size == 2 )
    funcs = evaluateFuncs( vars[ 0 ], vars[ 1 ] );
  else if ( size == 3 )
    funcs = evaluateFuncs( vars[ 0 ], vars[ 1 ], vars[ 2 ] );

  const std::tuple< double, double, std::complex< double > >&& nx = _amp.evaluate( bin );

  return ( std::get< 0 >( nx )                   +
           std::get< 1 >( nx ) * std::norm( vz ) +
           2.0 * vKappa * std::real( vz * std::get< 2 >( nx ) ) ) * funcs / _norm;
}




// No need to append an operator, since it can only be multiplication.
const Decay3BodyBin& Decay3BodyBin::operator*=( const Function& right ) throw( PdfException )
{
  // For this case, do not check that the function does not depend on any variables that the model does not.
  // The norm has to be calculated correctly, though.

  // Consider the function parameters as own ones.
  const std::map< std::string, Parameter >& parMap = right.getParMap();
  _parMap.insert( parMap.begin(), parMap.end() );

  // Append the function to the functions vector.
  _funcs.push_back( right );

  // Recompute the norm, since the pdf shape has changed under this operation.
  cache();

  return *this;
}



const Decay3BodyBin operator*( Decay3BodyBin left, const Function& right )
{
  // For this case, do not check that the function does not depend on any variables that the model does not.
  // The norm has to be calculated correctly, though.

  // Consider the function parameters as own ones.
  const std::map< std::string, Parameter >& parMap = right.getParMap();
  left._parMap.insert( parMap.begin(), parMap.end() );

  // Append the function to the functions vector.
  left._funcs.push_back( right );

  // Recompute the norm, since the pdf shape has changed under this operation.
  left.cache();

  return left;
}



const Decay3BodyBin operator*( const Function& left, Decay3BodyBin right )
{
  // For this case, do not check that the function does not depend on any variables that the model does not.
  // The norm has to be calculated correctly, though.

  // Consider the function parameters as own ones.
  const std::map< std::string, Parameter >& parMap = left.getParMap();
  right._parMap.insert( parMap.begin(), parMap.end() );

  // Append the function to the functions vector.
  right._funcs.push_back( left );

  // Recompute the norm, since the pdf shape has changed under this operation.
  right.cache();

  return right;
}

