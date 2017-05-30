
#include <cfit/models/genargusgauss.hh>
#include <cfit/math.hh>

#include <cfit/dataset.hh>

#include <cfit/random.hh>

GenArgusGauss::GenArgusGauss( const Variable&  x ,
                              const Parameter& c , const Parameter& chi, const Parameter& p,
                              const Parameter& mu, const Parameter& sigma )
  : _c( c ), _chi( chi ), _p( p ), _mu( mu ), _sigma( sigma ),
    _hasLower( false ), _hasUpper( false ), _lower( 0.0 ), _upper( 0.0 ),
    _normGenArgus( 0.0 ), _normGauss( 0.0 ), _norm( 0.0 ), _doCache( false ), _cacheIdx( 0 )
{
  push( x );

  push( c     );
  push( chi   );
  push( p     );
  push( mu    );
  push( sigma );

  cache();
}


GenArgusGauss::GenArgusGauss( const Variable&  x ,
                              const ParameterExpr& c , const ParameterExpr& chi, const ParameterExpr& p,
                              const ParameterExpr& mu, const ParameterExpr& sigma )
  : _c( c ), _chi( chi ), _p( p ), _mu( mu ), _sigma( sigma ),
    _hasLower( false ), _hasUpper( false ), _lower( 0.0 ), _upper( 0.0 ),
    _normGenArgus( 0.0 ), _normGauss( 0.0 ), _norm( 0.0 ), _doCache( false ), _cacheIdx( 0 )
{
  push( x );

  push( c     );
  push( chi   );
  push( p     );
  push( mu    );
  push( sigma );

  cache();
}


GenArgusGauss* GenArgusGauss::copy() const
{
  return new GenArgusGauss( *this );
}


double GenArgusGauss::c()     const { return _c    .evaluate(); }
double GenArgusGauss::chi()   const { return _chi  .evaluate(); }
double GenArgusGauss::p()     const { return _p    .evaluate(); }
double GenArgusGauss::mu()    const { return _mu   .evaluate(); }
double GenArgusGauss::sigma() const { return _sigma.evaluate(); }



void GenArgusGauss::setLowerLimit( const double& lower )
{
  if ( lower < 0 )
    throw PdfException( "Cannot set the lower limit of the generalized Argus distribution to anything smaller than 0." );

  _hasLower = true;
  _lower    = lower;

  cache();
}


void GenArgusGauss::setUpperLimit( const double& upper )
{
  if ( upper < 0.0 )
    throw PdfException( "Cannot set the upper limit of the generalized Argus distribution to anything smaller than 0." );

  _hasUpper = true;
  _upper    = upper;

  cache();
}


void GenArgusGauss::setLimits( const double& lower, const double& upper )
{
  if ( lower < 0.0 )
    throw PdfException( "Cannot set the lower limit of the generalized Argus distribution to anything smaller than 0." );

  if ( upper < 0.0 )
    throw PdfException( "Cannot set the upper limit of the generalized Argus distribution to anything smaller than 0." );

  _hasLower = true;
  _hasUpper = true;
  _lower    = lower;
  _upper    = upper;

  cache();
}


void GenArgusGauss::unsetLowerLimit()
{
  _hasLower = false;

  cache();
}


void GenArgusGauss::unsetUpperLimit()
{
  _hasUpper = false;

  cache();
}


void GenArgusGauss::unsetLimits()
{
  _hasLower = false;
  _hasUpper = false;

  cache();
}


void GenArgusGauss::cache()
{
  const double& vc     = c();
  const double& cSq    = std::pow( vc   , 2 );
  const double& chiSq  = std::pow( chi(), 2 );
  const double& pPlus1 = p() + 1.0;


  // For the specific case when chi = 0, the norm is
  //    c^2 / ( 2 ( p + 1 ) ) ( 1 - x^2 / c^2 )^( p + 1 )
  //    between upper and lower.
  if ( chiSq == 0.0 )
  {
    _normGenArgus = cSq / ( 2.0 * pPlus1 );
    return;
  }

  const double& chiPow = std::pow( chiSq, pPlus1 );

  // Since gamma_p( a, x ) is normalized to Gamma( a ), multiply by Gamma( p + 1 ).
  _normGenArgus  = cSq / ( 2.0 * chiPow ) * Math::gamma( pPlus1 );
  _normGenArgus *= ( Math::gamma_p( pPlus1, chiSq ) - Math::gamma_p( pPlus1, 0.0 ) );

  // Evaluate the norm of the Gaussian.
  _normGauss = sigma() * std::sqrt( 2.0 * M_PI );

  // If no lower or upper limit has been set, the pdf is already normalized to 1.
  if ( ( ! _hasLower ) && ( ! _hasUpper ) )
  {
    _norm = 1.0;
    return;
  }

  // Normalize the pdf if it has limits.
  _norm = 0.0;
  const unsigned& nbin = 300;
  const double& lower  = _hasLower ? std::max( _lower, 0.0 ) : 0.0;
  const double& upper  = _hasUpper ? std::min( _upper, vc  ) : vc;

  const double& dx = ( upper - lower ) / double( nbin );
  for ( double x = lower; x < upper; x += dx )
    _norm += genargusgauss( x );

  _norm *= dx;
}


// Cache the values of the pdf at every point in the dataset, if the parameters are fixed.
const std::map< unsigned, std::vector< double > > GenArgusGauss::cacheReal( const Dataset& data )
{
  // Determine whether the amplitudes should be cached, i.e. only if all their parameters are fixed.
  _doCache = true;
  for ( unsigned par = 0; par < 5; ++par )
    _doCache &= getPar( par ).isFixed();

  std::map< unsigned, std::vector< double > > cached;

  if ( ! _doCache )
    return cached;

  // Get an index for the cached complex amplitudes.
  _cacheIdx = _cacheIdxReal++;

  const std::string& varname = getVar( 0 ).name();
  const std::size_t& size = data.size();
  for ( std::size_t entry = 0; entry < size; ++entry )
    cached[ _cacheIdx ].push_back( evaluate( data.value( varname, entry ) ) );

  return cached;
}



const double GenArgusGauss::genargus( const double& x ) const
{
  const double& vc   = c();
  const double& vchi = chi();

  if ( _hasLower && ( x < _lower ) )
    return 0.0;

  if ( _hasUpper && ( x > _upper ) )
    return 0.0;

  if ( ( x < 0.0 ) || ( x > vc ) )
    return 0.0;

  const double& cSq   = std::pow( vc  , 2 );
  const double& chiSq = std::pow( vchi, 2 );

  const double& xSq = std::pow( x, 2 );

  const double& diff = 1.0 - xSq / cSq;

  return x * std::pow( diff, p() ) * std::exp( - chiSq * diff ) / _normGenArgus;
}


// Argus core with c = 1. I.e., Argus / sqrt( 1 - x^2 )
const double GenArgusGauss::genarguscore( const double& x ) const
{
  if ( ( x < 0.0 ) || ( x > 1.0 ) )
    return 0.0;

  const double& vchi  = chi();
  const double& chiSq = std::pow( vchi, 2 );
  const double& xSq   = std::pow( x, 2 );
  const double& diff  = 1.0 - xSq;

  return x * std::pow( diff, p() - 0.5 ) * std::exp( - chiSq * diff ) / _normGenArgus;
}


const double GenArgusGauss::gauss( const double& x ) const
{
  return std::exp( - 0.5 * pow( x - mu(), 2 ) / pow( sigma(), 2 ) ) / _normGauss;
}


const double GenArgusGauss::genargusgauss( const double& x ) const
{
  // Evaluate the convolution of the generalized Argus pdf with a Gaussian
  //    using Gauss-Chebyshev quadrature.
  const double& vc  = c();
  const double& cSq = std::pow( vc, 2.0 );

  // Determine the required polynomial degree as the number of Gaussian sigmas
  //    in the Argus range of definition (0,c).
  const unsigned& degree = 3.0 * vc / sigma();

  // Determine the maximum relevant root to be used. For index k > degree/2,
  //    the Chebyshev polynomial root is negative or zero, and the Argus pdf
  //    evaluates to zero at these points. There's no need to evaluate it there.
  const unsigned& maxRoot = degree / 2;

  double root   = 0.0;
  double weight = 0.0;
  double integ  = 0.0;
  for ( unsigned k = 1; k <= maxRoot; ++k )
  {
    // Evaluate the polynomial roots and weights.
    root   =           std::cos( ( M_PI * k ) / ( degree + 1.0 ) );
    weight = std::pow( std::sin( ( M_PI * k ) / ( degree + 1.0 ) ), 2 ) * M_PI / ( degree + 1.0 );

    // Accumulate the value of the integral.
    integ += weight * genarguscore( root ) * gauss( x - root * vc );
  }
  integ *= cSq;

  return integ;
}





const double GenArgusGauss::evaluate( const double& x ) const throw( PdfException )
{
  return genargusgauss( x ) / _norm;
}


const double GenArgusGauss::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  return evaluate( vars[ 0 ] );
}


const double GenArgusGauss::evaluate( const std::vector< double >&                 vars  ,
                                      const std::vector< double >&                 cacheR,
                                      const std::vector< std::complex< double > >& cacheC ) const throw( PdfException )
{
  if ( ! _doCache )
    return evaluate( vars );

  return cacheR[ _cacheIdx ];
}

void GenArgusGauss::setParExpr()
{
  _c    .setPars( _parMap );
  _chi  .setPars( _parMap );
  _p    .setPars( _parMap );
  _mu   .setPars( _parMap );
  _sigma.setPars( _parMap );

  // When parameters are set, no cached values are valid anymore.
  _areas.clear();
}


const double GenArgusGauss::area( const double& min, const double& max ) const throw( PdfException )
{
  std::pair< double, double > range = std::make_pair( min, max );

  // If any cached value can be used, use it.
  if ( _areas.count( range ) )
    return _areas.at( range );

  // Set the limits of integration.
  const double& xmin = _hasLower ? std::max( min, _lower ) : min;
  const double& xmax = _hasUpper ? std::min( max, _upper ) : max;

  // Define the precision.
  const unsigned& nbin = 300;
  const double& dx = ( xmax - xmin ) / double( nbin );

  // Do the integral.
  double retval = 0.0;
  for ( double x = xmin; x < xmax; x += dx )
    retval += evaluate( x );

  // Cache the calculation of the area for this range.
  _areas[ range ] = retval * dx;

  return retval * dx;
}


const double GenArgusGauss::generateArgus() const
{
  const double& vc     = c();
  const double& chiSq  = std::pow( chi(), 2 );
  const double& pPlus1 = p() + 1.0;

  // Generate a flat random number.
  const double& unif = Random::flat();

  std::map< std::string, double > gen;

  // Deal with the special case where chi = 0.
  if ( chiSq == 0.0 )
    return vc * std::sqrt( 1.0 - std::pow( unif, 1.0 / pPlus1 ) );

  // Define a useful term for the random number generation.
  const double& rndTerm = 1.0 - unif * Math::gamma_p( pPlus1, chiSq );

  // Generate a random number.
  return vc * std::sqrt( 1.0 - Math::invgamma_q( pPlus1, rndTerm ) / chiSq );
}


const std::map< std::string, double > GenArgusGauss::generate() const throw( PdfException )
{
  double x = 0.0;

  bool withinLimits = false;
  while ( ! withinLimits )
  {
    // Add an Argus distributed variable and a Gaussian variable.
    x = generateArgus() + Random::normal( mu(), sigma() );

    // Check if the result is within limits, if any.
    withinLimits = true;
    withinLimits &= ( ! _hasLower ) || ( x > _lower );
    withinLimits &= ( ! _hasUpper ) || ( x < _upper );
  }

  // Fill a map with the generated value.
  std::map< std::string, double > gen;
  gen[ getVar( 0 ).name() ] = x;

  return gen;
}
