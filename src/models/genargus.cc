
#include <cfit/models/genargus.hh>
#include <cfit/math.hh>

#include <cfit/random.hh>

GenArgus::GenArgus( const Variable& x, const Parameter& c, const Parameter& chi, const Parameter& p )
  : _hasLower( false ), _hasUpper( false ), _lower( 0.0 ), _upper( 0.0 )
{
  push( x );

  push( c   );
  push( chi );
  push( p   );

  cache();
}


GenArgus* GenArgus::copy() const
{
  return new GenArgus( *this );
}


double GenArgus::c()  const
{
  return getPar( 0 ).value();
}


double GenArgus::chi() const
{
  return getPar( 1 ).value();
}


double GenArgus::p() const
{
  return getPar( 2 ).value();
}


void GenArgus::setLowerLimit( const double& lower )
{
  if ( lower < 0 )
    throw PdfException( "Cannot set the lower limit of the generalized Argus distribution to anything smaller than 0." );

  _hasLower = true;
  _lower    = lower;

  cache();
}


void GenArgus::setUpperLimit( const double& upper )
{
  if ( upper < 0.0 )
    throw PdfException( "Cannot set the upper limit of the generalized Argus distribution to anything smaller than 0." );

  _hasUpper = true;
  _upper    = upper;

  cache();
}


void GenArgus::setLimits( const double& lower, const double& upper )
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


void GenArgus::unsetLowerLimit()
{
  _hasLower = false;

  cache();
}


void GenArgus::unsetUpperLimit()
{
  _hasUpper = false;

  cache();
}


void GenArgus::unsetLimits()
{
  _hasLower = false;
  _hasUpper = false;

  cache();
}


void GenArgus::cache()
{
  const double& vc    = c();
  const double& cSq    = std::pow( c()  , 2 );
  const double& chiSq  = std::pow( chi(), 2 );
  const double& pPlus1 = p() + 1.0;

  const double& lower = _hasLower ? std::max( _lower, 0.0 ) : 0.0;
  const double& upper = _hasUpper ? std::min( _upper, vc  ) : vc;

  const double& argmax = _hasLower ? 1.0 - std::pow( lower / vc, 2 ) : 1.0;
  const double& argmin = _hasUpper ? 1.0 - std::pow( upper / vc, 2 ) : 0.0;

  // For the specific case when chi = 0, the norm is
  //    c^2 / ( 2 ( p + 1 ) ) ( 1 - x^2 / c^2 )^( p + 1 )
  //    between upper and lower.
  if ( chiSq == 0.0 )
  {
    _norm = cSq / ( 2.0 * pPlus1 ) * ( std::pow( argmax, pPlus1 ) - std::pow( argmin, pPlus1 ) );
    return;
  }

  const double& chiPow = std::pow( chiSq, pPlus1 );

  // Since gamma_p( a, x ) is normalized to Gamma( a ), multiply by Gamma( p + 1 ).
  _norm  = cSq / ( 2.0 * chiPow ) * Math::gamma( pPlus1 );
  _norm *= ( Math::gamma_p( pPlus1, chiSq * argmax ) - Math::gamma_p( pPlus1, chiSq * argmin ) );
}


const double GenArgus::evaluate( const double& x ) const throw( PdfException )
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

  return x * std::pow( diff, p() ) * std::exp( - chiSq * diff ) / _norm;
}


const double GenArgus::evaluate() const throw( PdfException )
{
  return evaluate( getVar( 0 ).value() );
}


const double GenArgus::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  return evaluate( vars[ 0 ] );
}


const double GenArgus::area( const double& min, const double& max ) const throw( PdfException )
{
  const double& vc     = c();
  const double& cSq    = std::pow( c()  , 2 );
  const double& chiSq  = std::pow( chi(), 2 );
  const double& pPlus1 = p() + 1.0;

  const double& lower = _hasLower ? std::max( _lower, 0.0 ) : 0.0;
  const double& upper = _hasUpper ? std::min( _upper, vc  ) : vc;

  const double& xmin = std::max( min, lower );
  const double& xmax = std::min( max, upper );

  const double& argmax = _hasLower ? 1.0 - std::pow( xmin / vc, 2 ) : 1.0;
  const double& argmin = _hasUpper ? 1.0 - std::pow( xmax / vc, 2 ) : 0.0;

  // For the specific case when chi = 0, the norm is
  //    c^2 / ( 2 ( p + 1 ) ) ( 1 - x^2 / c^2 )^( p + 1 )
  //    between upper and lower.
  if ( chiSq == 0.0 )
    return cSq / ( 2.0 * pPlus1 ) * ( std::pow( argmax, pPlus1 ) - std::pow( argmin, pPlus1 ) );

  const double& chiPow = std::pow( chiSq, pPlus1 );

  // Since gamma_p( a, x ) is normalized to Gamma( a ), multiply by Gamma( p + 1 ).
  double retval = cSq / ( 2.0 * chiPow ) * Math::gamma( pPlus1 );
  retval *= ( Math::gamma_p( pPlus1, chiSq * argmax ) - Math::gamma_p( pPlus1, chiSq * argmin ) ) / _norm;

  return retval;
}


const std::map< std::string, double > GenArgus::generate() const throw( PdfException )
{
  const double& vc     = c();
  const double& cSq    = std::pow( c()  , 2 );
  const double& chiSq  = std::pow( chi(), 2 );
  const double& pPlus1 = p() + 1.0;

  const double& lower  = _hasLower ? std::max( _lower, 0.0 ) : 0.0;
  const double& argmax = _hasLower ? 1.0 - std::pow( lower / vc, 2 ) : 1.0;

  // Generate a flat random number.
  const double& unif = Random::flat();

  double genVal = 0.0;
  std::map< std::string, double > gen;

  // Deal with the special case where chi = 0.
  if ( chiSq == 0.0 )
  {
    genVal = vc * std::sqrt( 1.0 - std::pow( argmax - _norm * unif * 2.0 * pPlus1 / cSq, 1.0 / pPlus1 ) );
    gen[ getVar( 0 ).name() ] = genVal;

    return gen;
  }

  // Define useful terms for the random number generation.
  const double& minTerm = Math::gamma_q( pPlus1, chiSq * argmax );
  const double& rndTerm = _norm * unif * 2.0 * std::pow( chiSq, pPlus1 ) / cSq / std::tgamma( pPlus1 );

  // Generate a random number.
  genVal = vc * std::sqrt( 1.0 - Math::invgamma_q( pPlus1, minTerm + rndTerm ) / chiSq );
  gen[ getVar( 0 ).name() ] = genVal;

  return gen;
}

