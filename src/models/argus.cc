
#include <cfit/models/argus.hh>
#include <cfit/math.hh>

Argus::Argus( const Variable& x, const Parameter& c, const Parameter& chi )
  : _hasLower( false ), _hasUpper( false ), _lower( 0.0 ), _upper( 0.0 )
{
  push( x );

  push( c   );
  push( chi );

  cache();
}


Argus* Argus::copy() const
{
  return new Argus( *this );
}


double Argus::c()  const
{
  return getPar( 0 ).value();
}


double Argus::chi() const
{
  return getPar( 1 ).value();
}


void Argus::setLowerLimit( const double& lower )
{
  if ( lower < 0 )
    throw PdfException( "Cannot set the lower limit of the Argus distribution to anything smaller than 0." );

  _hasLower = true;
  _lower    = lower;
}


void Argus::setUpperLimit( const double& upper )
{
  if ( upper < 0.0 )
    throw PdfException( "Cannot set the upper limit of the generalized Argus distribution to anything smaller than 0." );

  _hasUpper = true;
  _upper    = upper;
}


void Argus::setLimits( const double& lower, const double& upper )
{
  if ( lower < 0 )
    throw PdfException( "Cannot set the lower limit of the Argus distribution to anything smaller than 0." );

  if ( upper < 0.0 )
    throw PdfException( "Cannot set the upper limit of the generalized Argus distribution to anything smaller than 0." );

  _hasLower = true;
  _hasUpper = true;
  _lower    = lower;
  _upper    = upper;
}


void Argus::unsetLowerLimit()
{
  _hasLower = false;
}


void Argus::unsetUpperLimit()
{
  _hasUpper = false;
}


void Argus::unsetLimits()
{
  _hasLower = false;
  _hasUpper = false;
}


void Argus::cache()
{
  const double& vc    = c();
  const double& cSq   = std::pow( c()  , 2 );
  const double& chiSq = std::pow( chi(), 2 );

  const double& lower = _hasLower ? std::max( _lower, 0.0 ) : 0.0;
  const double& upper = _hasUpper ? std::min( _upper, vc  ) : vc;

  const double& argmax = _hasLower ? 1.0 - std::pow( lower / vc, 2 ) : 1.0;
  const double& argmin = _hasUpper ? 1.0 - std::pow( upper / vc, 2 ) : 0.0;

  // For the specific case when chi = 0, the norm is
  //    c^2/3 ( 1 - x^2 / c^2 )^(3/2) between upper and lower.
  if ( chiSq == 0.0 )
  {
    _norm = cSq / 3.0 * ( std::pow( argmax, 1.5 ) - std::pow( argmin, 1.5 ) );
    return;
  }

  const double& chi3 = std::pow( chi(), 3 );

  // Since gamma_p( a, x ) is normalized to Gamma( a ), multiply by sqrt( pi / 2 ),
  //    which is Gamma( 3/2 ).
  _norm  = cSq / ( 2.0 * chi3 ) * std::sqrt( M_PI / 2.0 );
  _norm *= ( Math::gamma_p( 1.5, chiSq * argmax ) - Math::gamma_p( 1.5, chiSq * argmin ) );
}


const double Argus::evaluate( const double& x ) const throw( PdfException )
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

  return x * std::sqrt( diff ) * std::exp( - chiSq * diff ) / _norm;
}


const double Argus::evaluate() const throw( PdfException )
{
  return evaluate( getVar( 0 ).value() );
}


const double Argus::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  return evaluate( vars[ 0 ] );
}


const double Argus::area( const double& min, const double& max ) const throw( PdfException )
{
  const double& vc    = c();
  const double& cSq   = std::pow( c()  , 2 );
  const double& chiSq = std::pow( chi(), 2 );

  const double& lower = _hasLower ? std::max( _lower, 0.0 ) : 0.0;
  const double& upper = _hasUpper ? std::min( _upper, vc  ) : vc;

  const double& xmin = std::max( min, lower );
  const double& xmax = std::min( max, upper );

  const double& argmax = _hasLower ? 1.0 - std::pow( xmin / vc, 2 ) : 1.0;
  const double& argmin = _hasUpper ? 1.0 - std::pow( xmax / vc, 2 ) : 0.0;

  // For the specific case when chi = 0, the norm is
  //    c^2/3 ( 1 - x^2 / c^2 )^(3/2) between upper and lower.
  if ( chiSq == 0.0 )
    return cSq / 3.0 * ( std::pow( argmax, 1.5 ) - std::pow( argmin, 1.5 ) ) / _norm;

  const double& chi3 = std::pow( chi(), 3 );

  // Since gamma_p( a, x ) is normalized to Gamma( a ), multiply by sqrt( pi / 2 ),
  //    which is Gamma( 3/2 ).
  double retval = cSq / ( 2.0 * chi3 ) * std::sqrt( M_PI / 2.0 );
  retval *= ( Math::gamma_p( 1.5, chiSq * argmax ) - Math::gamma_p( 1.5, chiSq * argmin ) ) / _norm;

  return retval;
}

