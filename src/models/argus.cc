
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
  const double& cSq   = std::pow( c()  , 2 );
  const double& chiSq = std::pow( chi(), 2 );

  // For the specific case when chi = 0, the norm is c^2/3.
  if ( chiSq == 0.0 )
  {
    _norm = cSq / 3.0;
    return;
  }

  double argmax = chiSq;
  double argmin = 0.0;

  if ( _hasLower )
  {
    const double& lowerlimit = std::max( _lower, 0.0 );
    argmax = chiSq * ( 1.0 - std::pow( lowerlimit, 2 ) / cSq );
  }

  if ( _hasUpper )
  {
    const double& upperlimit = std::min( _upper, c() );
    argmin = chiSq * ( 1.0 - std::pow( upperlimit, 2 ) / cSq );
  }

  const double& chi3 = std::pow( chi(), 3 );

  // Since gamma_p( a, x ) is normalized to Gamma( a ), multiply by sqrt( pi / 2 ),
  //    which is Gamma( 3/2 ).
  _norm  = cSq / ( 2.0 * chi3 ) * std::sqrt( M_PI / 2.0 );
  _norm *= ( Math::gamma_p( 1.5, argmax ) - Math::gamma_p( 1.5, argmin ) );
}


double Argus::evaluate( double x ) const throw( PdfException )
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


double Argus::evaluate() const throw( PdfException )
{
  return evaluate( getVar( 0 ).value() );
}


double Argus::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  return evaluate( vars[ 0 ] );
}
