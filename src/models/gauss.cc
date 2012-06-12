
#include <cfit/models/gauss.hh>
#include <cfit/math.hh>


Gauss::Gauss( const Variable& x, const Parameter& mu, const Parameter& sigma )
  : _hasLower( false ), _hasUpper( false ), _lower( 0.0 ), _upper( 0.0 )
{
  push( x );

  push( mu    );
  push( sigma );

  cache();
}


double Gauss::mu()  const
{
  return getPar( 0 ).value();
}


double Gauss::sigma() const
{
  return getPar( 1 ).value();
}


void Gauss::setLowerLimit( const double& lower )
{
  _hasLower = true;
  _lower    = lower;

  cache();
}


void Gauss::setUpperLimit( const double& upper )
{
  _hasUpper = true;
  _upper    = upper;

  cache();
}


void Gauss::setLimits( const double& lower, const double& upper )
{
  _hasLower = true;
  _hasUpper = true;
  _lower    = lower;
  _upper    = upper;

  cache();
}


void Gauss::unsetLowerLimit()
{
  _hasLower = false;

  cache();
}


void Gauss::unsetUpperLimit()
{
  _hasUpper = false;

  cache();
}


void Gauss::unsetLimits()
{
  _hasLower = false;
  _hasUpper = false;

  cache();
}


void Gauss::cache()
{
  const double& vmu    = mu();
  const double& vsigma = sigma();
  const double& sqrt2  = std::sqrt( 2.0 );

  double argmin = 0.0;
  if ( _hasLower )
    argmin = 1.0 + Math::erf( ( _lower - vmu ) / ( vsigma * sqrt2 ) );

  double argmax = 2.0;
  if ( _hasUpper )
    argmax = 1.0 + Math::erf( ( _upper - vmu ) / ( vsigma * sqrt2 ) );

  const double& factor = vsigma * std::sqrt( M_PI / 2.0 );
  _norm = factor * ( argmax - argmin );
}


double Gauss::evaluate( double x ) const throw( PdfException )
{
  return std::exp( - 0.5 * pow( x - mu(), 2 ) / pow( sigma(), 2 ) ) / _norm;
}


double Gauss::evaluate() const throw( PdfException )
{
  return evaluate( getVar( 0 ).value() );
}


double Gauss::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  return evaluate( vars[ 0 ] );
}
