
#include <cfit/models/gauss.hh>

Gauss::Gauss( const Variable& x, const Parameter& mu, const Parameter& sigma )
  : _hasLower( false ), _hasUpper( false ), _lower( 0.0 ), _upper( 0.0 )
{
  push( x );

  push( mu    );
  push( sigma );
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
}


void Gauss::setUpperLimit( const double& upper )
{
  _hasUpper = true;
  _upper    = upper;
}


void Gauss::setLimits( const double& lower, const double& upper )
{
  _hasLower = true;
  _hasUpper = true;
  _lower    = lower;
  _upper    = upper;
}


void Gauss::unsetLowerLimit()
{
  _hasLower = false;
}


void Gauss::unsetUpperLimit()
{
  _hasUpper = false;
}


void Gauss::unsetLimits()
{
  _hasLower = false;
  _hasUpper = false;
}



double Gauss::evaluate() const throw( PdfException )
{
  return evaluate( getVar( 0 ).value() );
}


double Gauss::evaluate( double x ) const throw( PdfException )
{
  return 1. / ( sigma() * sqrt( 2. * M_PI ) ) * exp( - .5 * pow( x - mu(), 2 ) / pow( sigma(), 2 ) );
}


double Gauss::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  return evaluate( vars[ 0 ] );
}
