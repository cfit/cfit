
#include <cfit/models/gauss.hh>

Gauss::Gauss( const Variable& x, const Parameter& mu, const Parameter& sigma )
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
