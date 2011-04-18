
#include <cfit/models/gauss.hh>

Gauss::Gauss( const Variable& x, const Parameter& mean, const Parameter& sigma )
{
  push( x );

  push( mean  );
  push( sigma );
}


double Gauss::mean()  const
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
  return 1. / ( sigma() * sqrt( 2. * M_PI ) ) * exp( - .5 * pow( x - mean(), 2 ) / pow( sigma(), 2 ) );
}


double Gauss::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  return evaluate( vars[ 0 ] );
}
