
#include <cfit/math.hh>
#include <cfit/models/crystalball.hh>

CrystalBall::CrystalBall( const Variable& x,
			  const Parameter& mu, const Parameter& sigma, const Parameter& alpha, const Parameter& n )
{
  push( x );

  push( mu    );
  push( sigma );
  push( alpha );
  push( n     );

  cache();
}


double CrystalBall::mu()  const
{
  return getPar( 0 ).value();
}


double CrystalBall::sigma() const
{
  return getPar( 1 ).value();
}


double CrystalBall::alpha() const
{
  return getPar( 2 ).value();
}


double CrystalBall::n() const
{
  return getPar( 3 ).value();
}


double CrystalBall::evaluate() const throw( PdfException )
{
  return evaluate( getVar( 0 ).value() );
}


void CrystalBall::cache()
{
  _normCore = sigma() * std::sqrt( 2.0 * M_PI ) * ( 1. + Math::erf( alpha() / std::sqrt( 2.0 ) ) ) / 2.0;
  _normTail = sigma() / alpha() * n() / ( n() - 1.0 ) * std::exp( - std::pow( alpha(), 2 ) / 2.0 );

  _norm = _normCore + _normTail;
}


const double CrystalBall::core( const double& x ) const
{
  return std::exp( - std::pow( x - mu(), 2 ) / ( 2.0 * std::pow( sigma(), 2 ) ) ) / _norm;
}


const double CrystalBall::tail( const double& x ) const
{
  const double& vmu    = mu   ();
  const double& vsigma = sigma();
  const double& valpha = alpha();
  const double& vn     = n    ();

  const double& alphaSq = std::pow( valpha, 2 );

  const double& num = std::pow( vn, vn ) * std::exp( - alphaSq / 2.0 );
  const double& den = std::pow( vn - alphaSq - valpha * ( x - vmu ) / vsigma, vn );

  return num / den / _norm;
}


double CrystalBall::evaluate( double x ) const throw( PdfException )
{
  if ( x > mu() - alpha() * sigma() )
    return core( x );
  else
    return tail( x );
}


double CrystalBall::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  return evaluate( vars[ 0 ] );
}
