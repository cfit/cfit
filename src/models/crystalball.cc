
#include <cfit/math.hh>
#include <cfit/models/crystalball.hh>

CrystalBall::CrystalBall( const Variable& x,
			  const Parameter& mu, const Parameter& sigma, const Parameter& alpha, const Parameter& n )
  : _hasLower( false ), _hasUpper( false ), _lower( 0.0 ), _upper( 0.0 )
{
  push( x );

  push( mu    );
  push( sigma );
  push( alpha );
  push( n     );

  cache();
}


CrystalBall* CrystalBall::copy() const
{
  return new CrystalBall( *this );
}


void CrystalBall::setLowerLimit( const double& lower )
{
  _hasLower = true;
  _lower    = lower;

  cache();
}


void CrystalBall::setUpperLimit( const double& upper )
{
  _hasUpper = true;
  _upper    = upper;

  cache();
}


void CrystalBall::setLimits( const double& lower, const double& upper )
{
  _hasLower = true;
  _hasUpper = true;
  _lower    = lower;
  _upper    = upper;

  cache();
}


void CrystalBall::unsetLowerLimit()
{
  _hasLower = false;

  cache();
}


void CrystalBall::unsetUpperLimit()
{
  _hasUpper = false;

  cache();
}


void CrystalBall::unsetLimits()
{
  _hasLower = false;
  _hasUpper = false;

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


// Compute the area of the unnormalized pdf up to given value x.
// Kind of an unnormalized cdf.
const double CrystalBall::cumulativeNorm( const double& x ) const
{
  const double& vmu    = mu   ();
  const double& vsigma = sigma();
  const double& valpha = alpha();
  const double& vn     = n    ();

  const double& alphaSq = std::pow( valpha, 2 );

  if ( x < vmu - valpha * vsigma )
  {
    // Tail piece.
    const double& num = std::pow( vn, vn ) * std::exp( - alphaSq / 2.0 );
    const double& den = std::pow( vn - alphaSq - valpha * ( x - vmu ) / vsigma, vn - 1.0 );
    return vsigma / valpha / ( vn - 1.0 ) * num / den;
  }
  else
  {
    // Complete area of the tail.
    const double& normTail = vsigma / valpha * vn / ( vn - 1.0 ) * std::exp( - alphaSq / 2.0 );

    // Calculation of the area under the core (Gaussian) piece.
    const double& sqrt2    = std::sqrt( 2.0 );
    const double& sqrtpih  = std::sqrt( M_PI / 2.0 );
    return normTail + vsigma * sqrtpih * ( Math::erf( ( x - vmu ) / ( vsigma * sqrt2 ) ) + Math::erf( valpha / sqrt2 ) );
  }
}


void CrystalBall::cache()
{
  // Evaluate the area up to the lower limit (0 if it's -infinity).
  double areaLo = 0.0;
  if ( _hasLower )
    areaLo = cumulativeNorm( _lower );

  // Evaluate the area up to the upper limit (complete area if it's +infinity).
  double areaUp = 0.0;
  if ( _hasUpper )
    areaUp = cumulativeNorm( _upper );
  else
  {
    const double& normCore = sigma() * std::sqrt( 2.0 * M_PI ) * ( 1. + Math::erf( alpha() / std::sqrt( 2.0 ) ) ) / 2.0;
    const double& normTail = sigma() / alpha() * n() / ( n() - 1.0 ) * std::exp( - std::pow( alpha(), 2 ) / 2.0 );
    areaUp = normCore + normTail;
  }

  // Assign the value of the norm as the difference of areas between the upper and lower limits.
  _norm = areaUp - areaLo;
}


// Function to compute the value of the pdf at its core (Gaussian) part.
const double CrystalBall::core( const double& x ) const
{
  return std::exp( - std::pow( x - mu(), 2 ) / ( 2.0 * std::pow( sigma(), 2 ) ) ) / _norm;
}


// Function to compute the value of the pdf at its tail.
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


// Evaluate the function at a given point.
const double CrystalBall::evaluate( const double& x ) const throw( PdfException )
{
  if ( x > mu() - alpha() * sigma() )
    return core( x );
  else
    return tail( x );
}


const double CrystalBall::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  return evaluate( vars[ 0 ] );
}


const double CrystalBall::evaluate() const throw( PdfException )
{
  return evaluate( getVar( 0 ).value() );
}


const double CrystalBall::area( const double& min, const double& max ) const throw( PdfException )
{
  return ( cumulativeNorm( std::min( max, _upper ) ) - cumulativeNorm( std::max( min, _lower ) ) ) / _norm;
}
