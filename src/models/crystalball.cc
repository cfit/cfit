
#include <cfit/math.hh>
#include <cfit/models/crystalball.hh>

#include <cfit/random.hh>

CrystalBall::CrystalBall( const Variable& x,
			  const Parameter& mu, const Parameter& sigma, const Parameter& alpha, const Parameter& n )
  : _mu( mu ), _sigma( sigma ), _alpha( alpha ), _n( n ), _hasLower( false ), _hasUpper( false ), _lower( 0.0 ), _upper( 0.0 )
{
  push( x );

  push( mu    );
  push( sigma );
  push( alpha );
  push( n     );

  cache();
}


CrystalBall::CrystalBall( const Variable&      x    ,
                          const ParameterExpr& mu   ,
                          const ParameterExpr& sigma,
                          const ParameterExpr& alpha,
                          const ParameterExpr& n     )
  : _mu( mu ), _sigma( sigma ), _alpha( alpha ), _n( n ), _hasLower( false ), _hasUpper( false ), _lower( 0.0 ), _upper( 0.0 )
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
  return _mu.evaluate();
}


double CrystalBall::sigma() const
{
  return _sigma.evaluate();
}


double CrystalBall::alpha() const
{
  return _alpha.evaluate();
}


double CrystalBall::n() const
{
  return _n.evaluate();
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

  const double& chi = ( x - vmu ) / vsigma;

  const double& sqrt2    = std::sqrt( 2.0 );
  const double& sqrtpih  = std::sqrt( M_PI / 2.0 );

  if ( valpha > 0 )
  {
    if ( chi <= - valpha )
    {
      // Lower tail piece.
      const double& term1 = vsigma / valpha * vn / ( vn - 1.0 ) * std::exp( - alphaSq / 2.0 );
      const double& term2 = std::pow( vn / ( vn - alphaSq - valpha * chi ), vn - 1.0 );
      return term1 * term2;
    }
    else
    {
      // Complete area of the tail.
      const double& normTail = vsigma / valpha * vn / ( vn - 1.0 ) * std::exp( - alphaSq / 2.0 );

      // Calculation of the area under the core (Gaussian) piece.
      return normTail + vsigma * sqrtpih * ( Math::erf( chi / sqrt2 ) + Math::erf( valpha / sqrt2 ) );
    }
  }
  else // If alpha < 0.
  {
    if ( chi > - valpha )
    {
      const double& normCore = vsigma * sqrtpih * ( 1.0 + Math::erf( - valpha / sqrt2 ) );

      const double& term1 = vsigma / std::fabs( valpha ) * vn / ( vn - 1.0 ) * std::exp( - alphaSq / 2.0 );
      const double& term2 = 1.0 - std::pow( vn / ( vn - alphaSq - valpha * chi ), vn - 1.0 );

      return normCore + term1 * term2;
    }
    else
    {
      return vsigma * sqrtpih * ( 1.0 + Math::erf( chi / sqrt2 ) );
    }
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
    const double& absAlpha = std::fabs( alpha() );

    const double& normCore = sigma() * std::sqrt( 2.0 * M_PI ) * ( 1.0 + Math::erf( absAlpha / std::sqrt( 2.0 ) ) ) / 2.0;
    const double& normTail = sigma() / absAlpha * n() / ( n() - 1.0 ) * std::exp( - std::pow( alpha(), 2 ) / 2.0 );
    areaUp = normCore + normTail;
  }

  // Assign the value of the norm as the difference of areas between the upper and lower limits.
  _norm = areaUp - areaLo;
}


// Function to compute the value of the pdf at its core (Gaussian) part.
const double CrystalBall::core( const double& chi ) const
{
  return std::exp( - std::pow( chi, 2 ) / 2.0 ) / _norm;
}


// Function to compute the value of the pdf at its tail.
const double CrystalBall::tail( const double& chi ) const
{
  const double& valpha = alpha();
  const double& vn     = n    ();

  const double& alphaSq = std::pow( valpha, 2 );

  const double& term1 = std::exp( - alphaSq / 2.0 );
  const double& term2 = std::pow( vn / ( vn - alphaSq - std::fabs( valpha ) * chi ), vn );

  return term1 * term2 / _norm;
}


// Evaluate the function at a given point.
const double CrystalBall::evaluate( const double& x ) const throw( PdfException )
{
  if ( _hasLower && ( x < _lower ) )
    return 0.0;

  if ( _hasUpper && ( x > _upper ) )
    return 0.0;

  const double& sign = ( alpha() > 0 ) - ( alpha() < 0 );
  const double& chi  = sign * ( x - mu() ) / sigma();

  if ( chi < - std::fabs( alpha() ) )
    return tail( chi );

  return core( chi );
}


const double CrystalBall::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  return evaluate( vars[ 0 ] );
}


void CrystalBall::setParExpr()
{
  _mu   .setPars( _parMap );
  _sigma.setPars( _parMap );
  _alpha.setPars( _parMap );
  _n    .setPars( _parMap );
}


const double CrystalBall::area( const double& min, const double& max ) const throw( PdfException )
{
  const double& xmin = _hasLower ? std::max( min, _lower ) : min;
  const double& xmax = _hasUpper ? std::min( max, _upper ) : max;

  return ( cumulativeNorm( xmax ) - cumulativeNorm( xmin ) ) / _norm;
}


const std::map< std::string, double > CrystalBall::generate() const throw( PdfException )
{
  const double& vmu    = mu();
  const double& vn     = n();
  const double& vnm1   = vn - 1.0;
  const double& vsigma = sigma();
  const double& valpha = alpha();
  const double& sqrt2  = std::sqrt( 2.0 );

  // Evaluate the area up to the lower limit (0 if it's -infinity).
  double areaLo = 0.0;
  if ( valpha > 0 && _hasLower ) areaLo = cumulativeNorm(   _lower ) / _norm;
  if ( valpha < 0 && _hasUpper )
  {
    const double& absAlpha = std::fabs( alpha() );

    const double& normCore  = sigma() * std::sqrt( 2.0 * M_PI ) * ( 1.0 + Math::erf( absAlpha / std::sqrt( 2.0 ) ) ) / 2.0;
    const double& normTail  = sigma() / absAlpha * n() / ( n() - 1.0 ) * std::exp( - std::pow( valpha, 2 ) / 2.0 );
    const double& normTotal = normCore + normTail;
    areaLo = ( normTotal - cumulativeNorm( _upper ) ) / _norm;
  }

  // Cumulative pdf evaluated at the tail-core threshold for alpha > 0.
  const double& areaTh = vsigma / ( std::fabs( valpha ) * _norm ) * vn / vnm1 * std::exp( - std::pow( valpha, 2 ) / 2.0 );

  // Generate a flat random number.
  std::uniform_real_distribution< double > dist( 0.0, 1.0 );
  const double& unif = dist( Random::engine() ) + areaLo;

  // If random number is below the tail-core threshold, generate the tail.
  //    Otherwise, generate the (Gaussian) core.
  double genChi;
  if ( unif < areaTh )
    genChi = - std::fabs( valpha ) - vn / std::fabs( valpha ) * ( std::pow( areaTh / unif, 1.0 / vnm1 ) - 1.0 );
  else
    genChi = sqrt2 * Math::inverf( _norm / vsigma * sqrt( 2.0 / M_PI ) * ( unif - areaTh ) - std::erf( std::fabs( valpha ) / sqrt2 ) );

  genChi *= ( valpha > 0 ) - ( valpha < 0 );

  double genVal = vmu + vsigma * genChi;

  std::map< std::string, double > gen;
  gen[ getVar( 0 ).name() ] = genVal;

  return gen;
}

