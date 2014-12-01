
#include <cfit/math.hh>
#include <cfit/models/doublecrystalball.hh>

#include <cfit/random.hh>

DoubleCrystalBall::DoubleCrystalBall( const Variable& x,
                                      const Parameter& mu   , const Parameter& sigma,
                                      const Parameter& alpha, const Parameter& n,
                                      const Parameter& beta , const Parameter& m )
  : _hasLower( false ), _hasUpper( false ), _lower( 0.0 ), _upper( 0.0 )
{
  push( x );

  push( mu    );
  push( sigma );
  push( alpha );
  push( n     );
  push( beta  );
  push( m     );

  cache();
}



DoubleCrystalBall* DoubleCrystalBall::copy() const
{
  return new DoubleCrystalBall( *this );
}


void DoubleCrystalBall::setLowerLimit( const double& lower )
{
  _hasLower = true;
  _lower    = lower;

  cache();
}


void DoubleCrystalBall::setUpperLimit( const double& upper )
{
  _hasUpper = true;
  _upper    = upper;

  cache();
}


void DoubleCrystalBall::setLimits( const double& lower, const double& upper )
{
  _hasLower = true;
  _hasUpper = true;
  _lower    = lower;
  _upper    = upper;

  cache();
}


void DoubleCrystalBall::unsetLowerLimit()
{
  _hasLower = false;

  cache();
}


void DoubleCrystalBall::unsetUpperLimit()
{
  _hasUpper = false;

  cache();
}


void DoubleCrystalBall::unsetLimits()
{
  _hasLower = false;
  _hasUpper = false;

  cache();
}


double DoubleCrystalBall::mu()  const
{
  return getPar( 0 ).value();
}


double DoubleCrystalBall::sigma() const
{
  return getPar( 1 ).value();
}


double DoubleCrystalBall::alpha() const
{
  return getPar( 2 ).value();
}


double DoubleCrystalBall::n() const
{
  return getPar( 3 ).value();
}


double DoubleCrystalBall::beta() const
{
  return getPar( 4 ).value();
}


double DoubleCrystalBall::m() const
{
  return getPar( 5 ).value();
}


// Compute the area of the unnormalized pdf up to given value x.
// Kind of an unnormalized cdf.
const double DoubleCrystalBall::cumulativeNorm( const double& x ) const
{
  const double& vmu    = mu   ();
  const double& vsigma = sigma();
  const double& valpha = alpha();
  const double& vn     = n    ();
  const double& vbeta  = beta ();
  const double& vm     = m    ();

  const double& alphaSq = std::pow( valpha, 2 );
  const double& betaSq  = std::pow( vbeta , 2 );

  if ( x < vmu - valpha * vsigma )
  {
    // Lower tail piece.
    const double& num = std::pow( vn, vn ) * std::exp( - alphaSq / 2.0 );
    const double& den = std::pow( vn - alphaSq - valpha * ( x - vmu ) / vsigma, vn - 1.0 );
    return vsigma / valpha / ( vn - 1.0 ) * num / den;
  }
  else if ( x < vmu + vbeta * vsigma )
  {
    // Complete area of the lower tail.
    const double& normTail = vsigma / valpha * vn / ( vn - 1.0 ) * std::exp( - alphaSq / 2.0 );

    // Calculation of the area under the core (Gaussian) piece.
    const double& sqrt2   = std::sqrt( 2.0 );
    const double& sqrtpih = std::sqrt( M_PI / 2.0 );
    const double& erfs    = Math::erf( valpha / sqrt2 ) + Math::erf( ( x - vmu ) / ( vsigma * sqrt2 ) );

    return normTail + vsigma * sqrtpih * erfs;
  }
  else
  {
    const double& sqrt2 = std::sqrt( 2.0 );

    const double& normTail = vsigma / valpha * vn / ( vn - 1.0 ) * std::exp( - alphaSq / 2.0 );
    const double& erfs     = Math::erf( valpha / sqrt2 ) + Math::erf( vbeta / sqrt2 );
    const double& normCore = vsigma * std::sqrt( M_PI / 2.0 ) * erfs;

    const double& term1 = vsigma / vbeta * vm / ( vm - 1.0 ) * std::exp( -betaSq / 2.0 );
    const double& term2 = 1.0 - std::pow( vm / ( vm - betaSq + vbeta * ( x - vmu ) / vsigma ), vm - 1.0 );

    return normTail + normCore + term1 * term2;
  }
}


void DoubleCrystalBall::cache()
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
    const double& vsigma = sigma();
    const double& valpha = alpha();
    const double& vn     = n    ();
    const double& vbeta  = beta ();
    const double& vm     = m    ();

    const double& alphaSq = std::pow( valpha, 2 );
    const double& betaSq  = std::pow( vbeta , 2 );

    const double& sqrt2 = std::sqrt( 2.0 );

    const double& normTailLo = vsigma / valpha * vn / ( vn - 1.0 ) * std::exp( - alphaSq / 2.0 );
    const double& normTailUp = vsigma / vbeta  * vm / ( vm - 1.0 ) * std::exp( - betaSq  / 2.0 );
    const double& erfs       = Math::erf( valpha / sqrt2 ) + Math::erf( vbeta / sqrt2 );
    const double& normCore   = vsigma * std::sqrt( M_PI / 2.0 ) * erfs;

    areaUp = normCore + normTailLo + normTailUp;
  }

  // Assign the value of the norm as the difference of areas between the upper and lower limits.
  _norm = areaUp - areaLo;
}


// Function to compute the value of the pdf at its core (Gaussian) part.
const double DoubleCrystalBall::core( const double& x ) const
{
  return std::exp( - std::pow( x - mu(), 2 ) / ( 2.0 * std::pow( sigma(), 2 ) ) ) / _norm;
}


// Function to compute the value of the pdf at its lower tail.
const double DoubleCrystalBall::tailLo( const double& x ) const
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


// Function to compute the value of the pdf at its upper tail.
const double DoubleCrystalBall::tailUp( const double& x ) const
{
  const double& vmu    = mu   ();
  const double& vsigma = sigma();
  const double& vbeta  = beta ();
  const double& vm     = m    ();

  const double& betaSq = std::pow( vbeta, 2 );

  const double& num = std::pow( vm, vm ) * std::exp( - betaSq / 2.0 );
  const double& den = std::pow( vm - betaSq + vbeta * ( x - vmu ) / vsigma, vm );

  return num / den / _norm;
}


// Evaluate the function at a given point.
const double DoubleCrystalBall::evaluate( const double& x ) const throw( PdfException )
{
  if ( _hasLower && ( x < _lower ) )
    return 0.0;

  if ( _hasUpper && ( x > _upper ) )
    return 0.0;

  if ( x < mu() - alpha() * sigma() )
    return tailLo( x );

  if ( x > mu() + beta() * sigma() )
    return tailUp( x );

  return core( x );
}


const double DoubleCrystalBall::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  return evaluate( vars[ 0 ] );
}


const double DoubleCrystalBall::area( const double& min, const double& max ) const throw( PdfException )
{
  return ( cumulativeNorm( max ) - cumulativeNorm( min ) ) / _norm;
}



const std::map< std::string, double > DoubleCrystalBall::generate() const throw( PdfException )
{
  const double& vn     = n();
  const double& vnm1   = vn - 1.0;
  const double& vm     = m();
  const double& vmm1   = vm - 1.0;
  const double& vsigma = sigma();
  const double& valpha = alpha();
  const double& vbeta  = beta();
  const double& sqrt2  = std::sqrt( 2.0 );

  // Evaluate the area up to the lower limit (0 if it's -infinity).
  double areaLo = 0.0;
  if ( _hasLower )
    areaLo = cumulativeNorm( _lower ) / _norm;

  // Cumulative pdf evaluated at the tail-core threshold.
  const double& areaTh1 = vsigma / ( valpha * _norm ) * vn / vnm1 * std::exp( - std::pow( valpha, 2 ) / 2.0 );
  const double& areaTh2 = areaTh1 + vsigma / _norm * std::sqrt( M_PI / 2.0 ) * ( std::erf( valpha / sqrt2 ) + std::erf( vbeta / sqrt2 ) );

  // Useful factor used in the generation of numbers in the upper tail.
  const double& coefHi = vsigma / ( vbeta * _norm ) * vm / vmm1 * std::exp( - std::pow( vbeta, 2 ) / 2.0 );

  // Generate a flat random number.
  std::uniform_real_distribution< double > dist( 0.0, 1.0 );
  const double& unif = dist( Random::engine() ) + areaLo;

  double genVal = 0.0;

  // If random number is below the lower tail-core threshold, generate the lower tail.
  //    If it's above the upper tail-core threshold, generate the upper tail.
  //    Otherwise, generate the (Gaussian) core.
  if ( unif < areaTh1 )
    genVal = mu() - valpha * vsigma - vn * vsigma / valpha * ( std::pow( areaTh1 / unif                      , 1.0 / vnm1 ) - 1.0 );
  else if ( unif > areaTh2 )
    genVal = mu() + vbeta  * vsigma + vm * vsigma / vbeta  * ( std::pow( coefHi / ( areaTh2 + coefHi - unif ), 1.0 / vmm1 ) - 1.0 );
  else
    genVal = mu() + vsigma * sqrt2 * Math::inverf( _norm / vsigma * std::sqrt( 2.0 / M_PI ) * ( unif - areaTh1 ) -
                                                   std::erf( valpha / sqrt2 ) );

  // Fill the map with the generated value.
  std::map< std::string, double > gen;
  gen[ getVar( 0 ).name() ] = genVal;

  return gen;
}


