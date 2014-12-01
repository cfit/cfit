
#include <cfit/models/expogauss.hh>
#include <cfit/math.hh>

#include <cfit/random.hh>

ExpoGauss::ExpoGauss( const Variable&  x ,
                      const Parameter& gamma, const Parameter& mu, const Parameter& sigma )
  : _hasLower( false ), _hasUpper( false ), _lower( 0.0 ), _upper( 0.0 ),
    _normExpo( 0.0 ), _normGauss( 0.0 ), _norm( 0.0 )
{
  push( x );

  push( gamma );
  push( mu    );
  push( sigma );

  cache();
}


ExpoGauss* ExpoGauss::copy() const
{
  return new ExpoGauss( *this );
}


double ExpoGauss::gamma()  const
{
  return getPar( 0 ).value();
}


double ExpoGauss::mu() const
{
  return getPar( 1 ).value();
}


double ExpoGauss::sigma() const
{
  return getPar( 2 ).value();
}



void ExpoGauss::setLowerLimit( const double& lower )
{
  _hasLower = true;
  _lower    = lower;

  cache();
}


void ExpoGauss::setUpperLimit( const double& upper )
{
  _hasUpper = true;
  _upper    = upper;

  cache();
}


void ExpoGauss::setLimits( const double& lower, const double& upper )
{
  _hasLower = true;
  _hasUpper = true;
  _lower    = lower;
  _upper    = upper;

  cache();
}


void ExpoGauss::unsetLowerLimit()
{
  _hasLower = false;

  cache();
}


void ExpoGauss::unsetUpperLimit()
{
  _hasUpper = false;

  cache();
}


void ExpoGauss::unsetLimits()
{
  _hasLower = false;
  _hasUpper = false;

  cache();
}


void ExpoGauss::cache()
{
  // Evaluate the norm of the exponential.
  _normExpo = 1.0 / gamma();

  // Evaluate the norm of the Gaussian.
  _normGauss = sigma() * std::sqrt( 2.0 * M_PI );

  // If no lower or upper limit has been set, the pdf is already normalized to 1.
  if ( ( ! _hasLower ) && ( ! _hasUpper ) )
  {
    _norm = 1.0;
    return;
  }

  // Fast references to mu and sigma gaussian parameters.
  const double& vmu    = mu();
  const double& vsigma = sigma();

  // If the pdf has limits, define the integration technique.
  const unsigned& nbin = 300;

  // Approximate minus infinity to ( mu - 5 sigma ), which is 5 sigma away from the riseup of the pdf.
  const double& minusInf = vmu - 5.0 * vsigma;

  // If an upper limit has been set, integrate from the lower limit and assume it's
  //    mu - 5 sigma if it has not been set.
  if ( _hasUpper )
  {
    const double& lower = _hasLower ? _lower : minusInf;
    const double& dx    = ( _upper - lower ) / double( nbin );

    _norm = 0.0;
    for ( double x = lower; x < _upper; x += dx )
      _norm += expogauss( x );
    _norm *= dx;

    return;
  }

  // If the pdf has a lower limit but not an upper limit, normalize by integrating from
  //    5 sigmas from the left and subtract it from 1.
  if ( _hasLower && ( ! _hasUpper ) )
  {
    // If the lower limit is very far away from the riseup of the pdf, the norm is 1.
    if ( _lower < minusInf )
    {
      _norm = 1.0;
      return;
    }

    // Otherwise, integrate from minus infinity to the lower limit.
    const double& dx = ( _lower - minusInf ) / double( nbin );
    _norm = 0.0;
    for ( double x = minusInf; x < _lower; x += dx )
      _norm += expogauss( x );
    _norm = 1.0 - _norm * dx;

    return;
  }
}


const double ExpoGauss::expogauss( const double& x ) const
{
  if ( _hasLower && ( x < _lower ) )
    return 0.0;
  if ( _hasUpper && ( x > _upper ) )
    return 0.0;

  const double& vgamma = gamma();
  const double& vmu    = mu();
  const double& vsigma = sigma();
  const double& sqrt2  = std::sqrt( 2.0 );

  const double& vgammaSq = std::pow( vgamma, 2 );
  const double& vsigmaSq = std::pow( vsigma, 2 );

  const double& expo    = std::exp( - vgamma * ( x - vmu ) + vgammaSq * vsigmaSq / 2.0 );
  const double& erfFrac = ( 1.0 + std::erf( ( x - vmu - vgamma * vsigmaSq ) / ( vsigma * sqrt2 ) ) ) / 2.0;

  return expo * erfFrac / _normExpo;
}


const double ExpoGauss::evaluate( const double& x ) const throw( PdfException )
{
  return expogauss( x ) / _norm;
}


const double ExpoGauss::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  return evaluate( vars[ 0 ] );
}


const double ExpoGauss::area( const double& min, const double& max ) const throw( PdfException )
{
  // Set the limits of integration.
  const double& xmin = _hasLower ? std::max( min, _lower ) : min;
  const double& xmax = _hasUpper ? std::min( max, _upper ) : max;

  // Define the precision.
  const unsigned& nbin = 300;
  const double& dx = ( xmax - xmin ) / double( nbin );

  // Do the integral.
  double retval = 0.0;
  for ( double x = xmin; x < xmax; x += dx )
    retval += evaluate( x );

  return retval * dx;
}


const std::map< std::string, double > ExpoGauss::generate() const throw( PdfException )
{
  double x = 0.0;

  bool withinLimits = false;
  while ( ! withinLimits )
  {
    // Add an Argus distributed variable and a Gaussian variable.
    x = -std::log( Random::flat() ) / gamma() + Random::normal( mu(), sigma() );

    // Check if the result is within limits, if any.
    withinLimits = true;
    withinLimits &= ( ! _hasLower ) || ( x > _lower );
    withinLimits &= ( ! _hasUpper ) || ( x < _upper );
  }

  // Fill a map with the generated value.
  std::map< std::string, double > gen;
  gen[ getVar( 0 ).name() ] = x;

  return gen;
}
