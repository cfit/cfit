
#include <cmath>
#include <cfit/math.hh>

const double Math::erf( const double& x )
{
  // # save the sign of x
  // sign = 1 if x >= 0 else -1
  // x = abs(x)

  // constants
  const double& a1 =  0.254829592;
  const double& a2 = -0.284496736;
  const double& a3 =  1.421413741;
  const double& a4 = -1.453152027;
  const double& a5 =  1.061405429;
  const double& p  =  0.3275911;

  // A&S formula 7.1.26
  const double& t = 1.0 / ( 1.0 + p * std::abs( x ) );
  double approx = 1.0 - ( ( ( ( ( a5 * t + a4 ) * t ) + a3 ) * t + a2 ) * t + a1 ) * t * std::exp( -std::pow( x, 2 ) );

  if ( x < 0 )
    return -approx;

  return approx;
}

