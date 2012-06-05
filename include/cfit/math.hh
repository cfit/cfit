#ifndef __MATH_HH__
#define __MATH_HH__

class Math
{
private:
  static const double lngamma(       double  x );
public:
  static const double erf    ( const double& x );
  static const double erfc   ( const double& x );
  static const double gamma_p( const double& a, const double& x );
  static const double gamma_q( const double& a, const double& x );
};


#endif
