#ifndef __COEF_HH__
#define __COEF_HH__

#include <complex>

#include <cfit/parameter.hh>

class Amplitude;
class Resonance;

class Coef
{
private:
  Parameter _real;
  Parameter _imag;
public:
  Coef( const Parameter& real, const Parameter& imag )
    : _real( real ), _imag( imag )
  {}
  const Parameter& real() const { return _real; }
  const Parameter& imag() const { return _imag; }

  const std::complex< double > value() const;

  // Operations of coefficients with themselves.
  friend const Amplitude operator+( const Coef&      left, const Coef&      right );
  friend const Amplitude operator-( const Coef&      left, const Coef&      right );
  friend const Amplitude operator*( const Coef&      left, const Coef&      right );
  friend const Amplitude operator/( const Coef&      left, const Coef&      right );

  // Operations of constants and coefficients.
  friend const Amplitude operator+( const double&    left, const Coef&      right );
  friend const Amplitude operator-( const double&    left, const Coef&      right );
  friend const Amplitude operator*( const double&    left, const Coef&      right );
  friend const Amplitude operator/( const double&    left, const Coef&      right );

  friend const Amplitude operator+( const Coef&      left, const double&    right );
  friend const Amplitude operator-( const Coef&      left, const double&    right );
  friend const Amplitude operator*( const Coef&      left, const double&    right );
  friend const Amplitude operator/( const Coef&      left, const double&    right );

  // Operations of coefficients and resonances.
  friend const Amplitude operator+( const Coef&      left, const Resonance& right );
  friend const Amplitude operator-( const Coef&      left, const Resonance& right );
  friend const Amplitude operator*( const Coef&      left, const Resonance& right );

  friend const Amplitude operator+( const Resonance& left, const Coef&      right );
  friend const Amplitude operator-( const Resonance& left, const Coef&      right );
  friend const Amplitude operator*( const Resonance& left, const Coef&      right );
  friend const Amplitude operator/( const Resonance& left, const Coef&      right );

  // Operations with coefficients and amplitudes.
  friend const Amplitude operator+( const Coef&      left, const Amplitude& right );
  friend const Amplitude operator-( const Coef&      left, const Amplitude& right );
  friend const Amplitude operator*( const Coef&      left, const Amplitude& right );

  friend const Amplitude operator+( const Amplitude& left, const Coef&      right );
  friend const Amplitude operator-( const Amplitude& left, const Coef&      right );
  friend const Amplitude operator*( const Amplitude& left, const Coef&      right );
  friend const Amplitude operator/( const Amplitude& left, const Coef&      right );
};

#endif
