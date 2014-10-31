#ifndef __COEF_HH__
#define __COEF_HH__

#include <complex>

#include <cfit/parameter.hh>

class Amplitude;
class Resonance;
class Fvector;
class CoefExpr;

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

  void setValue( const double& re, const double& im )
  {
    _real.setValue( re );
    _imag.setValue( im );
  };

  const std::complex< double > value()   const;
  const bool                   isFixed() const;

  // Operations of coefficients with themselves.
  friend const CoefExpr operator+( const Coef&   left, const Coef&   right );
  friend const CoefExpr operator-( const Coef&   left, const Coef&   right );
  friend const CoefExpr operator*( const Coef&   left, const Coef&   right );
  friend const CoefExpr operator/( const Coef&   left, const Coef&   right );

  // Operations of constants and coefficients.
  friend const CoefExpr operator+( const double& left, const Coef&   right );
  friend const CoefExpr operator-( const double& left, const Coef&   right );
  friend const CoefExpr operator*( const double& left, const Coef&   right );
  friend const CoefExpr operator/( const double& left, const Coef&   right );

  friend const CoefExpr operator+( const Coef&   left, const double& right );
  friend const CoefExpr operator-( const Coef&   left, const double& right );
  friend const CoefExpr operator*( const Coef&   left, const double& right );
  friend const CoefExpr operator/( const Coef&   left, const double& right );


  // Binary coefficient operations.
  friend const CoefExpr pow      ( const Coef&     left, const std::complex< double >& right );
  friend const CoefExpr pow      ( const CoefExpr& left, const std::complex< double >& right );

  // Unary minus.
  friend const CoefExpr operator-( const Coef&     coef );
  friend const CoefExpr operator-( const CoefExpr& coef );

  // Unary coefficient operations.
  friend const CoefExpr exp      ( const Coef&     coef );
  friend const CoefExpr log      ( const Coef&     coef );
  friend const CoefExpr sin      ( const Coef&     coef );
  friend const CoefExpr cos      ( const Coef&     coef );
  friend const CoefExpr tan      ( const Coef&     coef );

  friend const CoefExpr exp      ( const CoefExpr& coef );
  friend const CoefExpr log      ( const CoefExpr& coef );
  friend const CoefExpr sin      ( const CoefExpr& coef );
  friend const CoefExpr cos      ( const CoefExpr& coef );
  friend const CoefExpr tan      ( const CoefExpr& coef );

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

  // Binary arithmetic operations with resonances and coeffiecient expressions.
  friend const Amplitude operator*( const CoefExpr&  left, const Resonance& right );
  friend const Amplitude operator*( const Resonance& left, const CoefExpr&  right );
  friend const Amplitude operator/( const Resonance& left, const CoefExpr&  right );

  friend const Amplitude operator*( const CoefExpr&  left, const Fvector&   right );
  friend const Amplitude operator*( const Fvector&   left, const CoefExpr&  right );
  friend const Amplitude operator/( const Fvector&   left, const CoefExpr&  right );

  friend const Amplitude operator*( const CoefExpr&  left, const Amplitude& right );
  friend const Amplitude operator*( const Amplitude& left, const CoefExpr&  right );
  friend const Amplitude operator/( const Amplitude& left, const CoefExpr&  right );
};

#endif
