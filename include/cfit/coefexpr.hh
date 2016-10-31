#ifndef __COEFEXPR_HH__
#define __COEFEXPR_HH__

#include <complex>
#include <vector>
#include <map>
#include <string>

#include <cfit/parameter.hh>
#include <cfit/parameterexpr.hh>
#include <cfit/coef.hh>
#include <cfit/operation.hh>


class Resonance;
class Fvector;
class Amplitude;

class CoefExpr
{
  friend class Amplitude;

private:
  std::vector< std::complex< double > > _ctnts;
  std::vector< Parameter >              _parms;
  std::vector< Coef >                   _coefs;
  std::string                           _expression;
  std::vector< Operation::Op >          _opers;

  void append( const std::complex< double >& val  );
  void append( const double&                 val  );
  void append( const Parameter&              par  );
  void append( const ParameterExpr&          expr );
  void append( const Coef&                   coef );
  void append( const CoefExpr&               expr );

  template< class L, class R >
  CoefExpr( const L& left, const R& right, const Operation::Op& oper )
  {
    append( left  );
    append( right );

    _expression += "b"; // b = binary operation.
    _opers.push_back( oper );
  }

  template< class T >
  CoefExpr( const T& var, const Operation::Op& oper )
  {
    append( var );

    _expression += "u"; // u = unary operation.
    _opers.push_back( oper );
  }

public:
  CoefExpr() {};

  template <class T>
  CoefExpr( const T& expr ) { append( expr ); }

  const std::map< std::string, Parameter > getPars() const;

  void setPars( const std::map< std::string, Parameter >& pars );
  void setPars( const std::map< std::string, double    >& pars );

  const std::complex< double > evaluate() const throw( PdfException );

  // Binary arithmetic operators.
  friend const CoefExpr operator+( const Coef&                   left, const Coef&                   right );
  friend const CoefExpr operator-( const Coef&                   left, const Coef&                   right );
  friend const CoefExpr operator*( const Coef&                   left, const Coef&                   right );
  friend const CoefExpr operator/( const Coef&                   left, const Coef&                   right );

  friend const CoefExpr operator+( const std::complex< double >& left, const Coef&                   right );
  friend const CoefExpr operator-( const std::complex< double >& left, const Coef&                   right );
  friend const CoefExpr operator*( const std::complex< double >& left, const Coef&                   right );
  friend const CoefExpr operator/( const std::complex< double >& left, const Coef&                   right );

  friend const CoefExpr operator+( const Coef&                   left, const std::complex< double >& right );
  friend const CoefExpr operator-( const Coef&                   left, const std::complex< double >& right );
  friend const CoefExpr operator*( const Coef&                   left, const std::complex< double >& right );
  friend const CoefExpr operator/( const Coef&                   left, const std::complex< double >& right );

  friend const CoefExpr operator+( const CoefExpr&               left, const CoefExpr&               right );
  friend const CoefExpr operator-( const CoefExpr&               left, const CoefExpr&               right );
  friend const CoefExpr operator*( const CoefExpr&               left, const CoefExpr&               right );
  friend const CoefExpr operator/( const CoefExpr&               left, const CoefExpr&               right );

  friend const CoefExpr operator+( const std::complex< double >& left, const CoefExpr&               right );
  friend const CoefExpr operator-( const std::complex< double >& left, const CoefExpr&               right );
  friend const CoefExpr operator*( const std::complex< double >& left, const CoefExpr&               right );
  friend const CoefExpr operator/( const std::complex< double >& left, const CoefExpr&               right );

  friend const CoefExpr operator+( const CoefExpr&               left, const std::complex< double >& right );
  friend const CoefExpr operator-( const CoefExpr&               left, const std::complex< double >& right );
  friend const CoefExpr operator*( const CoefExpr&               left, const std::complex< double >& right );
  friend const CoefExpr operator/( const CoefExpr&               left, const std::complex< double >& right );

  friend const CoefExpr operator+( const Coef&                   left, const CoefExpr&               right );
  friend const CoefExpr operator-( const Coef&                   left, const CoefExpr&               right );
  friend const CoefExpr operator*( const Coef&                   left, const CoefExpr&               right );
  friend const CoefExpr operator/( const Coef&                   left, const CoefExpr&               right );

  friend const CoefExpr operator+( const CoefExpr&               left, const Coef&                   right );
  friend const CoefExpr operator-( const CoefExpr&               left, const Coef&                   right );
  friend const CoefExpr operator*( const CoefExpr&               left, const Coef&                   right );
  friend const CoefExpr operator/( const CoefExpr&               left, const Coef&                   right );


  // Operations of parameters and parameter expressions with complex numbers yield coefficient expressions.
  friend const CoefExpr operator+( const std::complex< double >& left, const Parameter&              right );
  friend const CoefExpr operator-( const std::complex< double >& left, const Parameter&              right );
  friend const CoefExpr operator*( const std::complex< double >& left, const Parameter&              right );
  friend const CoefExpr operator/( const std::complex< double >& left, const Parameter&              right );

  friend const CoefExpr operator+( const Parameter&              left, const std::complex< double >& right );
  friend const CoefExpr operator-( const Parameter&              left, const std::complex< double >& right );
  friend const CoefExpr operator*( const Parameter&              left, const std::complex< double >& right );
  friend const CoefExpr operator/( const Parameter&              left, const std::complex< double >& right );

  friend const CoefExpr operator+( const std::complex< double >& left, const ParameterExpr&          right );
  friend const CoefExpr operator-( const std::complex< double >& left, const ParameterExpr&          right );
  friend const CoefExpr operator*( const std::complex< double >& left, const ParameterExpr&          right );
  friend const CoefExpr operator/( const std::complex< double >& left, const ParameterExpr&          right );

  friend const CoefExpr operator+( const ParameterExpr&          left, const std::complex< double >& right );
  friend const CoefExpr operator-( const ParameterExpr&          left, const std::complex< double >& right );
  friend const CoefExpr operator*( const ParameterExpr&          left, const std::complex< double >& right );
  friend const CoefExpr operator/( const ParameterExpr&          left, const std::complex< double >& right );




  // Operations with parameters and coefficients.
  friend const CoefExpr operator+( const Parameter&     left, const Coef&          right );
  friend const CoefExpr operator-( const Parameter&     left, const Coef&          right );
  friend const CoefExpr operator*( const Parameter&     left, const Coef&          right );
  friend const CoefExpr operator/( const Parameter&     left, const Coef&          right );

  friend const CoefExpr operator+( const Coef&          left, const Parameter&     right );
  friend const CoefExpr operator-( const Coef&          left, const Parameter&     right );
  friend const CoefExpr operator*( const Coef&          left, const Parameter&     right );
  friend const CoefExpr operator/( const Coef&          left, const Parameter&     right );

  // Operations with parameter expressions and coefficients.
  friend const CoefExpr operator+( const ParameterExpr& left, const Coef&          right );
  friend const CoefExpr operator-( const ParameterExpr& left, const Coef&          right );
  friend const CoefExpr operator*( const ParameterExpr& left, const Coef&          right );
  friend const CoefExpr operator/( const ParameterExpr& left, const Coef&          right );

  friend const CoefExpr operator+( const Coef&          left, const ParameterExpr& right );
  friend const CoefExpr operator-( const Coef&          left, const ParameterExpr& right );
  friend const CoefExpr operator*( const Coef&          left, const ParameterExpr& right );
  friend const CoefExpr operator/( const Coef&          left, const ParameterExpr& right );

  // Operations with parameters and coefficient expressions
  friend const CoefExpr operator+( const Parameter&     left, const CoefExpr&      right );
  friend const CoefExpr operator-( const Parameter&     left, const CoefExpr&      right );
  friend const CoefExpr operator*( const Parameter&     left, const CoefExpr&      right );
  friend const CoefExpr operator/( const Parameter&     left, const CoefExpr&      right );

  friend const CoefExpr operator+( const CoefExpr&      left, const Parameter&     right );
  friend const CoefExpr operator-( const CoefExpr&      left, const Parameter&     right );
  friend const CoefExpr operator*( const CoefExpr&      left, const Parameter&     right );
  friend const CoefExpr operator/( const CoefExpr&      left, const Parameter&     right );

  // Operations with parameter expressions and coefficient expressions.
  friend const CoefExpr operator+( const ParameterExpr& left, const CoefExpr&      right );
  friend const CoefExpr operator-( const ParameterExpr& left, const CoefExpr&      right );
  friend const CoefExpr operator*( const ParameterExpr& left, const CoefExpr&      right );
  friend const CoefExpr operator/( const ParameterExpr& left, const CoefExpr&      right );

  friend const CoefExpr operator+( const CoefExpr&      left, const ParameterExpr& right );
  friend const CoefExpr operator-( const CoefExpr&      left, const ParameterExpr& right );
  friend const CoefExpr operator*( const CoefExpr&      left, const ParameterExpr& right );
  friend const CoefExpr operator/( const CoefExpr&      left, const ParameterExpr& right );

  // Operations of constants and coefficients.
  friend const CoefExpr operator+( const double&    left, const Coef&      right );
  friend const CoefExpr operator-( const double&    left, const Coef&      right );
  friend const CoefExpr operator*( const double&    left, const Coef&      right );
  friend const CoefExpr operator/( const double&    left, const Coef&      right );

  friend const CoefExpr operator+( const Coef&      left, const double&    right );
  friend const CoefExpr operator-( const Coef&      left, const double&    right );
  friend const CoefExpr operator*( const Coef&      left, const double&    right );
  friend const CoefExpr operator/( const Coef&      left, const double&    right );

  // Operations of constants and coefficients.
  friend const CoefExpr operator+( const double&    left, const CoefExpr&  right );
  friend const CoefExpr operator-( const double&    left, const CoefExpr&  right );
  friend const CoefExpr operator*( const double&    left, const CoefExpr&  right );
  friend const CoefExpr operator/( const double&    left, const CoefExpr&  right );

  friend const CoefExpr operator+( const CoefExpr&  left, const double&    right );
  friend const CoefExpr operator-( const CoefExpr&  left, const double&    right );
  friend const CoefExpr operator*( const CoefExpr&  left, const double&    right );
  friend const CoefExpr operator/( const CoefExpr&  left, const double&    right );



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
  friend const CoefExpr tanh     ( const Coef&     coef );
  friend const CoefExpr atanh    ( const Coef&     coef );

  friend const CoefExpr exp      ( const CoefExpr& coef );
  friend const CoefExpr log      ( const CoefExpr& coef );
  friend const CoefExpr sin      ( const CoefExpr& coef );
  friend const CoefExpr cos      ( const CoefExpr& coef );
  friend const CoefExpr tan      ( const CoefExpr& coef );
  friend const CoefExpr tanh     ( const CoefExpr& coef );
  friend const CoefExpr atanh    ( const CoefExpr& coef );

  // Binary arithmetic operations with pdfs.
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
