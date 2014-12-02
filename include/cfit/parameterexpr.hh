#ifndef __PARAMETEREXPR_HH__
#define __PARAMETEREXPR_HH__

#include <complex>
#include <vector>
#include <map>
#include <string>

#include <cfit/parameter.hh>
#include <cfit/operation.hh>


class PdfModel;
class PdfExpr;
class CoefExpr;
class Function;

class ParameterExpr
{
  friend class CoefExpr;
  friend class PdfExpr;
  friend class Function;

private:
  std::vector< double >        _ctnts;
  std::vector< Parameter >     _parms;
  std::string                  _expression;
  std::vector< Operation::Op > _opers;

  void append( const double&        ctnt );
  void append( const Parameter&     parm );
  void append( const ParameterExpr& expr );
  void append( const Operation::Op& oper );

  template< class L, class R >
  ParameterExpr( const L& left, const R& right, const Operation::Op& oper )
  {
    append( left  );
    append( right );

    _expression += "b"; // b = binary operation.
    _opers.push_back( oper );
  }

  template< class T >
  ParameterExpr( const T& var, const Operation::Op& oper )
  {
    append( var );

    _expression += "u"; // u = unary operation.
    _opers.push_back( oper );
  }

public:
  ParameterExpr() {};

  // Return a map from parameter names to parameters for all parameters used in the expression.
  const std::map< std::string, Parameter > getPars() const;

  const ParameterExpr& operator+=( const ParameterExpr& right );
  const ParameterExpr& operator-=( const ParameterExpr& right );
  const ParameterExpr& operator*=( const ParameterExpr& right );
  const ParameterExpr& operator/=( const ParameterExpr& right );

  const ParameterExpr& operator+=( const Parameter&     right );
  const ParameterExpr& operator-=( const Parameter&     right );
  const ParameterExpr& operator*=( const Parameter&     right );
  const ParameterExpr& operator/=( const Parameter&     right );

  const ParameterExpr& operator+=( const double&        right );
  const ParameterExpr& operator-=( const double&        right );
  const ParameterExpr& operator*=( const double&        right );
  const ParameterExpr& operator/=( const double&        right );

  // Binary arithmetic operators.
  friend const ParameterExpr operator+( const Parameter&     left, const Parameter&     right );
  friend const ParameterExpr operator-( const Parameter&     left, const Parameter&     right );
  friend const ParameterExpr operator*( const Parameter&     left, const Parameter&     right );
  friend const ParameterExpr operator/( const Parameter&     left, const Parameter&     right );

  friend const ParameterExpr operator+( const double&        left, const Parameter&     right );
  friend const ParameterExpr operator-( const double&        left, const Parameter&     right );
  friend const ParameterExpr operator*( const double&        left, const Parameter&     right );
  friend const ParameterExpr operator/( const double&        left, const Parameter&     right );

  friend const ParameterExpr operator+( const Parameter&     left, const double&        right );
  friend const ParameterExpr operator-( const Parameter&     left, const double&        right );
  friend const ParameterExpr operator*( const Parameter&     left, const double&        right );
  friend const ParameterExpr operator/( const Parameter&     left, const double&        right );

  friend const ParameterExpr operator+( const ParameterExpr& left, const ParameterExpr& right );
  friend const ParameterExpr operator-( const ParameterExpr& left, const ParameterExpr& right );
  friend const ParameterExpr operator*( const ParameterExpr& left, const ParameterExpr& right );
  friend const ParameterExpr operator/( const ParameterExpr& left, const ParameterExpr& right );

  friend const ParameterExpr operator+( const double&        left, const ParameterExpr& right );
  friend const ParameterExpr operator-( const double&        left, const ParameterExpr& right );
  friend const ParameterExpr operator*( const double&        left, const ParameterExpr& right );
  friend const ParameterExpr operator/( const double&        left, const ParameterExpr& right );

  friend const ParameterExpr operator+( const ParameterExpr& left, const double&        right );
  friend const ParameterExpr operator-( const ParameterExpr& left, const double&        right );
  friend const ParameterExpr operator*( const ParameterExpr& left, const double&        right );
  friend const ParameterExpr operator/( const ParameterExpr& left, const double&        right );

  friend const ParameterExpr operator+( const Parameter&     left, const ParameterExpr& right );
  friend const ParameterExpr operator-( const Parameter&     left, const ParameterExpr& right );
  friend const ParameterExpr operator*( const Parameter&     left, const ParameterExpr& right );
  friend const ParameterExpr operator/( const Parameter&     left, const ParameterExpr& right );

  friend const ParameterExpr operator+( const ParameterExpr& left, const Parameter&     right );
  friend const ParameterExpr operator-( const ParameterExpr& left, const Parameter&     right );
  friend const ParameterExpr operator*( const ParameterExpr& left, const Parameter&     right );
  friend const ParameterExpr operator/( const ParameterExpr& left, const Parameter&     right );

  friend const ParameterExpr pow      ( const Parameter&     left, const double&        right );
  friend const ParameterExpr pow      ( const ParameterExpr& left, const double&        right );

  // Unary minus.
  friend const ParameterExpr operator-( const Parameter&     par );
  friend const ParameterExpr operator-( const ParameterExpr& par );

  // Unary parameter operations.
  friend const ParameterExpr exp      ( const Parameter&     par );
  friend const ParameterExpr log      ( const Parameter&     par );
  friend const ParameterExpr sin      ( const Parameter&     par );
  friend const ParameterExpr cos      ( const Parameter&     par );
  friend const ParameterExpr tan      ( const Parameter&     par );

  friend const ParameterExpr exp      ( const ParameterExpr& par );
  friend const ParameterExpr log      ( const ParameterExpr& par );
  friend const ParameterExpr sin      ( const ParameterExpr& par );
  friend const ParameterExpr cos      ( const ParameterExpr& par );
  friend const ParameterExpr tan      ( const ParameterExpr& par );

  friend const CoefExpr operator+( const std::complex< double >& left, const ParameterExpr& right );
  friend const CoefExpr operator-( const std::complex< double >& left, const ParameterExpr& right );
  friend const CoefExpr operator*( const std::complex< double >& left, const ParameterExpr& right );
  friend const CoefExpr operator/( const std::complex< double >& left, const ParameterExpr& right );

  friend const CoefExpr operator+( const ParameterExpr& left, const std::complex< double >& right );
  friend const CoefExpr operator-( const ParameterExpr& left, const std::complex< double >& right );
  friend const CoefExpr operator*( const ParameterExpr& left, const std::complex< double >& right );
  friend const CoefExpr operator/( const ParameterExpr& left, const std::complex< double >& right );

  // Binary arithmetic operations with pdfs.
  friend const PdfExpr operator*( const ParameterExpr& left, const PdfModel&      right );
  friend const PdfExpr operator*( const PdfModel&      left, const ParameterExpr& right );
  friend const PdfExpr operator/( const PdfModel&      left, const ParameterExpr& right );

  friend const PdfExpr operator*( const ParameterExpr& left, const PdfExpr&       right );
  friend const PdfExpr operator*( const PdfExpr&       left, const ParameterExpr& right );
  friend const PdfExpr operator/( const PdfExpr&       left, const ParameterExpr& right );


  // Operations with variables as the left operand.
  friend const Function operator+( const Variable&      left, const ParameterExpr& right );
  friend const Function operator-( const Variable&      left, const ParameterExpr& right );
  friend const Function operator*( const Variable&      left, const ParameterExpr& right );
  friend const Function operator/( const Variable&      left, const ParameterExpr& right );

  // Operations with variables as the right operand.
  friend const Function operator+( const ParameterExpr& left, const Variable&      right );
  friend const Function operator-( const ParameterExpr& left, const Variable&      right );
  friend const Function operator*( const ParameterExpr& left, const Variable&      right );
  friend const Function operator/( const ParameterExpr& left, const Variable&      right );

  // Operations with functions as the left operand.
  friend const Function operator+( const Function&      left, const ParameterExpr& right );
  friend const Function operator-( const Function&      left, const ParameterExpr& right );
  friend const Function operator*( const Function&      left, const ParameterExpr& right );
  friend const Function operator/( const Function&      left, const ParameterExpr& right );

  // Operations with functions as the right operand.
  friend const Function operator+( const ParameterExpr& left, const Function&      right );
  friend const Function operator-( const ParameterExpr& left, const Function&      right );
  friend const Function operator*( const ParameterExpr& left, const Function&      right );
  friend const Function operator/( const ParameterExpr& left, const Function&      right );
};

#endif
