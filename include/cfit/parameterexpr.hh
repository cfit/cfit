#ifndef __PARAMETEREXPR_HH__
#define __PARAMETEREXPR_HH__

#include <vector>
#include <string>

#include <cfit/parameter.hh>
#include <cfit/operation.hh>


class PdfModel;
class Pdf;

class ParameterExpr
{
  friend class Pdf;

private:
  std::vector< double >        _ctnts;
  std::vector< Parameter >     _pars;
  std::string                  _expression;
  std::vector< Operation::Op > _opers;

  void append( const double&        val  );
  void append( const Parameter&     par  );
  void append( const ParameterExpr& expr );

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

  // Binary arithmetic operations with pdfs.
  friend const Pdf           operator*( const ParameterExpr& left, const PdfModel&      right );
  friend const Pdf           operator*( const PdfModel&      left, const ParameterExpr& right );
  friend const Pdf           operator/( const PdfModel&      left, const ParameterExpr& right );

  friend const Pdf           operator*( const ParameterExpr& left, const Pdf&           right );
  friend const Pdf           operator*( const Pdf&           left, const ParameterExpr& right );
  friend const Pdf           operator/( const Pdf&           left, const ParameterExpr& right );

};

#endif
