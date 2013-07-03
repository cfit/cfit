#ifndef __PARAMETER_HH__
#define __PARAMETER_HH__

#include <string>

class PdfModel;
class PdfExpr;
class ParameterExpr;
class Variable;
class Function;

class Parameter
{
  friend class PdfExpr;

private:
  std::string _name;
  double      _value; // Should contain the initial value.
  double      _error; // Should contain the initial error.
  bool        _isFixed;
  bool        _isBlind;
  double      _lower;
  double      _upper;
  bool        _hasLimits;

public:
  Parameter() {}

  Parameter( std::string name, double value = 0., double error = 0. )
    : _name( name ), _value( value ), _error( error ),
      _isFixed( false ), _isBlind( false ),
      _lower( 0. ), _upper( 0. ), _hasLimits( false )
  {}

  // Setters.
  void set ( double value, double error = -1. )
  {
    _value = value;
    if ( error >= 0. )
      _error = error;
  }
  void setValue( double value ) { _value   = value; }
  void setError( double error ) { _error   = error; }
  void fix()                    { _isFixed = true;  }
  void release()                { _isFixed = false; }
  void blind()                  { _isBlind = true;  }
  void unblind()                { _isBlind = false; }
  void setLimits( double lower, double upper )
  {
    _lower     = lower;
    _upper     = upper;
    _hasLimits = true;
  }

  // Getters.
  std::string name()       const { return   _name;      }
  double      value()      const { return   _value;     }
  double      error()      const { return   _error;     }
  bool        isFixed()    const { return   _isFixed;   }
  bool        isReleased() const { return ! _isFixed;   }
  bool        isBlind()    const { return   _isBlind;   }
  bool        hasLimits()  const { return   _hasLimits; }
  double      upper()      const { return   _upper;     }
  double      upperLimit() const { return   _upper;     }
  double      lower()      const { return   _lower;     }
  double      lowerLimit() const { return   _lower;     }

  // Arithmetic operators.
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

//   friend const ParameterExpr operator+( const ParameterExpr& left, const ParameterExpr& right );
//   friend const ParameterExpr operator-( const ParameterExpr& left, const ParameterExpr& right );
//   friend const ParameterExpr operator*( const ParameterExpr& left, const ParameterExpr& right );
//   friend const ParameterExpr operator/( const ParameterExpr& left, const ParameterExpr& right );

//   friend const ParameterExpr operator+( const double&        left, const ParameterExpr& right );
//   friend const ParameterExpr operator-( const double&        left, const ParameterExpr& right );
//   friend const ParameterExpr operator*( const double&        left, const ParameterExpr& right );
//   friend const ParameterExpr operator/( const double&        left, const ParameterExpr& right );

//   friend const ParameterExpr operator+( const ParameterExpr& left, const double&        right );
//   friend const ParameterExpr operator-( const ParameterExpr& left, const double&        right );
//   friend const ParameterExpr operator*( const ParameterExpr& left, const double&        right );
//   friend const ParameterExpr operator/( const ParameterExpr& left, const double&        right );

  friend const ParameterExpr operator+( const ParameterExpr& left, const Parameter&     right );
  friend const ParameterExpr operator-( const ParameterExpr& left, const Parameter&     right );
  friend const ParameterExpr operator*( const ParameterExpr& left, const Parameter&     right );
  friend const ParameterExpr operator/( const ParameterExpr& left, const Parameter&     right );

  friend const ParameterExpr operator+( const Parameter&     left, const ParameterExpr& right );
  friend const ParameterExpr operator-( const Parameter&     left, const ParameterExpr& right );
  friend const ParameterExpr operator*( const Parameter&     left, const ParameterExpr& right );
  friend const ParameterExpr operator/( const Parameter&     left, const ParameterExpr& right );


  // Operations with variables as the left operand.
  friend const Function operator+( const Variable&      left, const Parameter&     right );
  friend const Function operator-( const Variable&      left, const Parameter&     right );
  friend const Function operator*( const Variable&      left, const Parameter&     right );
  friend const Function operator/( const Variable&      left, const Parameter&     right );

  // Operations with variables as the right operand.
  friend const Function operator+( const Parameter&     left, const Variable&      right );
  friend const Function operator-( const Parameter&     left, const Variable&      right );
  friend const Function operator*( const Parameter&     left, const Variable&      right );
  friend const Function operator/( const Parameter&     left, const Variable&      right );


  // Operations with functions as the left operand.
  friend const Function operator+( const Function&      left, const Parameter&     right );
  friend const Function operator-( const Function&      left, const Parameter&     right );
  friend const Function operator*( const Function&      left, const Parameter&     right );
  friend const Function operator/( const Function&      left, const Parameter&     right );

  // Operations with functions as the right operand.
  friend const Function operator+( const Parameter&     left, const Function&      right );
  friend const Function operator-( const Parameter&     left, const Function&      right );
  friend const Function operator*( const Parameter&     left, const Function&      right );
  friend const Function operator/( const Parameter&     left, const Function&      right );


  friend const ParameterExpr pow      ( const Parameter&     left, const double&        right );
  friend const ParameterExpr pow      ( const ParameterExpr& left, const double&        right );

  friend const ParameterExpr operator-( const Parameter&     par );
  friend const ParameterExpr operator-( const ParameterExpr& par );

  friend const ParameterExpr exp      ( const Parameter&     par );
  friend const ParameterExpr log      ( const Parameter&     par );
  friend const ParameterExpr sin      ( const Parameter&     par );
  friend const ParameterExpr cos      ( const Parameter&     par );
  friend const ParameterExpr tan      ( const Parameter&     par );

  friend const PdfExpr       operator*( const Parameter&     left, const PdfModel&      right );
  friend const PdfExpr       operator*( const PdfModel&      left, const Parameter&     right );
  friend const PdfExpr       operator/( const PdfModel&      left, const Parameter&     right );

  friend const PdfExpr       operator*( const Parameter&     left, const PdfExpr&       right );
  friend const PdfExpr       operator*( const PdfExpr&       left, const Parameter&     right );
  friend const PdfExpr       operator/( const PdfExpr&       left, const Parameter&     right );
};

#endif
