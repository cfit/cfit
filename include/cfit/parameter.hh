#ifndef __PARAMETER_HH__
#define __PARAMETER_HH__

class Parameter
{
private:
  std::string _name;
  double      _value; // Should contain the initial value.
  double      _error; // Should contain the initial error.
  bool        _isFixed;
  double      _lower;
  double      _upper;
  bool        _hasLimits;

public:
  Parameter() {}

  Parameter( std::string name, double value = 0., double error = 0. )
    : _name( name ), _value( value ), _error( error ), _isFixed( false ),
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
  bool        hasLimits()  const { return   _hasLimits; }
  double      upper()      const { return   _upper;     }
  double      upperLimit() const { return   _upper;     }
  double      lower()      const { return   _lower;     }
  double      lowerLimit() const { return   _lower;     }

  // Arithmetic operators.
  friend const double operator+( const Parameter& l, const Parameter& r ) { return l._value + r._value; }
  friend const double operator-( const Parameter& l, const Parameter& r ) { return l._value - r._value; }
  friend const double operator*( const Parameter& l, const Parameter& r ) { return l._value * r._value; }
  friend const double operator/( const Parameter& l, const Parameter& r ) { return l._value / r._value; }

  friend const double operator+( const double& l, const Parameter& r ) { return l + r._value; }
  friend const double operator-( const double& l, const Parameter& r ) { return l - r._value; }
  friend const double operator*( const double& l, const Parameter& r ) { return l * r._value; }
  friend const double operator/( const double& l, const Parameter& r ) { return l / r._value; }

  friend const double operator+( const Parameter& l, const double& r ) { return l._value + r; }
  friend const double operator-( const Parameter& l, const double& r ) { return l._value - r; }
  friend const double operator*( const Parameter& l, const double& r ) { return l._value * r; }
  friend const double operator/( const Parameter& l, const double& r ) { return l._value / r; }
};

#endif
