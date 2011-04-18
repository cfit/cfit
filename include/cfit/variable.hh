#ifndef __VARIABLE_HH__
#define __VARIABLE_HH__

#include <string>

class Variable
{
private:
  std::string _name;
  double      _value;
  double      _error;

public:
  Variable() {}

  Variable( std::string name, double value = 0., double error = 0. )
    : _name( name ), _value( value ), _error( error )
  {}

  void set( double value, double error = -1. )
  {
    _value = value;
    if ( error >= 0. )
      _error = error;
  }
  void setValue( double value ) { _value = value; }
  void setError( double error ) { _error = error; }

  // Getters.
  std::string name()  const { return _name;  }
  double      value() const { return _value; }
  double      error() const { return _error; }
};


#endif
