#ifndef __FUNCTION_HH__
#define __FUNCTION_HH__

#include <sstream>

#include <string>
#include <vector>
#include <map>

#include <cfit/exceptions.hh>
#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/parameterexpr.hh>
#include <cfit/operation.hh>


class FunctionMinimum;

class Function
{
private:
  std::map< std::string, Variable  > _varMap;
  std::map< std::string, Parameter > _parMap;

  std::string                  _expression;
  std::vector< Operation::Op > _opers;
  std::vector< double        > _ctnts;
  std::vector< std::string   > _varbs;
  std::vector< std::string   > _parms;

  void append( const double&        ctnt );
  void append( const Variable&      var  );
  void append( const Parameter&     par  );
  void append( const ParameterExpr& expr );
  void append( const Function&      func );

  void clear();

  template< class L, class R >
  Function( const L& left, const R& right, const Operation::Op& oper )
  {
    append( left  );
    append( right );

    _expression += "b"; // b = binary operation.
    _opers.push_back( oper );
  }

  template< class T >
  Function( const T& var, const Operation::Op& oper )
  {
    append( var );

    _expression += "u"; // u = unary operation.
    _opers.push_back( oper );
  }

public:
  Function() {};

  // Constructor from other objects.
  // arg could be a variable, parameter, parameter expression, or constant.
  template< class T >
  explicit Function( const T& arg )
  {
    append( arg );
  }


  // Assignment operators.
  // arg could be a variable, parameter, parameter expression, or constant.
  template< class T >
  inline const Function& operator=( const T& arg )
  {
    clear();
    append( arg );

    return *this;
  }

  // Setters.
  void setPar( const std::string& name, const double& val, const double& err = -1. ) throw( PdfException );

  void setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException );
  void setPars( const FunctionMinimum&                    pars ) throw( PdfException );

  // Getters.
  const std::map< std::string, Variable  >& getVarMap() const { return _varMap; }
  const std::map< std::string, Parameter >& getParMap() const { return _parMap; }

  const bool dependsOn( const std::string& varName ) const { return _varMap.count( varName ); };

  double evaluate( const std::map< std::string, double >& varMap ) const throw( PdfException );

  // Assignment operators.
  template< class T > const Function& operator+=( const T& arg );
  template< class T > const Function& operator-=( const T& arg );
  template< class T > const Function& operator*=( const T& arg );
  template< class T > const Function& operator/=( const T& arg );

  friend const Function pow( const Function& left, const double& right );
  friend const Function pow( const Variable& left, const double& right );

  // Binary operations that need access to this class.
  // Operations with variables with themselves.
  friend const Function operator+( const Variable&      left, const Variable&      right );
  friend const Function operator-( const Variable&      left, const Variable&      right );
  friend const Function operator*( const Variable&      left, const Variable&      right );
  friend const Function operator/( const Variable&      left, const Variable&      right );

  // Operations with variables as the left operand.
  friend const Function operator+( const Variable&      left, const double&        right );
  friend const Function operator-( const Variable&      left, const double&        right );
  friend const Function operator*( const Variable&      left, const double&        right );
  friend const Function operator/( const Variable&      left, const double&        right );

  friend const Function operator+( const Variable&      left, const Parameter&     right );
  friend const Function operator-( const Variable&      left, const Parameter&     right );
  friend const Function operator*( const Variable&      left, const Parameter&     right );
  friend const Function operator/( const Variable&      left, const Parameter&     right );

  friend const Function operator+( const Variable&      left, const ParameterExpr& right );
  friend const Function operator-( const Variable&      left, const ParameterExpr& right );
  friend const Function operator*( const Variable&      left, const ParameterExpr& right );
  friend const Function operator/( const Variable&      left, const ParameterExpr& right );

  // Operations with variables as the right operand.
  friend const Function operator+( const double&        left, const Variable&      right );
  friend const Function operator-( const double&        left, const Variable&      right );
  friend const Function operator*( const double&        left, const Variable&      right );
  friend const Function operator/( const double&        left, const Variable&      right );

  friend const Function operator+( const Parameter&     left, const Variable&      right );
  friend const Function operator-( const Parameter&     left, const Variable&      right );
  friend const Function operator*( const Parameter&     left, const Variable&      right );
  friend const Function operator/( const Parameter&     left, const Variable&      right );

  friend const Function operator+( const ParameterExpr& left, const Variable&      right );
  friend const Function operator-( const ParameterExpr& left, const Variable&      right );
  friend const Function operator*( const ParameterExpr& left, const Variable&      right );
  friend const Function operator/( const ParameterExpr& left, const Variable&      right );

  // Operations with functions with themselves.
  friend const Function operator+( const Function&      left, const Function&      right );
  friend const Function operator-( const Function&      left, const Function&      right );
  friend const Function operator*( const Function&      left, const Function&      right );
  friend const Function operator/( const Function&      left, const Function&      right );

  // Operations with functions as the left operand.
  friend const Function operator+( const Function&      left, const double&        right );
  friend const Function operator-( const Function&      left, const double&        right );
  friend const Function operator*( const Function&      left, const double&        right );
  friend const Function operator/( const Function&      left, const double&        right );

  friend const Function operator+( const Function&      left, const Variable&      right );
  friend const Function operator-( const Function&      left, const Variable&      right );
  friend const Function operator*( const Function&      left, const Variable&      right );
  friend const Function operator/( const Function&      left, const Variable&      right );

  friend const Function operator+( const Function&      left, const Parameter&     right );
  friend const Function operator-( const Function&      left, const Parameter&     right );
  friend const Function operator*( const Function&      left, const Parameter&     right );
  friend const Function operator/( const Function&      left, const Parameter&     right );

  friend const Function operator+( const Function&      left, const ParameterExpr& right );
  friend const Function operator-( const Function&      left, const ParameterExpr& right );
  friend const Function operator*( const Function&      left, const ParameterExpr& right );
  friend const Function operator/( const Function&      left, const ParameterExpr& right );

  // Operations with functions as the right operand.
  friend const Function operator+( const double&        left, const Function&      right );
  friend const Function operator-( const double&        left, const Function&      right );
  friend const Function operator*( const double&        left, const Function&      right );
  friend const Function operator/( const double&        left, const Function&      right );

  friend const Function operator+( const Variable&      left, const Function&      right );
  friend const Function operator-( const Variable&      left, const Function&      right );
  friend const Function operator*( const Variable&      left, const Function&      right );
  friend const Function operator/( const Variable&      left, const Function&      right );

  friend const Function operator+( const Parameter&     left, const Function&      right );
  friend const Function operator-( const Parameter&     left, const Function&      right );
  friend const Function operator*( const Parameter&     left, const Function&      right );
  friend const Function operator/( const Parameter&     left, const Function&      right );

  friend const Function operator+( const ParameterExpr& left, const Function&      right );
  friend const Function operator-( const ParameterExpr& left, const Function&      right );
  friend const Function operator*( const ParameterExpr& left, const Function&      right );
  friend const Function operator/( const ParameterExpr& left, const Function&      right );
};




template< class T >
inline const Function& Function::operator+=( const T& arg )
{
  if ( _expression.empty() )
    append( arg );
  else
  {
    append( arg );
    _expression += "b"; // b = binary operation.
    _opers.push_back( Operation::plus );
  }

  return *this;
}


template< class T >
inline const Function& Function::operator-=( const T& arg )
{
  if ( _expression.empty() )
  {
    append( arg );
    _expression += "u"; // If it's the first element, assume unary minus.
    _opers.push_back( Operation::minus );
  }
  else
  {
    append( arg );
    _expression += "b"; // b = binary operation.
    _opers.push_back( Operation::minus );
  }

  return *this;
}


template< class T >
inline const Function& Function::operator*=( const T& arg )
{
  if ( _expression.empty() )
    return *this;

  append( arg );
  _expression += "b"; // b = binary operation.
  _opers.push_back( Operation::mult );

  return *this;
}


template< class T >
inline const Function& Function::operator/=( const T& arg )
{
  if ( _expression.empty() )
    return *this;

  append( arg );
  _expression += "b"; // b = binary operation.
  _opers.push_back( Operation::div );

  return *this;
}

#endif
