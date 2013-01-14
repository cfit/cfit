#ifndef __OPERATION_HH__
#define __OPERATION_HH__

#include <string>
#include <sstream>

#include <cfit/exceptions.hh>

class Operation
{
public:
  enum Op { plus, minus, mult, div, pow, exp, log, sin, cos, tan };
  static std::string tostring( Op val )
  {
    if ( val == plus  ) return "plus";
    if ( val == minus ) return "minus";
    if ( val == mult  ) return "mult";
    if ( val == div   ) return "div";
    if ( val == pow   ) return "pow";
    if ( val == exp   ) return "exp";
    if ( val == log   ) return "log";
    if ( val == sin   ) return "sin";
    if ( val == cos   ) return "cos";
    if ( val == tan   ) return "tan";

    std::ostringstream str;
    str << val;
    return str.str();
  }

  // Binary operations.
  template < class T >
  static T operate( const T& x, const T& y, const Operation::Op& oper ) throw( PdfException );

  // Unary operations.
  template < class T >
  static T operate( const T& x,             const Operation::Op& oper ) throw( PdfException );

};


// Binary operations.
template < class T >
inline T Operation::operate( const T& x, const T& y, const Operation::Op& oper ) throw( PdfException )
{
  if ( oper == Operation::plus )
    return x + y;
  if ( oper == Operation::minus )
    return x - y;
  if ( oper == Operation::mult )
    return x * y;
  if ( oper == Operation::div )
    return x / y;
  if ( oper == Operation::pow )
    return std::pow( x, y );

  throw PdfException( std::string( "Parse error: unknown binary operation " ) + Operation::tostring( oper ) + "." );
}


// Unary operations.
template < class T >
inline T Operation::operate( const T& x, const Operation::Op& oper ) throw( PdfException )
{
  if ( oper == Operation::minus )
    return -x;
  if ( oper == Operation::exp )
    return std::exp( x );
  if ( oper == Operation::log )
    return std::log( x );
  if ( oper == Operation::sin )
    return std::sin( x );
  if ( oper == Operation::cos )
    return std::cos( x );
  if ( oper == Operation::tan )
    return std::tan( x );

  throw PdfException( std::string( "Parse error: unknown unary operation " ) + Operation::tostring( oper ) + "." );
}


#endif
