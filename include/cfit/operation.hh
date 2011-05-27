#ifndef __OPERATION_HH__
#define __OPERATION_HH__

#include <string>
#include <sstream>

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
};

#endif
