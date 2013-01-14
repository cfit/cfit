#ifndef __EXCEPTIONS_HH__
#define __EXCEPTIONS_HH__

#include <string>
#include <exception>



class MinimizerException : public std::exception
{
private:
  std::string _what;
public:
  MinimizerException( const std::string& str )
    : _what( str )
  {}
  ~MinimizerException() throw() {}
  const char* what() const throw() { return _what.c_str(); }
};


class PdfException : public std::exception
{
private:
  std::string _what;
  bool _critical;
public:
  PdfException( std::string str )
  {
    _what = str;
  }
  ~PdfException() throw() {}
  const char* what() const throw() { return _what.c_str(); }
};

#endif
