#ifndef __PDFEXCEPTION_HH__
#define __PDFEXCEPTION_HH__

#include <string>
#include <exception>


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
