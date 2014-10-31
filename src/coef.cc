
#include <cfit/coef.hh>

const std::complex< double > Coef::value() const
{
  return std::complex< double >( _real.value(), _imag.value() );
}

const bool Coef::isFixed() const
{
  return _real.isFixed() && _imag.isFixed();
}
