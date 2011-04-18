#ifndef __MINIMIZER_HH__
#define __MINIMIZER_HH__

#include <vector>

#include <Minuit/FCNBase.h>
#include <Minuit/FunctionMinimum.h>

#include <cfit/variable.hh>
#include <cfit/dataset.hh>
#include <cfit/pdfbase.hh>
#include <cfit/minimizer.hh>

class Minimizer : public FCNBase
{
protected:
  PdfBase&        _pdf;
  const Dataset&  _data;

public:
  Minimizer( PdfBase& pdf, const Dataset& data )
    : _pdf( pdf ), _data( data )
  {}

  double up() const = 0;
  double operator()( const std::vector<double>& par ) const throw( PdfException ) = 0;

  FunctionMinimum minimize() const;
};

#endif
