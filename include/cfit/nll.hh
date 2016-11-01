#ifndef __NLL_HH__
#define __NLL_HH__

#include <vector>

#include <cfit/minimizer.hh>
#include <cfit/pdfbase.hh>
#include <cfit/dataset.hh>
#include <cfit/exceptions.hh>

class Nll : public Minimizer
{
public:
  Nll( const PdfModel& pdf, const Dataset& data );
  Nll( const PdfExpr&  pdf, const Dataset& data );

  Nll( const Nll& nll );

  Nll* copy() const { return new Nll( *this ); }

  double operator()( const std::vector<double>& par ) const throw( PdfException );
};

#endif
