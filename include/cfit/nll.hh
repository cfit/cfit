#ifndef __NLL_HH__
#define __NLL_HH__

#include <vector>

#include <cfit/minimizer.hh>
#include <cfit/pdfbase.hh>
#include <cfit/dataset.hh>
#include <cfit/pdfexception.hh>

class Nll : public Minimizer
{
public:
  Nll( PdfBase& pdf, const Dataset& data );

  double up() const { return 1.; }
  double operator()( const std::vector<double>& par ) const throw( PdfException );
};

#endif
