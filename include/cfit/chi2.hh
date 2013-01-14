#ifndef __CHI2_HH__
#define __CHI2_HH__

#include <vector>

#include <cfit/minimizer.hh>
#include <cfit/pdfbase.hh>
#include <cfit/variable.hh>
#include <cfit/dataset.hh>
#include <cfit/exceptions.hh>

class Chi2 : public Minimizer
{
private:
  const Variable& _y;

public:
  Chi2( PdfBase& pdf, const Variable& y, const Dataset& data );

  double up() const { return 1.; }
  double operator()( const std::vector<double>& par ) const throw( PdfException );
};

#endif
