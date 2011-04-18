#ifndef __GAUSS_HH__
#define __GAUSS_HH__

#include <vector>
#include <cmath>

#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/pdfmodel.hh>


class Gauss : public PdfModel
{
public:
  Gauss( const Variable& x, const Parameter& mean, const Parameter& sigma );

  // Getters.
  double mean()  const;
  double sigma() const;

  double evaluate(                                   ) const throw( PdfException );
  double evaluate( double x                          ) const throw( PdfException );
  double evaluate( const std::vector< double >& vars ) const throw( PdfException );
};

#endif
