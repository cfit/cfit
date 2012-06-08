#ifndef __EXPONENTIAL_HH__
#define __EXPONENTIAL_HH__

#include <vector>
#include <cmath>

#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/pdfmodel.hh>


class Exponential : public PdfModel
{
private:
  bool   _hasLower;
  bool   _hasUpper;
  double _lower;
  double _upper;

  double _norm;

public:
  Exponential( const Variable& x, const Parameter& gamma );

  // Getters.
  double gamma() const;

  void setLowerLimit  ( const double& lower );
  void setUpperLimit  ( const double& upper );
  void setLimits      ( const double& lower, const double& upper );
  void unsetLowerLimit();
  void unsetUpperLimit();
  void unsetLimits    ();

  void cache();

  double evaluate(                                   ) const throw( PdfException );
  double evaluate( double x                          ) const throw( PdfException );
  double evaluate( const std::vector< double >& vars ) const throw( PdfException );
};

#endif
