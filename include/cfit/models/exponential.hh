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

  Exponential* copy() const;

  // Getters.
  double gamma() const;

  void setLowerLimit  ( const double& lower );
  void setUpperLimit  ( const double& upper );
  void setLimits      ( const double& lower, const double& upper );
  void unsetLowerLimit();
  void unsetUpperLimit();
  void unsetLimits    ();

  void cache();

  const double evaluate(                                   ) const throw( PdfException );
  const double evaluate( const double& x                   ) const throw( PdfException );
  const double evaluate( const std::vector< double >& vars ) const throw( PdfException );

  const double area    ( const double& min, const double& max ) const throw( PdfException );

  const std::map< std::string, double > generate()              const throw( PdfException );
};

#endif
