#ifndef __ARGUS_HH__
#define __ARGUS_HH__

#include <vector>
#include <cmath>

#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/pdfmodel.hh>


class Argus : public PdfModel
{
private:
  ParameterExpr _c;
  ParameterExpr _chi;

  bool   _hasLower;
  bool   _hasUpper;
  double _lower;
  double _upper;

  double _norm;

  void setParExpr();

public:
  Argus( const Variable& x, const Parameter&     c, const Parameter&     chi );
  Argus( const Variable& x, const ParameterExpr& c, const ParameterExpr& chi );

  Argus* copy() const;

  // Getters.
  double c()   const;
  double chi() const;

  void setLowerLimit  ( const double& lower );
  void setUpperLimit  ( const double& upper );
  void setLimits      ( const double& lower, const double& upper );
  void unsetLowerLimit();
  void unsetUpperLimit();
  void unsetLimits    ();

  void cache();

  const double evaluate( const double& x                   ) const throw( PdfException );
  const double evaluate( const std::vector< double >& vars ) const throw( PdfException );

  const double area( const double& min, const double& max ) const throw( PdfException );
};

#endif
