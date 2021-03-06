#ifndef __GENARGUS_HH__
#define __GENARGUS_HH__

#include <vector>
#include <cmath>

#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/pdfmodel.hh>


class GenArgus : public PdfModel
{
private:
  ParameterExpr _c;
  ParameterExpr _chi;
  ParameterExpr _p;

  bool   _hasLower;
  bool   _hasUpper;
  double _lower;
  double _upper;

  double _norm;

  void setParExpr();

public:
  GenArgus( const Variable& x, const Parameter&     c, const Parameter&     chi, const Parameter&     p );
  GenArgus( const Variable& x, const ParameterExpr& c, const ParameterExpr& chi, const ParameterExpr& p );

  GenArgus* copy() const;

  // Getters.
  double c()   const;
  double chi() const;
  double p()   const;

  void setLowerLimit  ( const double& lower );
  void setUpperLimit  ( const double& upper );
  void setLimits      ( const double& lower, const double& upper );
  void unsetLowerLimit();
  void unsetUpperLimit();
  void unsetLimits    ();

  void cache();

  const double evaluate( const double& x                   ) const throw( PdfException );
  const double evaluate( const std::vector< double >& vars ) const throw( PdfException );

  const double area    ( const double& min, const double& max ) const throw( PdfException );

  const std::map< std::string, double > generate()              const throw( PdfException );
};

#endif
