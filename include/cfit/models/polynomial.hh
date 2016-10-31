#ifndef __POLYNOMIAL_HH__
#define __POLYNOMIAL_HH__

#include <vector>
#include <cmath>

#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/pdfmodel.hh>


class Polynomial : public PdfModel
{
private:
  std::vector< ParameterExpr > _coefs;

  bool   _hasLower;
  bool   _hasUpper;
  double _lower;
  double _upper;

  double _norm;

  void setParExpr();

public:
  Polynomial( const Variable& x, const Parameter& c1 );
  Polynomial( const Variable& x, const Parameter& c1, const Parameter& c2 );
  Polynomial( const Variable& x, const std::vector< Parameter >& coefs );

  Polynomial( const Variable& x, const ParameterExpr& c1 );
  Polynomial( const Variable& x, const ParameterExpr& c1, const ParameterExpr& c2 );
  Polynomial( const Variable& x, const std::vector< ParameterExpr >& coefs );

  Polynomial* copy() const;

  // Getters.
  double coef( const unsigned& index ) const;

  void setLowerLimit  ( const double& lower );
  void setUpperLimit  ( const double& upper );
  void setLimits      ( const double& lower, const double& upper );

  void cache();

  const double evaluate( const double& x                   ) const throw( PdfException );
  const double evaluate( const std::vector< double >& vars ) const throw( PdfException );

  const double area    ( const double& min, const double& max ) const throw( PdfException );
};

#endif
