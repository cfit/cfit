#ifndef __EXPOGAUSS_HH__
#define __EXPOGAUSS_HH__

#include <vector>
#include <cmath>

#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/pdfmodel.hh>


class ExpoGauss : public PdfModel
{
private:
  ParameterExpr _gamma;
  ParameterExpr _mu;
  ParameterExpr _sigma;

  bool   _hasLower;
  bool   _hasUpper;
  double _lower;
  double _upper;

  double _normExpo;
  double _normGauss;

  double _norm;

  const double expogauss( const double& x ) const;

  void setParExpr();

public:
  ExpoGauss( const Variable& x,
             const Parameter&     gamma, const Parameter&     mu, const Parameter&     sigma );

  ExpoGauss( const Variable& x,
             const ParameterExpr& gamma, const ParameterExpr& mu, const ParameterExpr& sigma );

  ExpoGauss* copy() const;

  // Getters.
  double gamma() const;
  double mu()    const;
  double sigma() const;

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
