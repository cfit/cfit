#ifndef __CRYSTALBALL_HH__
#define __CRYSTALBALL_HH__

#include <vector>
#include <cmath>

#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/pdfmodel.hh>


class CrystalBall : public PdfModel
{
private:
  bool   _hasLower;
  bool   _hasUpper;
  double _lower;
  double _upper;

  double _norm;

  const double cumulativeNorm( const double& x ) const;

  const double core( const double& x ) const;
  const double tail( const double& x ) const;

public:
  CrystalBall( const Variable& x, const Parameter& mu, const Parameter& sigma, const Parameter& alpha, const Parameter& n );

  CrystalBall* copy() const;

  // Getters.
  double mu()    const;
  double sigma() const;
  double alpha() const;
  double n()     const;

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

  const std::map< std::string, double > generate() const throw( PdfException );
};

#endif
