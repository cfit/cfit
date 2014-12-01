#ifndef __GENARGUSGAUSS_HH__
#define __GENARGUSGAUSS_HH__

#include <vector>
#include <cmath>

#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/pdfmodel.hh>


class GenArgusGauss : public PdfModel
{
private:
  bool   _hasLower;
  bool   _hasUpper;
  double _lower;
  double _upper;

  double _normGenArgus;
  double _normGauss;

  double _norm;

  const double genargus     ( const double& x ) const;
  const double gauss        ( const double& x ) const;
  const double genargusgauss( const double& x ) const;

  const double genarguscore ( const double& x ) const;

  const double generateArgus() const;

public:
  GenArgusGauss( const Variable& x,
                 const Parameter& c , const Parameter& chi, const Parameter& p,
                 const Parameter& mu, const Parameter& sigma );

  GenArgusGauss* copy() const;

  // Getters.
  double c()     const;
  double chi()   const;
  double p()     const;
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
