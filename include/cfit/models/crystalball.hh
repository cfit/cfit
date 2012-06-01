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
  double _normCore;
  double _normTail;
  double _norm;

  const double core( const double& x ) const;
  const double tail( const double& x ) const;

public:
  CrystalBall( const Variable& x, const Parameter& mu, const Parameter& sigma, const Parameter& alpha, const Parameter& n );

  // Getters.
  double mu()    const;
  double sigma() const;
  double alpha() const;
  double n()     const;

  void cache();

  double evaluate(                                   ) const throw( PdfException );
  double evaluate( double x                          ) const throw( PdfException );
  double evaluate( const std::vector< double >& vars ) const throw( PdfException );
};

#endif
