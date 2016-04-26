#ifndef __BINNEDAMPLITUDE_HH__
#define __BINNEDAMPLITUDE_HH__

#include <vector>
#include <string>
#include <complex>
#include <map>

#include <cfit/operation.hh>
#include <cfit/exceptions.hh>
#include <cfit/parameter.hh>
#include <cfit/coef.hh>
#include <cfit/coefexpr.hh>
#include <cfit/resonance.hh>
#include <cfit/fvector.hh>
#include <cfit/binning.hh>


class PhaseSpace;

class BinnedAmplitude
{
private:
  std::map< std::string, Parameter > _parMap;

  unsigned _nbins;

  // Vectors of the T_{+b}, T_{-b} parameters and X_b coefficients.
  std::vector< Parameter > _tpb;
  std::vector< Parameter > _tmb;
  std::vector< CoefExpr  > _xb;

public:
  BinnedAmplitude() = default;

  BinnedAmplitude( const std::vector< Parameter >& tpb,
                   const std::vector< Parameter >& tmb,
                   const std::vector< CoefExpr  >& xb  );

  ~BinnedAmplitude() {}

  // const bool isFixed() const;

  const std::map< std::string, Parameter >& getPars() const { return _parMap; };
  void setPars( const std::map< std::string, Parameter >& pars );

  const std::vector< Parameter >& tpb() const { return _tpb; }
  const std::vector< Parameter >& tmb() const { return _tmb; }
  const std::vector< CoefExpr  >& xb()  const { return _xb;  }

  const double                 getT( const int& bin ) const
  {
    if ( std::abs( bin ) > _nbins )
      throw PdfException( "BinnedAmplitude: requesting invalid bin" );

    return ( ( bin > 0 ) ? _tpb[ bin - 1 ] : _tmb[ std::abs( bin ) - 1 ] ).value();
  }

  const std::complex< double > getX( const int& bin ) const
  {
    if ( std::abs( bin ) > _nbins )
      throw PdfException( "BinnedAmplitude: requesting invalid bin" );

    return ( bin > 0 ) ? _xb[ bin - 1 ].evaluate() : std::conj( _xb[ std::abs( bin ) - 1 ].evaluate() );
  }


  const std::tuple< double, double, std::complex< double > > evaluate( const int& bin ) const
  {
    const double&&                 tpb = getT(   bin );
    const double&&                 tmb = getT( - bin );
    const std::complex< double >&& txb = std::sqrt( tpb * tmb ) * getX( bin );
    return std::tuple< double, double, std::complex< double > >( tpb, tmb, txb );
  }


  // // Assignment operations.
  // const BinnedAmplitude& operator= ( const double&                 ctnt );
  // const BinnedAmplitude& operator= ( const std::complex< double >& ctnt );
  // const BinnedAmplitude& operator= ( const Coef&                   coef );
  // const BinnedAmplitude& operator= ( const CoefExpr&               expr );

  // const BinnedAmplitude& operator= ( const BinnedAmplitude&        right );

  // const BinnedAmplitude& operator*=( const double&                 right );
  // const BinnedAmplitude& operator/=( const double&                 right );

  // const BinnedAmplitude& operator*=( const std::complex< double >& right );
  // const BinnedAmplitude& operator/=( const std::complex< double >& right );

  // const BinnedAmplitude& operator*=( const Coef&                   right );
  // const BinnedAmplitude& operator/=( const Coef&                   right );

  // const BinnedAmplitude& operator*=( const CoefExpr&               right );
  // const BinnedAmplitude& operator/=( const CoefExpr&               right );



  // // Operations with constants and amplitudes.
  // friend const BinnedAmplitude operator*( const double&    left, const BinnedAmplitude& right );

  // friend const BinnedAmplitude operator*( const BinnedAmplitude& left, const double&    right );
  // friend const BinnedAmplitude operator/( const BinnedAmplitude& left, const double&    right );


  // // Operations with coefficients and amplitudes.
  // friend const BinnedAmplitude operator*( const Coef&      left, const BinnedAmplitude& right );

  // friend const BinnedAmplitude operator*( const BinnedAmplitude& left, const Coef&      right );
  // friend const BinnedAmplitude operator/( const BinnedAmplitude& left, const Coef&      right );


  // // Operations with coefficient expressions and amplitudes.
  // friend const BinnedAmplitude operator*( const CoefExpr&  left, const BinnedAmplitude& right );

  // friend const BinnedAmplitude operator*( const BinnedAmplitude& left, const CoefExpr&  right );
  // friend const BinnedAmplitude operator/( const BinnedAmplitude& left, const CoefExpr&  right );
};

#endif
