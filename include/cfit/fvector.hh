#ifndef __FVECTOR_HH__
#define __FVECTOR_HH__

#include <map>
#include <vector>
#include <complex>
#include <functional>

#include <cfit/parameter.hh>
#include <cfit/coef.hh>
#include <cfit/pdfexception.hh>
#include <cfit/matrix.hh>

class Coef;
class Amplitude;
class PhaseSpace;

class Fvector
{
protected:
  unsigned  _resoA;      // Indices of the resonant particles.
  unsigned  _resoB;      //
  unsigned  _noRes;      // Index of the non-resonant particle.

  std::vector< double > _m0;  // K-matrix poles.
  Matrix< double >      _g0;  // Base residue functions.
  Matrix< double >      _fSc; // K-matrix background terms.

  double                _s0sc; // Pole of the resonant term (typically outside the phase space).
  double                _s0A;  // Adler pole (typically outside the phase space).
  double                _sA;   // Adler parameter.

  std::vector< std::pair< std::string, std::string > > _beta; // Names of real and imag parts.
  std::vector< std::pair< std::string, std::string > > _fPr;  // Name of fPr elements.
  std::string                                          _s0pr; // Name of s0pr.

  static const double& _mPi;

  std::map< const std::string, Parameter > _parMap;
  std::vector< std::string >               _parOrder;

  std::complex< double > rho   ( const double& mCh, const double& mSqAB ) const;
  std::complex< double > rho4pi(                    const double& mSqAB ) const;
  std::complex< double > rho   ( const int& index , const double& mSqAB ) const;

  double getPar( const std::string& name ) const
    {
      return _parMap.find( name )->second.value();
    }

  std::complex< double > getCoef( const std::pair< std::string, std::string >& name ) const
    {
      return std::complex< double >( _parMap.find( name.first  )->second.value(),
                                     _parMap.find( name.second )->second.value() );
    }

public:
  template <class T>
  Fvector( const T&                     resoA, const T&                resoB,
           const std::vector< double >& m0   , const Matrix< double >& g0   ,
           const Matrix< double >&      fSc  , const double&           s0sc ,
           const double&                s0A  , const double&           sA   ,
           const std::vector< Coef >&   beta ,
           const std::vector< Coef >&   fPr  , const Parameter&        s0pr  )
    : _m0( m0 ), _g0( g0 ), _fSc( fSc ), _s0sc( s0sc ), _s0A( s0A ), _sA( sA )
  {
    _resoA  = resoA;
    _resoB  = resoB;
    _noRes  = 6 - resoA - resoB;

    pushBeta( beta );
    pushfPr ( fPr  );
    pushS0pr( s0pr );
  }

  ~Fvector() {};

  void pushBeta( const std::vector< Coef >& beta );
  void pushfPr ( const std::vector< Coef >& fPr  );
  void pushS0pr( const Parameter&           s0pr );

  void setPars( const std::map< std::string, Parameter >& pars );

  // AB is the resonant pair, with A the first and B the second particle in the pair.
  //    Order is only relevant for the sign of the Zemach angular term for l = 1.
  // m2ij functions below select the squared invariant mass according to the given
  //    _reso1 and _reso2 values.
  double                 m2AB               ( const double& mSq12, const double& mSq13, const double& mSq23 ) const;
  double                 m2AC               ( const double& mSq12, const double& mSq13, const double& mSq23 ) const;
  double                 m2BC               ( const double& mSq12, const double& mSq13, const double& mSq23 ) const;

  static double          kallen             ( const double& x    , const double& y    , const double& z     );

  // Momentum of a resonant particle in the rest frame of the resonant pair.
  double                 q                  ( const PhaseSpace& ps, const double& mSqAB )                     const;
  double                 qSq                ( const PhaseSpace& ps, const double& mSqAB )                     const;

  // Phase space factor, equal to 2q/m.
  double                 rho                ( const PhaseSpace& ps, const double& mSqAB )                     const;
  std::complex< double > evaluate           ( const PhaseSpace& ps,
                                              const double&     mSqAB,
                                              const double&     mSqAC,
                                              const double&     mSqBC )                                       const;

  virtual std::complex< double > propagator( const PhaseSpace& ps, const double& mSqAB ) const;
  virtual Fvector*               copy()                                                  const
    {
      return new Fvector( *this );
    }

  // Operations of resonances with themselves.
  friend const Amplitude operator+( const Fvector&   left, const Fvector&   right );
  friend const Amplitude operator-( const Fvector&   left, const Fvector&   right );

  // Operations of constants and resonances.
  friend const Amplitude operator+( const double&    left, const Fvector&   right );
  friend const Amplitude operator-( const double&    left, const Fvector&   right );
  friend const Amplitude operator*( const double&    left, const Fvector&   right );

  friend const Amplitude operator+( const Fvector&   left, const double&    right );
  friend const Amplitude operator-( const Fvector&   left, const double&    right );
  friend const Amplitude operator*( const Fvector&   left, const double&    right );
  friend const Amplitude operator/( const Fvector&   left, const double&    right );

  // Operations of coefficients and resonances.
  friend const Amplitude operator+( const Coef&      left, const Fvector&   right );
  friend const Amplitude operator-( const Coef&      left, const Fvector&   right );
  friend const Amplitude operator*( const Coef&      left, const Fvector&   right );

  friend const Amplitude operator+( const Fvector&   left, const Coef&      right );
  friend const Amplitude operator-( const Fvector&   left, const Coef&      right );
  friend const Amplitude operator*( const Fvector&   left, const Coef&      right );
  friend const Amplitude operator/( const Fvector&   left, const Coef&      right );

  // Operations with resonances and amplitudes.
  friend const Amplitude operator+( const Fvector&   left, const Amplitude& right );
  friend const Amplitude operator-( const Fvector&   left, const Amplitude& right );

  friend const Amplitude operator+( const Amplitude& left, const Fvector&   right );
  friend const Amplitude operator-( const Amplitude& left, const Fvector&   right );
};



template <>
inline Fvector::Fvector( const char&                  resoA, const char&             resoB,
                         const std::vector< double >& m0   , const Matrix< double >& g0   ,
                         const Matrix< double >&      fSc  , const double&           s0sc ,
                         const double&                s0A  , const double&           sA   ,
                         const std::vector< Coef >&   beta ,
                         const std::vector< Coef >&   fPr  , const Parameter&        s0pr  )
  : _m0( m0 ), _g0( g0 ), _fSc( fSc ), _s0sc( s0sc ), _s0A( s0A ), _sA( sA )
{
  _resoA  = std::tolower( resoA ) - 'a' + 1;
  _resoB  = std::tolower( resoB ) - 'a' + 1;
  _noRes  = 6 - resoA - resoB;

  pushBeta( beta );
  pushfPr ( fPr  );
  pushS0pr( s0pr );
}


#endif
