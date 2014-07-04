#ifndef __RESONANCE_HH__
#define __RESONANCE_HH__

#include <map>
#include <vector>
#include <complex>

#include <cfit/parameter.hh>
#include <cfit/exceptions.hh>

class Coef;
class Amplitude;
class PhaseSpace;

class Resonance
{
  friend class Amplitude;
protected:
  unsigned _resoA;      // Indices of the resonant particles.
  unsigned _resoB;      //
  unsigned _noRes;      // Index of the non-resonant particle.
  int      _l;          // Angular momentum quantum number.

  bool     _helicity;   // Use helicity formalism for angular distribution, instead of Zemach.
  bool     _twoBW;      // Use two Blatt-Weisskopf centrifugal terms, instead of one.

  std::map< const std::string, Parameter > _parMap;
  std::vector< std::string >               _parOrder;
public:
  template <class T>
  Resonance( const T&         resoA, const T&         resoB,
             const Parameter& mass , const Parameter& width,
             const Parameter& r    , const int&       l      )
  {
    _resoA  = resoA;
    _resoB  = resoB;
    _noRes  = 6 - resoA - resoB;
    _l      = l;

    _helicity = false;
    _twoBW    = false;

    push( mass  );
    push( width );
    push( r     );
  }

  virtual ~Resonance() {};

  void push( const Parameter& par );

  void useHelicity( const bool helicity = true ) { _helicity = helicity; }
  void useTwoBW   ( const bool twoBW    = true ) { _twoBW    = twoBW;    }

  // For resonances with larger number of parameters, be able to get them by index.
  //    Important: the zeroth extra parameter is the 3rd element in the vector.
  double getPar( const unsigned index ) const throw( PdfException );

  const double mass()   const { return _parMap.find( _parOrder[ 0 ] )->second.value(); }
  const double m()      const { return _parMap.find( _parOrder[ 0 ] )->second.value(); }
  const double width()  const { return _parMap.find( _parOrder[ 1 ] )->second.value(); }
  const double r()      const { return _parMap.find( _parOrder[ 2 ] )->second.value(); }
  const double radius() const { return _parMap.find( _parOrder[ 2 ] )->second.value(); }

  const double mSq()    const { return std::pow( m(), 2 );                             }
  const double mGamma() const { return m() * width();                                  }

  void setPars( const std::map< std::string, Parameter >& pars );

  // AB is the resonant pair, with A the first and B the second particle in the pair.
  //    Order is only relevant for the sign of the Zemach angular term for l = 1.
  // m2ij functions below select the squared invariant mass according to the given
  //    _reso1 and _reso2 values.
  double                 m2AB                ( const double& mSq12, const double& mSq13, const double& mSq23 ) const;
  double                 m2AC                ( const double& mSq12, const double& mSq13, const double& mSq23 ) const;
  double                 m2BC                ( const double& mSq12, const double& mSq13, const double& mSq23 ) const;

  static double          kallen              ( const double& x    , const double& y    , const double& z     );

  // Momentum of a resonant particle in the rest frame of the resonant pair.
  double                 q                   ( const PhaseSpace& ps, const double& mSqAB )                     const;
  double                 qSq                 ( const PhaseSpace& ps, const double& mSqAB )                     const;

  // Momentum of the non-resonant particle in the rest frame of the resonant pair.
  double                 p                   ( const PhaseSpace& ps, const double& mSqAB )                     const;
  double                 pSq                 ( const PhaseSpace& ps, const double& mSqAB )                     const;

  // Phase space factor, equal to 2q/m.
  double                 rho                 ( const PhaseSpace& ps, const double& mSqAB )                     const;
  double                 runningWidth        ( const PhaseSpace& ps, const double& mSqAB )                     const;
  double                 blattWeisskopfPrime ( const PhaseSpace& ps, const double& mSqAB )                     const;
  double                 blattWeisskopfPrimeP( const PhaseSpace& ps, const double& mSqAB )                     const;
  double                 blattWeisskopf      ( const PhaseSpace& ps, const double& mSqAB )                     const;
  double                 zemach              ( const PhaseSpace& ps,
                                               const double&     mSqAB,
                                               const double&     mSqAC,
                                               const double&     mSqBC )                                       const;
  double                 helicity            ( const PhaseSpace& ps,
                                               const double&     mSqAB,
                                               const double&     mSqAC,
                                               const double&     mSqBC )                                       const;

  std::complex< double > evaluate            ( const PhaseSpace& ps,
                                               const double&     mSqAB,
                                               const double&     mSqAC,
                                               const double&     mSqBC )                                       const;

  virtual std::complex< double > propagator( const PhaseSpace& ps, const double& mSqAB ) const = 0;
  virtual Resonance*             copy()                                                  const = 0;

  // Operations of resonances with themselves.
  friend const Amplitude operator+( const Resonance& left, const Resonance& right );
  friend const Amplitude operator-( const Resonance& left, const Resonance& right );

  // Operations of constants and resonances.
  friend const Amplitude operator+( const double&    left, const Resonance& right );
  friend const Amplitude operator-( const double&    left, const Resonance& right );
  friend const Amplitude operator*( const double&    left, const Resonance& right );

  friend const Amplitude operator+( const Resonance& left, const double&    right );
  friend const Amplitude operator-( const Resonance& left, const double&    right );
  friend const Amplitude operator*( const Resonance& left, const double&    right );
  friend const Amplitude operator/( const Resonance& left, const double&    right );

  // Operations of coefficients and resonances.
  friend const Amplitude operator+( const Coef&      left, const Resonance& right );
  friend const Amplitude operator-( const Coef&      left, const Resonance& right );
  friend const Amplitude operator*( const Coef&      left, const Resonance& right );

  friend const Amplitude operator+( const Resonance& left, const Coef&      right );
  friend const Amplitude operator-( const Resonance& left, const Coef&      right );
  friend const Amplitude operator*( const Resonance& left, const Coef&      right );
  friend const Amplitude operator/( const Resonance& left, const Coef&      right );

  // Operations with resonances and amplitudes.
  friend const Amplitude operator+( const Resonance& left, const Amplitude& right );
  friend const Amplitude operator-( const Resonance& left, const Amplitude& right );

  friend const Amplitude operator+( const Amplitude& left, const Resonance& right );
  friend const Amplitude operator-( const Amplitude& left, const Resonance& right );
};

template <>
inline Resonance::Resonance( const char&      resoA, const char&      resoB,
                             const Parameter& mass , const Parameter& width,
                             const Parameter& r    , const int&       l      )
{
  _resoA  = std::tolower( resoA ) - 'a' + 1;
  _resoB  = std::tolower( resoB ) - 'a' + 1;
  _noRes  = 6 - resoA - resoB;
  _l      = l;

  _helicity = false;
  _twoBW    = false;

  push( mass  );
  push( width );
  push( r     );
}

#endif
