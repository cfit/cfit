#ifndef __RESONANCE_HH__
#define __RESONANCE_HH__

#include <complex>
#include <cfit/parameter.hh>

class Coef;
class Amplitude;
class PhaseSpace;

class Resonance
{
protected:
  unsigned  _resoA;      // Indices of the resonant particles.
  unsigned  _resoB;      //
  unsigned  _noRes;      // Index of the non-resonant particle.
  int       _l;          // Angular momentum quantum number.

  Parameter _mass;       // Squared mass of the resonant propagator.
  Parameter _width;      // Width of the resonant propagator.
  Parameter _R;          // Radius of the Blatt-Weisskopf centrifugal barrier factor.
public:
  virtual ~Resonance() {};

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

  // Phase space factor, equal to 2q/m.
  double                 rho                ( const PhaseSpace& ps, const double& mSqAB )                     const;
  double                 runningWidth       ( const PhaseSpace& ps, const double& mSqAB )                     const;
  double                 blattWeisskopfPrime( const PhaseSpace& ps, const double& mSqAB )                     const;
  double                 blattWeisskopf     ( const PhaseSpace& ps, const double& mSqAB )                     const;
  double                 zemach             ( const PhaseSpace& ps,
					      const double&     mSqAB,
					      const double&     mSqAC,
					      const double&     mSqBC )                                       const;
  std::complex< double > evaluate           ( const PhaseSpace& ps,
					      const double&     mSqAB,
					      const double&     mSqAC,
					      const double&     mSqBC )                                       const;

  virtual std::complex< double > propagator( const PhaseSpace& ps, const double& mSqAB ) const = 0;



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

#endif
