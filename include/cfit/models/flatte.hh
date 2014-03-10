#ifndef __FLATTE_HH__
#define __FLATTE_HH__

#include <cfit/parameter.hh>
#include <cfit/resonance.hh>

class PhaseSpace;

// Flatté resonance.
class Flatte : public Resonance
{
public:
  template <class T>
  Flatte( const T&         resoA ,
          const T&         resoB ,
          const Parameter& mass  ,
          const Parameter& width ,
          const Parameter& r     ,
          const Parameter& gamma1,
          const Parameter& gamma2,
          const Parameter& m02a  ,
          const Parameter& m02b  ,
          const int&       l      )
    : Resonance( resoA, resoB, mass, width, r, l )
    {
      push( gamma1 ); // Gamma coefficients of each channel.
      push( gamma2 );
      push( m02a   ); // Mass of the second channel resonance.
      push( m02b   );
    }

  double gamma1() const { return getPar( 0 ); }
  double gamma2() const { return getPar( 1 ); }
  double m02a()   const { return getPar( 2 ); }
  double m02b()   const { return getPar( 3 ); }

  double gamma1Sq() const { return std::pow( gamma1(), 2 ); }
  double gamma2Sq() const { return std::pow( gamma2(), 2 ); }
  double m02aSq()   const { return std::pow( m02a()  , 2 ); }
  double m02bSq()   const { return std::pow( m02b()  , 2 ); }

  std::complex< double > propagator( const PhaseSpace& ps, const double& mSqAB ) const;
  Flatte*                copy()                                                  const;
};

#endif
