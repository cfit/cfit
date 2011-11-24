#ifndef __GLASS_HH__
#define __GLASS_HH__

#include <cfit/parameter.hh>
#include <cfit/resonance.hh>

class PhaseSpace;

// Generalized Lass resonance.
class GLass : public Resonance
{
public:
  template <class T>
  GLass( const T&         resoA, const T&         resoB,
         const Parameter& mass , const Parameter& width,
         const Parameter& r    ,
         const Parameter& lassR, const Parameter& lassB,
         const Parameter& phiR , const Parameter& phiB ,
         const Parameter& lassr, const Parameter& lassa,
         const int&       l      )
    : Resonance( resoA, resoB, mass, width, r, l )
    {
      push( lassR );
      push( lassB );
      push( phiR  );
      push( phiB  );
      push( lassr );
      push( lassa );
    }

  double lassR() const { return getPar( 0 ); }
  double lassB() const { return getPar( 1 ); }
  double phiR()  const { return getPar( 2 ); }
  double phiB()  const { return getPar( 3 ); }
  double lassr() const { return getPar( 4 ); }
  double lassa() const { return getPar( 5 ); }

  std::complex< double > propagator( const PhaseSpace& ps, const double& mSqAB ) const;
  GLass*                 copy()                                                  const;
};

#endif
