#ifndef __RELBREITWIGNER_HH__
#define __RELBREITWIGNER_HH__

#include <cctype> // For the tolower function.

#include <cfit/parameter.hh>
#include <cfit/resonance.hh>

class PhaseSpace;

class RelBreitWigner : public Resonance
{
public:
  template <class T>
  RelBreitWigner( const T&         resoA, const T&         resoB,
                  const Parameter& mass , const Parameter& width,
                  const Parameter& r    , const int&       l      )
    : Resonance( resoA, resoB, mass, width, r, l )
    {}

  RelBreitWigner( const RelBreitWigner& right )
    : Resonance( right )
    {}

  std::complex< double > propagator( const PhaseSpace& ps, const double& mSqAB ) const;
  RelBreitWigner*        copy()                                                  const;
};

#endif
