#ifndef __GOUNARISSAKURAI_HH__
#define __GOUNARISSAKURAI_HH__

#include <cfit/parameter.hh>
#include <cfit/resonance.hh>

class PhaseSpace;

class GounarisSakurai : public Resonance
{
private:
  double gsf     ( const PhaseSpace& ps, const double& mSq12 ) const;
  double gsh     ( const PhaseSpace& ps, const double& mSq12 ) const;
  double gshprime( const PhaseSpace& ps, const double& mSq12 ) const;

public:
  template <class T>
  GounarisSakurai( const T&         resoA, const T&         resoB,
                   const Parameter& mass , const Parameter& width,
                   const Parameter& r    , const int&       l      )
    : Resonance( resoA, resoB, mass, width, r, l )
    {}

  GounarisSakurai( const GounarisSakurai& right )
    : Resonance( right )
    {}

  std::complex< double > propagator( const PhaseSpace& ps, const double& mSqAB ) const;
  GounarisSakurai*       copy()                                                  const;
};

#endif
