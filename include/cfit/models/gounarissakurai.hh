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

  // Decide whether to use the BaBar buggy or corrected GS propagator.
  bool _buggy;

public:
  template <class T>
  GounarisSakurai( const T&         resoA, const T&         resoB,
                   const Parameter& mass , const Parameter& width,
                   const Parameter& r    , const int&       l    ,
                   const bool       buggy = true )
    : Resonance( resoA, resoB, mass, width, r, l ), _buggy( buggy )
    {}

  GounarisSakurai( const GounarisSakurai& right )
    : Resonance( right ), _buggy( right._buggy )
    {}

  std::complex< double > propagator( const PhaseSpace& ps, const double& mSqAB ) const;
  GounarisSakurai*       copy()                                                  const;
};

#endif
