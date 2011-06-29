#ifndef __RELBREITWIGNER_HH__
#define __RELBREITWIGNER_HH__

#include <cctype> // For the tolower function.

#include <cfit/parameter.hh>
#include <cfit/resonance.hh>

class PhaseSpace;

class RelBreitWigner : public Resonance
{
public:
  RelBreitWigner( const unsigned&  resoA, const unsigned&  resoB,
		  const Parameter& mass , const Parameter& width,
		  const Parameter& R    , const int& l            )
  {
    _resoA  = resoA;
    _resoB  = resoB;
    _noRes  = 6 - resoA - resoB;
    _mass   = mass;
    _width  = width;
    _R      = R;
    _l      = l;
  }

  RelBreitWigner( const int&       resoA, const int&       resoB,
		  const Parameter& mass , const Parameter& width,
		  const Parameter& R    , const int& l            )
  {
    _resoA  = resoA;
    _resoB  = resoB;
    _noRes  = 6 - resoA - resoB;
    _mass   = mass;
    _width  = width;
    _R      = R;
    _l      = l;
  }

  RelBreitWigner( const char&      resoA, const char&      resoB,
		  const Parameter& mass , const Parameter& width,
		  const Parameter& R    , const int& l            )
  {
    _resoA  = std::tolower( resoA ) - 'a' + 1;
    _resoB  = std::tolower( resoB ) - 'a' + 1;
    _noRes  = 6 - resoA - resoB;
    _mass   = mass;
    _width  = width;
    _R      = R;
    _l      = l;
  }

  std::complex< double > propagator( const PhaseSpace& ps, const double& mSqAB ) const;
};

#endif
