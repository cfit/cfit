#ifndef __PHASESPACE_HH__
#define __PHASESPACE_HH__

#include <cmath>
#include <exception>

class PhaseSpace
{
private:
  double _mMother; // Mass of the mother.
  double _m1;      //
  double _m2;      // Masses of the three daughters.
  double _m3;      //

  double _mSqMother; // Mass of the mother.
  double _mSq1;      //
  double _mSq2;      // Masses of the three daughters.
  double _mSq3;      //

  double _mSqSum;

  static const double kallen( const double& x, const double& y, const double& z );

public:
  PhaseSpace();
  PhaseSpace( const double& mMother, const double& m1, const double& m2, const double& m3 );

  // Getters.
  const double mMother()   const { return _mMother;   }
  const double m1()        const { return _m1;        }
  const double m2()        const { return _m2;        }
  const double m3()        const { return _m3;        }

  const double mSqMother() const { return _mSqMother; }
  const double mSq1()      const { return _mSq1;      }
  const double mSq2()      const { return _mSq2;      }
  const double mSq3()      const { return _mSq3;      }

  const double mSqSum()    const { return _mSqSum;    }

  const double m  ( unsigned index ) const;
  const double mSq( unsigned index ) const;

  const double mSq13min( const double& mSq12 ) const;
  const double mSq13max( const double& mSq12 ) const;
  const double mSq23min( const double& mSq12 ) const;
  const double mSq23max( const double& mSq12 ) const;

  const double mSq12min() const { return std::pow( _m1      + _m2, 2 ); };
  const double mSq12max() const { return std::pow( _mMother - _m3, 2 ); };
  const double mSq13min() const { return std::pow( _m1      + _m3, 2 ); };
  const double mSq13max() const { return std::pow( _mMother - _m2, 2 ); };
  const double mSq23min() const { return std::pow( _m2      + _m3, 2 ); };
  const double mSq23max() const { return std::pow( _mMother - _m1, 2 ); };

  const double mSqMin( const unsigned& index ) const
  {
    if ( index == 0 ) return mSq12min();
    if ( index == 1 ) return mSq13min();
    if ( index == 2 ) return mSq23min();

    throw std::exception();
  };

  const double mSqMax( const unsigned& index ) const
  {
    if ( index == 0 ) return mSq12max();
    if ( index == 1 ) return mSq13max();
    if ( index == 2 ) return mSq23max();

    throw std::exception();
  };

  // Check if the kinematically allowed region contains a given point.
  bool contains( const double& mSq12, const double& mSq13, const double& mSq23 ) const;
  bool contains( const double& mSq12, const double& mSq13                      ) const;
};

#endif
