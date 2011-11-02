#ifndef __PHASESPACE_HH__
#define __PHASESPACE_HH__

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

  static const double kallen( const double& x, const double& y, const double& z );

  const double mSq13min( const double& mSq12 ) const;
  const double mSq13max( const double& mSq12 ) const;
  const double mSq23min( const double& mSq12 ) const;
  const double mSq23max( const double& mSq12 ) const;

public:
  PhaseSpace() {};
  PhaseSpace( const double& mMother, const double& m1, const double& m2, const double& m3 )
    : _mMother  ( mMother           ), _m1  ( m1      ), _m2  ( m2      ), _m3  ( m3      ),
      _mSqMother( mMother * mMother ), _mSq1( m1 * m1 ), _mSq2( m2 * m2 ), _mSq3( m3 * m3 )
  {}

  // Getters.
  const double mMother()   const { return _mMother;   }
  const double m1()        const { return _m1;        }
  const double m2()        const { return _m2;        }
  const double m3()        const { return _m3;        }

  const double mSqMother() const { return _mSqMother; }
  const double mSq1()      const { return _mSq1;      }
  const double mSq2()      const { return _mSq2;      }
  const double mSq3()      const { return _mSq3;      }

  const double m  ( unsigned index ) const;
  const double mSq( unsigned index ) const;

  // Check if the kinematically allowed region contains a given point.
  bool contains( const double& mSq12, const double& mSq13, const double& mSq23 ) const;
};

#endif
