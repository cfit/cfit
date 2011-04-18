#ifndef __DALITZRESO_HH__
#define __DALITZRESO_HH__

#include <complex>

class DalitzReso
{
private:
  char   _pair1; // Resonant particles. Order is important because the angular
  char   _pair2; //    momentum term is sign sensitive for L = 1.
  char   _other; // Non resonant particle.

  int    _L;     // Angular momentum.

  double _mSqM;  // Squared mass of the mother particle.
  double _mSq1;  // Squared masses of the particles that decay.
  double _mSq2;
  double _mSq3;

  double _m1;    // sqrt( _mSq1 )

  double _m0;    // Mass of the resonant propagator.
  double _g0;    // Width of the resonant propagator.
  double _R;     // Blatt-Weisskopf radius.

  double _mSq0;  // Squared mass of the resonant propagator.

  double _p0;    // Equivalent momentum of the resonance.
  double _pSq0;
  double _rho0;  // R^2 p0^2

  double _gsd;   // Gounaris-Sakurai constant d factor. Compute it only once.

  double _lassa;    // LASS parameters. Only set by specific lass setters,
  double _lassr;    //     since they are a lot.
  double _lassB;
  double _lassPhiB;
  double _lassR;
  double _lassPhiR;

  static double kallen( double x, double y, double z );

  double p2             ( const double& mSq12 ) const { return kallen( mSq12, _mSq1, _mSq2 ) / ( 4. * mSq12 ); }
  double p              ( const double& mSq12 ) const { return std::sqrt( p2( mSq12 ) );                       }
  double rho            ( const double& mSq12 ) const { return std::pow( _R, 2 ) * p2( mSq12 );                }
  double blattWeisskopf ( const double& mSq12 ) const;
  double runningWidth   ( const double& mSq12 ) const;
  double angular        ( const double& mSq12, const double& mSq13, const double& mSq23 ) const;
  double mSquared12     ( const double& m2AB , const double& m2AC , const double& m2BC  ) const;
  double mSquared13     ( const double& m2AB , const double& m2AC , const double& m2BC  ) const;
  double mSquared23     ( const double& m2AB , const double& m2AC , const double& m2BC  ) const;

  // Utility functions for the computation of the Gounaris-Sakurai propagator.
  double gsf      ( const double& mSq12 ) const;
  double gsh      ( const double& mSq12 ) const;
  double gshprime ( const double& mSq12 ) const;

public:
  DalitzReso() {}
  DalitzReso( char pair1, char pair2, int L,
	      double mM = 1.8645, double mA = .49767, double mB = .139570, double mC = .139570 );

  // Setters.
  void setReso  ( double m0, double g0, double R = 1.5 );

  // LASS setters.
  void setLassa   ( double a    ) { _lassa    = a;    }
  void setLassr   ( double r    ) { _lassr    = r;    }
  void setLassB   ( double B    ) { _lassB    = B;    }
  void setLassPhiB( double phiB ) { _lassPhiB = phiB; }
  void setLassR   ( double R    ) { _lassR    = R;    }
  void setLassPhiR( double phiR ) { _lassPhiR = phiR; }

  std::complex< double > rbw ( const double& m2AB, const double& m2AC, const double& m2BC ) const; // Relativistic Breit-Wigner.
  std::complex< double > gs  ( const double& m2AB, const double& m2AC, const double& m2BC ) const; // Gounaris-Sakurai.
  std::complex< double > km  ( const double& m2AB, const double& m2AC, const double& m2BC ) const; // K-matrix.
  std::complex< double > lass( const double& m2AB, const double& m2AC, const double& m2BC ) const; // LASS.
};

#endif
