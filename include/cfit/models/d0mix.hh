#ifndef __DALITZD0MIX_HH__
#define __DALITZD0MIX_HH__

#include <vector>
#include <complex>

#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/pdfmodel.hh>

#include <cfit/models/d0mix/dalitzReso.hh>

class DalitzD0mix : public PdfModel
{
private:
  // Masses of the particles involved.
  static const double& _mPi;
  static const double& _mKs;
  static const double& _mD0;

  // Overload the setPars function to set the resonance parameters.
  void setPars( const std::vector< double >& pars ) throw( PdfException );

  // Parameter of CPV in the mixing (q/p).
  std::complex< double > _qp;

  DalitzReso* _reso; // Resonances.

  // Booleans to control if direct and crossed matrices should be computed at each iteration.
  bool  _propagatorsFixed; // True if all the propagators are fixed.
  bool  _integralsFixed;   // True if all amplitudes and propagators are fixed.
  bool* _unfixed;          // Propagators that may vary in minuit iterations.
  bool* _conjOfPrevious;   // True if a resonance is the conjugate of the previous one (e.g. K*+ after K*-).

  // Direct and crossed matrices to compute the norm of the pdf with splitting.
  std::complex< double >** _Id;
  std::complex< double >** _Ix;

  // Integrals over the Dalitz plot.
  std::complex< double > _i1;
  std::complex< double > _iChi;

  // Value of the norm (cached after setPars).
  double _norm;

  // Variables of the model.
  double m2AB() const { return getVar( 0 ).value(); }
  double m2AC() const { return getVar( 1 ).value(); }
  double m2BC() const { return getVar( 2 ).value(); }
  double t()    const { return getVar( 3 ).value(); }
//double m2BC() const { return pow( _mD0, 2 ) + pow( _mKs, 2 ) + 2. * pow( _mPi, 2 ) - m2AB() - m2AC(); }
//double t()    const { return getVar( 2 ).value(); }

  // Parameters of the model.
  double invGamma      () const { return getPar(  0 ).value(); }
  double x             () const { return getPar(  1 ).value(); }
  double y             () const { return getPar(  2 ).value(); }

  double mKstarc       () const { return getPar(  3 ).value(); }
  double wKstarc       () const { return getPar(  4 ).value(); }
  double mRho          () const { return getPar(  5 ).value(); }
  double wRho          () const { return getPar(  6 ).value(); }
  double mOmega        () const { return getPar(  7 ).value(); }
  double wOmega        () const { return getPar(  8 ).value(); }
  double mF0_980       () const { return getPar(  9 ).value(); }
  double wF0_980       () const { return getPar( 10 ).value(); }
  double mF0_1370      () const { return getPar( 11 ).value(); }
  double wF0_1370      () const { return getPar( 12 ).value(); }
  double mF2_1270      () const { return getPar( 13 ).value(); }
  double wF2_1270      () const { return getPar( 14 ).value(); }
  double mK0starc_1430 () const { return getPar( 15 ).value(); }
  double wK0starc_1430 () const { return getPar( 16 ).value(); }
  double mK2starc_1430 () const { return getPar( 17 ).value(); }
  double wK2starc_1430 () const { return getPar( 18 ).value(); }
  double mSigma        () const { return getPar( 19 ).value(); }
  double wSigma        () const { return getPar( 20 ).value(); }
  double mSigma2       () const { return getPar( 21 ).value(); }
  double wSigma2       () const { return getPar( 22 ).value(); }
  double mKstarm_1680  () const { return getPar( 23 ).value(); }
  double wKstarm_1680  () const { return getPar( 24 ).value(); }

  double rKstarm       () const { return getPar( 25 ).value(); }
  double iKstarm       () const { return getPar( 26 ).value(); }
  double rKstarp       () const { return getPar( 27 ).value(); }
  double iKstarp       () const { return getPar( 28 ).value(); }
  double rRho          () const { return getPar( 29 ).value(); }
  double iRho          () const { return getPar( 30 ).value(); }
  double rOmega        () const { return getPar( 31 ).value(); }
  double iOmega        () const { return getPar( 32 ).value(); }
  double rF0_980       () const { return getPar( 33 ).value(); }
  double iF0_980       () const { return getPar( 34 ).value(); }
  double rF0_1370      () const { return getPar( 35 ).value(); }
  double iF0_1370      () const { return getPar( 36 ).value(); }
  double rF2_1270      () const { return getPar( 37 ).value(); }
  double iF2_1270      () const { return getPar( 38 ).value(); }
  double rK0starm_1430 () const { return getPar( 39 ).value(); }
  double iK0starm_1430 () const { return getPar( 40 ).value(); }
  double rK0starp_1430 () const { return getPar( 41 ).value(); }
  double iK0starp_1430 () const { return getPar( 42 ).value(); }
  double rK2starm_1430 () const { return getPar( 43 ).value(); }
  double iK2starm_1430 () const { return getPar( 44 ).value(); }
  double rK2starp_1430 () const { return getPar( 45 ).value(); }
  double iK2starp_1430 () const { return getPar( 46 ).value(); }
  double rSigma        () const { return getPar( 47 ).value(); }
  double iSigma        () const { return getPar( 48 ).value(); }
  double rSigma2       () const { return getPar( 49 ).value(); }
  double iSigma2       () const { return getPar( 50 ).value(); }
  double rKstarm_1680  () const { return getPar( 51 ).value(); }
  double iKstarm_1680  () const { return getPar( 52 ).value(); }
  double rNonReson     () const { return getPar( 53 ).value(); }
  double iNonReson     () const { return getPar( 54 ).value(); }

  // Utility functions.
  static double kallen( double m1, double m2, double m3 );
  static double m2ACmin( double m2AB );
  static double m2ACmax( double m2AB );

  std::complex< double > dalitzKsPiPi( double m2AB, double m2AC, double m2BC ) const;

  std::complex< double > evaluateReso( int index, double m2AB, double m2AC, double m2BC ) const;
  std::complex< double > coef        ( int index ) const;
  void cacheIntegralsMatrix( bool all = false ) const;
  void cacheDalitzIntegrals();

  //std::pair< Complex, Complex > matrixElements( int x, int y ) const;

  // < f | H | D^0 (t) > = 1/2 * [ ( 1 + \chi_f ) * A_f * e_1(gt) + ( 1 - \chi_f ) * A_f * e_2(gt) ]
  std::complex< double > e1( const double& gt ) const;
  std::complex< double > e2( const double& gt ) const;

public:
  DalitzD0mix( const std::vector< Variable >& vars, const std::vector< Parameter >& pars );

  static const double& mPi() { return _mPi; }
  static const double& mKs() { return _mKs; }
  static const double& mD0() { return _mD0; }

  void cache()
  {
    // Compute the terms of the matrices of integrals (Id and Ix) that need recalculation.
    if ( ! _propagatorsFixed )
      cacheIntegralsMatrix();

    // Compute the integrals I_1 and I_\chi if the amplitudes are constant.
    if ( ! _integralsFixed )
      cacheDalitzIntegrals();

    _norm = norm();
  }
  double evaluate()                                    const throw( PdfException );
  double evaluate( const std::vector< double >& vars ) const throw( PdfException );

  double norm() const;
};

#endif
