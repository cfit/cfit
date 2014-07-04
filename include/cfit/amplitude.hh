#ifndef __AMPLITUDE_HH__
#define __AMPLITUDE_HH__

#include <vector>
#include <string>
#include <complex>
#include <map>

#include <cfit/operation.hh>
#include <cfit/exceptions.hh>
#include <cfit/coef.hh>
#include <cfit/coefexpr.hh>
#include <cfit/resonance.hh>
#include <cfit/fvector.hh>

class PhaseSpace;

class Amplitude
{
private:
  std::map< std::string, Parameter > _parMap;

  std::vector< std::complex< double > > _ctnts;
  std::vector< Parameter              > _parms;
  std::vector< Coef                   > _coefs;
  std::vector< Resonance*             > _resos;
  std::vector< Fvector                > _fvecs;
  std::vector< Operation::Op          > _opers;
  std::string                           _expression;

  // Clean up the content of the Amplitude containers.
  void clear();

  void append( const double&                 ctnt );
  void append( const std::complex< double >& ctnt );
  void append( const Parameter&              parm );
  void append( const Coef&                   coef );
  void append( const CoefExpr&               expr );
  void append( const Resonance&              reso );
  void append( const Fvector&                fvec );
  void append( const Amplitude&              ampl );
  void append( const Operation::Op&          oper );

  // Constructor to be called by binary operators.
  template< class L, class R >
  Amplitude( const L& left, const R& right, const Operation::Op& oper )
  {
    append( left  );
    append( right );
    append( oper  );
  }

public:
  Amplitude() {};

  Amplitude( const double& ctnt )
  {
    append( ctnt );
  }

  Amplitude( const std::complex< double >& ctnt )
  {
    append( ctnt );
  }

  Amplitude( const Coef& coef )
  {
    append( coef );
  }

  Amplitude( const CoefExpr& expr )
  {
    append( expr );
  }

  Amplitude( const Resonance& reso )
  {
    append( reso );
  }

  Amplitude( const Fvector& vec )
  {
    append( vec );
  }

  Amplitude( const Amplitude& amp )
  {
    append( amp );
  }

  ~Amplitude()
  {
    // Delete all the pointers to resonance, since they have been allocated
    //    by the copy() function of each resonance.
    typedef std::vector< Resonance* >::iterator rIter;
    for ( rIter res = _resos.begin(); res != _resos.end(); ++res )
      delete *res;
  }

  void useHelicity( const bool helicity = true );
  void useTwoBW   ( const bool twoBW    = true );

  const std::map< std::string, Parameter >& getPars() const { return _parMap; };
  void setPars( const std::map< std::string, Parameter >& pars );

  std::complex< double > evaluate( const PhaseSpace& ps,
				   const double&     mSq12,
				   const double&     mSq13,
				   const double&     mSq23 ) const throw( PdfException );

  // Assignment operations.
  const Amplitude& operator= ( const double&                 ctnt );
  const Amplitude& operator= ( const std::complex< double >& ctnt );
  const Amplitude& operator= ( const Coef&                   coef );
  const Amplitude& operator= ( const CoefExpr&               expr );

  const Amplitude& operator= ( const Resonance&              right );
  const Amplitude& operator= ( const Fvector&                right );
  const Amplitude& operator= ( const Amplitude&              right );

  const Amplitude& operator+=( const double&                 right );
  const Amplitude& operator-=( const double&                 right );
  const Amplitude& operator*=( const double&                 right );
  const Amplitude& operator/=( const double&                 right );

  const Amplitude& operator+=( const std::complex< double >& right );
  const Amplitude& operator-=( const std::complex< double >& right );
  const Amplitude& operator*=( const std::complex< double >& right );
  const Amplitude& operator/=( const std::complex< double >& right );

  const Amplitude& operator+=( const Coef&                   right );
  const Amplitude& operator-=( const Coef&                   right );
  const Amplitude& operator*=( const Coef&                   right );
  const Amplitude& operator/=( const Coef&                   right );

  const Amplitude& operator+=( const CoefExpr&               right );
  const Amplitude& operator-=( const CoefExpr&               right );
  const Amplitude& operator*=( const CoefExpr&               right );
  const Amplitude& operator/=( const CoefExpr&               right );

  const Amplitude& operator+=( const Resonance&              right );
  const Amplitude& operator-=( const Resonance&              right );

  const Amplitude& operator+=( const Fvector&                right );
  const Amplitude& operator-=( const Fvector&                right );

  const Amplitude& operator+=( const Amplitude&              right );
  const Amplitude& operator-=( const Amplitude&              right );

  friend const Amplitude operator+( const Resonance& left, const Resonance& right );
  friend const Amplitude operator-( const Resonance& left, const Resonance& right );

  friend const Amplitude operator+( const Fvector&   left, const Fvector&   right );
  friend const Amplitude operator-( const Fvector&   left, const Fvector&   right );

  friend const Amplitude operator+( const Amplitude& left, const Amplitude& right );
  friend const Amplitude operator-( const Amplitude& left, const Amplitude& right );


  // Operations of constants and resonances.
  friend const Amplitude operator+( const double&    left, const Resonance& right );
  friend const Amplitude operator-( const double&    left, const Resonance& right );
  friend const Amplitude operator*( const double&    left, const Resonance& right );

  friend const Amplitude operator+( const Resonance& left, const double&    right );
  friend const Amplitude operator-( const Resonance& left, const double&    right );
  friend const Amplitude operator*( const Resonance& left, const double&    right );
  friend const Amplitude operator/( const Resonance& left, const double&    right );


  // Operations with constants and amplitudes.
  friend const Amplitude operator+( const double&    left, const Amplitude& right );
  friend const Amplitude operator-( const double&    left, const Amplitude& right );
  friend const Amplitude operator*( const double&    left, const Amplitude& right );

  friend const Amplitude operator+( const Amplitude& left, const double&    right );
  friend const Amplitude operator-( const Amplitude& left, const double&    right );
  friend const Amplitude operator*( const Amplitude& left, const double&    right );
  friend const Amplitude operator/( const Amplitude& left, const double&    right );


  // Operations of coefficients and resonances.
  friend const Amplitude operator+( const Coef&      left, const Resonance& right );
  friend const Amplitude operator-( const Coef&      left, const Resonance& right );
  friend const Amplitude operator*( const Coef&      left, const Resonance& right );

  friend const Amplitude operator+( const Resonance& left, const Coef&      right );
  friend const Amplitude operator-( const Resonance& left, const Coef&      right );
  friend const Amplitude operator*( const Resonance& left, const Coef&      right );
  friend const Amplitude operator/( const Resonance& left, const Coef&      right );


  // Operations with coefficients and amplitudes.
  friend const Amplitude operator+( const Coef&      left, const Amplitude& right );
  friend const Amplitude operator-( const Coef&      left, const Amplitude& right );
  friend const Amplitude operator*( const Coef&      left, const Amplitude& right );

  friend const Amplitude operator+( const Amplitude& left, const Coef&      right );
  friend const Amplitude operator-( const Amplitude& left, const Coef&      right );
  friend const Amplitude operator*( const Amplitude& left, const Coef&      right );
  friend const Amplitude operator/( const Amplitude& left, const Coef&      right );


  // Operations of coefficient expressions and resonances.
  friend const Amplitude operator+( const CoefExpr&  left, const Resonance& right );
  friend const Amplitude operator-( const CoefExpr&  left, const Resonance& right );
  friend const Amplitude operator*( const CoefExpr&  left, const Resonance& right );

  friend const Amplitude operator+( const Resonance& left, const CoefExpr&  right );
  friend const Amplitude operator-( const Resonance& left, const CoefExpr&  right );
  friend const Amplitude operator*( const Resonance& left, const CoefExpr&  right );
  friend const Amplitude operator/( const Resonance& left, const CoefExpr&  right );


  // Operations with coefficient expressions and amplitudes.
  friend const Amplitude operator+( const CoefExpr&  left, const Amplitude& right );
  friend const Amplitude operator-( const CoefExpr&  left, const Amplitude& right );
  friend const Amplitude operator*( const CoefExpr&  left, const Amplitude& right );

  friend const Amplitude operator+( const Amplitude& left, const CoefExpr&  right );
  friend const Amplitude operator-( const Amplitude& left, const CoefExpr&  right );
  friend const Amplitude operator*( const Amplitude& left, const CoefExpr&  right );
  friend const Amplitude operator/( const Amplitude& left, const CoefExpr&  right );


  // Operations with resonances and amplitudes.
  friend const Amplitude operator+( const Resonance& left, const Amplitude& right );
  friend const Amplitude operator-( const Resonance& left, const Amplitude& right );

  friend const Amplitude operator+( const Amplitude& left, const Resonance& right );
  friend const Amplitude operator-( const Amplitude& left, const Resonance& right );


  // Operations with Fvector components and amplitudes.
  friend const Amplitude operator+( const Fvector&   left, const Amplitude& right );
  friend const Amplitude operator-( const Fvector&   left, const Amplitude& right );

  friend const Amplitude operator+( const Amplitude& left, const Fvector&   right );
  friend const Amplitude operator-( const Amplitude& left, const Fvector&   right );
};

#endif
