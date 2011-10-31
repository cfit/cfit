#ifndef __AMPLITUDE_HH__
#define __AMPLITUDE_HH__

#include <vector>
#include <string>
#include <complex>
#include <map>

#include <cfit/operation.hh>
#include <cfit/pdfexception.hh>
#include <cfit/coef.hh>
#include <cfit/resonance.hh>

class PhaseSpace;

class Amplitude
{
private:
  std::vector< std::complex< double > > _ctnts;
  std::vector< Coef                   > _coefs;
  std::vector< Resonance*             > _resos;
  std::vector< Operation::Op          > _opers;
  std::string                           _expression;

  void append( const double&                 ctnt );
  void append( const std::complex< double >& ctnt );
  void append( const Coef&                   coef );
  void append( const Resonance&              reso );
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

  std::complex< double > operate( const std::complex< double >& x,
				  const std::complex< double >& y,
				  const Operation::Op&          oper ) const throw( PdfException );
public:
  Amplitude() {};
  Amplitude( const Resonance& reso )
  {
    append( reso );
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

  void setPars( const std::map< std::string, Parameter >& pars );

  std::complex< double > evaluate( const PhaseSpace& ps,
				   const double&     mSq12,
				   const double&     mSq13,
				   const double&     mSq23 ) const throw( PdfException );

  // Assignment operations.
  const Amplitude& operator= ( const Resonance& right );
  const Amplitude& operator= ( const Amplitude& right );

  const Amplitude& operator+=( const double&    right );
  const Amplitude& operator-=( const double&    right );
  const Amplitude& operator*=( const double&    right );
  const Amplitude& operator/=( const double&    right );

  const Amplitude& operator+=( const Coef&      right );
  const Amplitude& operator-=( const Coef&      right );
  const Amplitude& operator*=( const Coef&      right );
  const Amplitude& operator/=( const Coef&      right );

  const Amplitude& operator+=( const Resonance& right );
  const Amplitude& operator-=( const Resonance& right );

  const Amplitude& operator+=( const Amplitude& right );
  const Amplitude& operator-=( const Amplitude& right );

  // Operations of objects with themselves.
  friend const Amplitude operator+( const Coef&      left, const Coef&      right );
  friend const Amplitude operator-( const Coef&      left, const Coef&      right );
  friend const Amplitude operator*( const Coef&      left, const Coef&      right );
  friend const Amplitude operator/( const Coef&      left, const Coef&      right );

  friend const Amplitude operator+( const Resonance& left, const Resonance& right );
  friend const Amplitude operator-( const Resonance& left, const Resonance& right );

  friend const Amplitude operator+( const Amplitude& left, const Amplitude& right );
  friend const Amplitude operator-( const Amplitude& left, const Amplitude& right );


  // Operations of constants and coefficients.
  friend const Amplitude operator+( const double&    left, const Coef&      right );
  friend const Amplitude operator-( const double&    left, const Coef&      right );
  friend const Amplitude operator*( const double&    left, const Coef&      right );
  friend const Amplitude operator/( const double&    left, const Coef&      right );

  friend const Amplitude operator+( const Coef&      left, const double&    right );
  friend const Amplitude operator-( const Coef&      left, const double&    right );
  friend const Amplitude operator*( const Coef&      left, const double&    right );
  friend const Amplitude operator/( const Coef&      left, const double&    right );


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


  // Operations with resonances and amplitudes.
  friend const Amplitude operator+( const Resonance& left, const Amplitude& right );
  friend const Amplitude operator-( const Resonance& left, const Amplitude& right );

  friend const Amplitude operator+( const Amplitude& left, const Resonance& right );
  friend const Amplitude operator-( const Amplitude& left, const Resonance& right );
};

#endif
