#ifndef __MATRIX_HH__
#define __MATRIX_HH__

#include <iostream> // TEMPORAL

#include <cmath>
#include <vector>
#include <sstream>


template <class T> class Matrix
{
private:
  T** _mat;
  int _range;

  Matrix* min       ( int row, int col ) const;
  T*      operator[]( int row          ) const { return _mat[ row ]; }
public:
  Matrix( int range = 0 );
  Matrix( const Matrix< T >& mat );
  ~Matrix();

  inline void setRange( int range );

  T&          operator()( int row, int col ) const { return _mat[ row ][ col ]; }
  T           det       ()                   const;
  Matrix      inverse   ()                   const;
  std::string dump      ()                   const;

  const Matrix< T >& operator= ( const Matrix< T >&       right );
  const Matrix< T >& operator+=( const Matrix< T >&       right );
  const Matrix< T >& operator*=( const Matrix< T >&       right );
  const Matrix< T >& operator*=( const T&                 right );

  template <class M> friend const Matrix< M > operator*( const Matrix< M >& left, const Matrix< M >& right );
  template <class M> friend const Matrix< M > operator*( const Matrix< M >& left, const M&           right );
  template <class M> friend const Matrix< M > operator*( const M&           left, const Matrix< M >& right );

  template <class M>
  friend const Matrix< std::complex< M > > operator*( const Matrix< M >&       left, const std::complex< M >& right );
  template <class M>
  friend const Matrix< std::complex< M > > operator*( const std::complex< M >& left, const Matrix< M >&       right );

  template <class M> friend const Matrix< M > operator/( const Matrix< M >& left, const Matrix< M >& right );
  template <class M> friend const Matrix< M > operator/( const Matrix< M >& left, const M&           right );
  template <class M> friend const Matrix< M > operator/( const M&           left, const Matrix< M >& right );
};





template < class T >
Matrix< T >::Matrix( int range )
  : _range( range )
{
  // If range is specified, allocate memory for desired matrix.
  if ( range )
  {
    _mat = new T*[ range ];
    for ( int row = 0; row < range; ++row )
      _mat[ row ] = new T[ range ];
  }

  // Set the matrix elements to zero.
  for ( int row = 0; row < _range; row++ )
    for ( int col = 0; col < _range; col++ )
      _mat[ row ][ col ] = 0.;
};


template < class T >
Matrix< T >::Matrix( const Matrix< T >& mat )
  : _range( mat._range )
{
  // If the range is changed, delete any previous matrix stored
  //    and allocate elements with the newly specified range.
  if ( _range )
  {
    _mat = new T*[ _range ];
    for ( int row = 0; row < _range; ++row )
      _mat[ row ] = new T[ _range ];
  }

  // Set the matrix elements to the elements at the original matrix.
  for ( int row = 0; row < _range; row++ )
    for ( int col = 0; col < _range; col++ )
      _mat[ row ][ col ] = mat._mat[ row ][ col ];
}



template <class T> Matrix< T >::~Matrix()
{
  if ( _range )
  {
    for( int row = 0; row < _range; ++row )
      delete[] _mat[ row ];
    delete[] _mat;
  }
}



template <class T> inline void Matrix< T >::setRange( int range )
{
  // If the range is changed, delete any previous matrix stored
  //    and allocate elements with the newly specified range.
  if ( _range != range )
  {
    if ( _range )
    {
      for ( int row = 0; row < _range; ++row )
        delete[] _mat[ row ];
      delete[] _mat;
    }

    _mat = new T*[ range ];
    for ( int row = 0; row < range; ++row )
      _mat[ row ] = new T[ range ];

    // Set the new range.
    _range = range;
  }

  // Since user is willing to change the range, reset the matrix elements.
  for ( int row = 0; row < _range; row++ )
    for ( int col = 0; col < _range; col++ )
      _mat[ row ][ col ] = 0.;
}



template <class T> std::string Matrix< T >::dump() const
{
  std::ostringstream str;

  for ( int row = 0; row < _range; row++ )
  {
    str << "|";
    for ( int col = 0; col < _range; col++ )
      str << "\t" << _mat[ row ][ col ];
    str << "\t|" << std::endl;
  }

  return str.str();
}


template <class T> T Matrix< T >::det() const
{
  if ( _range == 1 )
    return _mat[ 0 ][ 0 ];

  // There's no need to define the range 2 determinant manually, but it may
  //    speed up the calculation.
  if ( _range == 2 )
    return _mat[ 0 ][ 0 ] * _mat[ 1 ][ 1 ] - _mat[ 0 ][ 1 ] * _mat[ 1 ][ 0 ];

  T sum = 0.;

  for ( int col = 0; col < _range; col++ )
  {
    Matrix< T >* minor = min( 0, col );
    sum += std::pow( -1., col ) * _mat[ 0 ][ col ] * minor->det();
    delete minor;
  }

  return sum;
}


// Returns the minor at (i, j).
template <class T> Matrix< T >* Matrix< T >::min( int row, int col ) const
{
  Matrix< T >* minor = new Matrix< T >();
  minor->setRange( _range - 1 );

  int minIndex = 0;

  for ( int r = 0; r < _range; r++ )
    for ( int c = 0; c < _range; c++ )
      if ( ( r != row ) && ( c != col ) )
      {
        (*minor)( minIndex / ( _range - 1 ), minIndex % ( _range - 1 ) ) = _mat[ r ][ c ];
        minIndex++;
      }

  return minor;
}


template <class T> Matrix< T > Matrix< T >::inverse() const
{
  Matrix< T > inv;
  inv.setRange( _range );

  T determinant = det();

  if ( determinant == 0. )
  {
    std::cerr << "This matrix has a null determinant and cannot be inverted. Returning zero matrix." << std::endl;
    for ( int row = 0; row < _range; row++ )
      for ( int col = 0; col < _range; col++ )
        inv( row, col ) = 0.;
    return inv;
  }

  for ( int row = 0; row < _range; row++ )
    for ( int col = 0; col < _range; col++ )
    {
      Matrix< T >* minor = min( row, col );
      inv._mat[col][row] = std::pow( -1., row + col ) * minor->det() / determinant;
      delete minor;
    }

  return inv;
}



template < class T >
const Matrix< T >& Matrix< T >::operator=( const Matrix< T >& mat )
{
  _range = mat._range;

  // Allocate elements with the given matrix range.
  if ( _range )
  {
    _mat = new T*[ _range ];
    for ( int row = 0; row < _range; ++row )
      _mat[ row ] = new T[ _range ];
  }

  // Set the matrix elements to the elements at the original matrix.
  for ( int row = 0; row < _range; row++ )
    for ( int col = 0; col < _range; col++ )
      _mat[ row ][ col ] = mat._mat[ row ][ col ];

  return *this;
}




template <class T>
const Matrix< T >& Matrix< T >::operator+=( const Matrix< T >& right )
{
  // Chech that the matrices have the correct range.
  if ( _range != right._range )
  {
    std::cerr << "These matrices cannot be added." << std::endl;
    return *this;
  }

  // Do the sum.
  for ( int row = 0; row < _range; row++ )
    for ( int col = 0; col < right._range; col++ )
      _mat[ row ][ col ] += right._mat[ row ][ col ];

  return *this;
}


template <class T>
const Matrix< T >& Matrix< T >::operator*=( const Matrix< T >& right )
{
  Matrix< T > left = *this;

  // Initialize the matrix.
  for ( int row = 0; row < left._range; row++ )
    for ( int col = 0; col < right._range; col++ )
      _mat[ row ][ col ] = 0.;

  // Do the product.
  for ( int row = 0; row < left._range; row++ )
    for ( int col = 0; col < right._range; col++ )
      for ( int line = 0; line < right._range; line++ )
	_mat[ row ][ col ] += left._mat[ row ][ line ] * right._mat[ line ][ col ];

  return *this;
}



template <class T>
const Matrix< T >& Matrix< T >::operator*=( const T& right )
{
  // Do the product.
  for ( int row = 0; row < _range; ++row )
    for ( int col = 0; col < _range; ++col )
      _mat[ row ][ col ] *= right;

  return *this;
}



template <class T>
const Matrix< T > operator*( const Matrix< T >& left, const Matrix< T >& right )
{
  Matrix< T > mat;

  // Chech that the matrices have the correct range.
  if ( left._range != right._range )
  {
    std::cerr << "These matrices cannot be multiplied." << std::endl;
    return mat;
  }

  mat.setRange( left._range );

  // Initialize the elements of the matrix.
  for ( int row = 0; row < left._range; row++ )
    for ( int col = 0; col < right._range; col++ )
      mat[ row ][ col ] = 0;

  // Do the product.
  for ( int row = 0; row < left._range; row++ )
    for ( int col = 0; col < right._range; col++ )
      for ( int line = 0; line < right._range; line++ )
	mat[ row ][ col ] += left._mat[ row ][ line ] * right._mat[ line ][ col ];

  return mat;
}


template <class T>
const Matrix< T > operator*( const Matrix< T >& left, const T& right )
{
  Matrix< T > mat;

  mat.setRange( left._range );

  for ( int row = 0; row < left._range; row++ )
    for ( int col = 0; col < left._range; col++ )
      mat[ row ][ col ] = left._mat[ row ][ col ] * right;

  return mat;
}


template <class T>
const Matrix< T > operator*( const T& left, const Matrix< T >& right )
{
  Matrix< T > mat;

  mat.setRange( right._range );

  for ( int row = 0; row < right._range; row++ )
    for ( int col = 0; col < right._range; col++ )
      mat[ row ][ col ] = left * right._mat[ row ][ col ];

  return mat;
}


template <class T>
const Matrix< std::complex< T > > operator*( const Matrix< T >& left, const std::complex< T >& right )
{
  Matrix< std::complex< T > > mat;

  mat.setRange( left._range );

  for ( int row = 0; row < left._range; row++ )
    for ( int col = 0; col < left._range; col++ )
      mat[ row ][ col ] = right * left._mat[ row ][ col ];

  return mat;
}


template <class T>
const Matrix< std::complex< T > > operator*( const std::complex< T >& left, const Matrix< T >& right )
{
  Matrix< std::complex< T > > mat;

  mat.setRange( right._range );

  for ( int row = 0; row < right._range; row++ )
    for ( int col = 0; col < right._range; col++ )
      mat[ row ][ col ] = left * right._mat[ row ][ col ];

  return mat;
}


template <class T>
const Matrix< T > operator/( const Matrix< T >& left, const Matrix< T >& right )
{
  Matrix< T > mat;

  // Chech that the matrices have the correct range.
  if ( left._range != right._range )
  {
    std::cerr << "These matrices cannot be multiplied." << std::endl;
    return mat;
  }

  mat.setRange( left._range );

  const Matrix< T >& invRight = right.inverse();

  // Initialize the elements of the matrix.
  for ( int row = 0; row < left._range; row++ )
    for ( int col = 0; col < right._range; col++ )
      mat[ row ][ col ] = 0;

  // Do the product.
  for ( int row = 0; row < left._range; row++ )
    for ( int col = 0; col < right._range; col++ )
      for ( int line = 0; line < right._range; line++ )
	mat[ row ][ col ] += left._mat[ row ][ line ] * invRight._mat[ line ][ col ];

  return mat;
}


template <class T>
const Matrix< T > operator/( const Matrix< T >& left, const T& right )
{
  Matrix< T > mat;

  mat.setRange( left._range );

  for ( int row = 0; row < left._range; row++ )
    for ( int col = 0; col < left._range; col++ )
      mat[ row ][ col ] = left._mat[ row ][ col ] / right;

  return mat;
}


template <class T>
const Matrix< T > operator/( const T& left, const Matrix< T >& right )
{
  Matrix< T > mat;

  mat.setRange( right._range );

  const Matrix< T >& invRight = right.inverse();

  for ( int row = 0; row < right._range; row++ )
    for ( int col = 0; col < right._range; col++ )
      mat[ row ][ col ] = left * invRight._mat[ row ][ col ];

  return mat;
}



#endif
