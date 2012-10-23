
#include <algorithm>
#include <functional>

#include <cfit/models/polynomial.hh>
#include <cfit/math.hh>


Polynomial::Polynomial( const Variable& x, const Parameter& c1 )
  : _hasLower( false ), _hasUpper( false ), _lower( 0.0 ), _upper( 0.0 )
{
  push( x );

  push( c1 );
}


Polynomial::Polynomial( const Variable& x, const Parameter& c1, const Parameter& c2 )
  : _hasLower( false ), _hasUpper( false ), _lower( 0.0 ), _upper( 0.0 )
{
  push( x );

  push( c1 );
  push( c2 );
}


Polynomial::Polynomial( const Variable& x, const std::vector< Parameter >& coefs )
  : _hasLower( false ), _hasUpper( false ), _lower( 0.0 ), _upper( 0.0 )
{
  push( x );

//   std::for_each( coefs.begin(), coefs.end(), std::bind1st( std::mem_fun_ref( &Polynomial::push ), this ) );

  typedef std::vector< Parameter >::const_iterator pIter;
  for ( pIter par = coefs.begin(); par != coefs.end(); ++par )
    push( *par );
}



Polynomial* Polynomial::copy() const
{
  return new Polynomial( *this );
}


double Polynomial::coef( const unsigned& index )  const
{
  return getPar( index ).value();
}


void Polynomial::setLowerLimit( const double& lower )
{
  _hasLower = true;
  _lower    = lower;

  // Run cache if both upper and lower limits are defined.
  if ( _hasUpper )
    cache();
}


void Polynomial::setUpperLimit( const double& upper )
{
  _hasUpper = true;
  _upper    = upper;

  // Run cache if both upper and lower limits are defined.
  if ( _hasLower )
    cache();
}


void Polynomial::setLimits( const double& lower, const double& upper )
{
  _hasLower = true;
  _hasUpper = true;
  _lower    = lower;
  _upper    = upper;

  cache();
}


void Polynomial::cache()
{
  if ( ! _hasLower || ! _hasUpper )
    throw PdfException( "Cannot evaluate polynomial without upper and lower limits defined." );

  unsigned order = _parOrder.size();

  _norm = ( _upper - _lower );
  for ( unsigned ord = 0; ord < order; ++ord )
    _norm += coef( ord ) * ( std::pow( _upper, ord + 2 ) - std::pow( _lower, ord + 2 ) ) / ( ord + 2 );
}


const double Polynomial::evaluate( const double& x ) const throw( PdfException )
{
  if ( ! _hasLower || ! _hasUpper )
    throw PdfException( "Cannot evaluate polynomial without upper and lower limits defined." );

  unsigned order = _parOrder.size();

  // ( 1 + Sum( x^k ) ) / norm.
  double val = 1.0; // Constant term.
  for ( unsigned ord = 0; ord < order; ++ord )
    val += coef( ord ) * std::pow( x, ord + 1 );

  return val / _norm;
}


const double Polynomial::evaluate() const throw( PdfException )
{
  return evaluate( getVar( 0 ).value() );
}


const double Polynomial::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  return evaluate( vars[ 0 ] );
}
