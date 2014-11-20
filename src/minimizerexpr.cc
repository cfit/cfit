
#include <iostream>

#include <algorithm>

#include <Minuit/MnMigrad.h>

#include <cfit/functors.hh>
#include <cfit/minimizerexpr.hh>


double MinimizerExpr::up() const throw( MinimizerException )
{
  if ( _up < 0.0 )
    throw MinimizerException( "The minimizer variation that specifies the sigma level of uncertainties must be positive." );

  return _up;
}


double MinimizerExpr::operator()( const std::vector< double >& pars ) const throw( PdfException )
{
  if ( _minimizers.empty() )
    throw PdfException( "Minimizer expression does not contain any minimizer." );

  if ( pars.size() != _parMap.size() )
    throw PdfException( "Number of parameters passed does not match number of required arguments." );

  // Set the local values of the parameters to match them with their names.
  //    Unfortunately, _parMap cannot be modified because this function MUST be const
  //    (it is a minuit requirement). A non-member parMap has to be built.
  std::map< std::string, Parameter > parMap;
  typedef std::map< std::string, Parameter >::const_iterator pIter;
  int index = 0;
  for ( pIter par = _parMap.begin(); par != _parMap.end(); ++par )
    parMap[ par->first ].setValue( pars[ index++ ] );

  // Set each minimizer's parameter vector from the just set local values of the parameters.
  typedef std::vector< const Minimizer* >::const_iterator mIter;
  double total = 0.;
  for ( mIter mmzr = _minimizers.begin(); mmzr != _minimizers.end(); ++mmzr )
  {
    const std::map< std::string, Parameter >& mParMap = (*mmzr)->pdf().getPars();
    std::vector< double > mPars;
    for ( pIter par = mParMap.begin(); par != mParMap.end(); ++par )
      mPars.push_back( parMap.find( par->first )->second.value() );

    total += (**mmzr)( mPars );
  }

  if ( _verbose )
    std::cout << "total = " << total << std::endl;

  return total;
}



FunctionMinimum MinimizerExpr::minimize() const
{
  // Work with Minuit user defined parameters.
  MnUserParameters upar;

  typedef std::map< std::string, Parameter >::const_iterator pIter;

  // Set the Minuit parameters' name, value and uncertainty.
  for ( pIter par = _parMap.begin(); par != _parMap.end(); ++par )
  {
    upar.add( par->first.c_str(), par->second.value(), par->second.error() );

    // Fix the parameters that are set to be fixed.
    if ( par->second.isFixed() )
      upar.fix( par->first.c_str() );

    // Set the blinding if requested.
    if ( par->second.isBlind() )
      upar.blind( par->first.c_str() );

    // Set the limits for those parameters that have some.
    if ( par->second.hasLimits() )
      upar.setLimits( par->first.c_str(), par->second.lower(), par->second.upper() );
  }

  MnMigrad migrad( *this, upar );

  return migrad();
}



MinimizerExpr& MinimizerExpr::operator=( const Minimizer& right )
{
  // If this MinimizerExpr already contains some stuff, it should not be deallocated,
  //    since MinimizerExpr does not own any of the minimizers it uses.

  // Just assign everything from the given minimizer information.

  _up = right.up();

  // Append all the parameters to the local _parMap container.
  _parMap = right.pdf().getPars();

  // Append the given minimizer.
  _minimizers.clear();
  _minimizers.push_back( &right );

  return *this;
}



MinimizerExpr& MinimizerExpr::operator+=( const Minimizer& right )
{
  if ( ( ! _minimizers.empty() ) && ( _up != right.up() ) )
    throw PdfException( "Cannot add two minimizers that do not have a common up value." );

  // If the expression is being initialized right here, set its up value.
  if ( _minimizers.empty() )
    _up = right.up();

  const std::map< std::string, Parameter >& rPars = right.pdf().getPars();

  // Append all the parameters to the local _parMap container.
  _parMap.insert( rPars.begin(), rPars.end() );

  // Append the given minimizer.
  _minimizers.push_back( &right );

  return *this;
}


MinimizerExpr& MinimizerExpr::operator+=( const MinimizerExpr& right )
{
  if ( ( ! _minimizers.empty() ) && ( _up != right._up ) )
    throw PdfException( "Cannot add two minimizers that do not have a common up value." );

  // If the expression is being initialized right here, set its up value.
  if ( _minimizers.empty() )
    _up = right.up();

  // Append all the parameters to the local _parMap container.
  _parMap.insert( right._parMap.begin(), right._parMap.end() );

  // Append all the given minimizers.
  _minimizers.insert( _minimizers.end(), right._minimizers.begin(), right._minimizers.end() );

  return *this;
}



MinimizerExpr operator+( const Minimizer& left, const Minimizer& right )
{
  if ( left.up() != right.up() )
    throw PdfException( "Cannot add two minimizers that do not have a common up value." );

  // Minimizer expression to be returned.
  MinimizerExpr total;

  // If both left and right minimizers share the same up value, this is what the
  //    minimizer expression will use.
  total._up = left.up();

  // Find the parameters that the given arguments depend on.
  const std::map< std::string, Parameter >& lPars = left .pdf().getPars();
  const std::map< std::string, Parameter >& rPars = right.pdf().getPars();

  // Append all the parameters to the local _parMap container.
  total._parMap.insert( lPars.begin(), lPars.end() );
  total._parMap.insert( rPars.begin(), rPars.end() );

  // Append the two given minimizers.
  total._minimizers.push_back( &left  );
  total._minimizers.push_back( &right );

  return total;
}

