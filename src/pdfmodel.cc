
#include <vector>
#include <algorithm>
#include <functional>

#include <Minuit/FunctionMinimum.h>
#include <Minuit/MnUserParameters.h>
#include <Minuit/MinuitParameter.h>

#include <cfit/exceptions.hh>
#include <cfit/pdfmodel.hh>
#include <cfit/functors.hh>



// Add a variable to the variables map.
void PdfModel::push( const Variable& var )
{
  if ( _varMap.find( var.name() ) != _varMap.end() )
    throw PdfException( "Variable \"" + var.name() + "\" does already exist in this model." );

  _varMap[ var.name() ] = var;
  _varOrder.push_back( var.name() );
}


// Add a parameter to the parameters map and add its name to the ordering vector.
void PdfModel::push( const Parameter& par )
{
  if ( _parMap.find( par.name() ) != _parMap.end() )
    throw PdfException( "Parameter \"" + par.name() + "\" does already exist in this model." );

  _parMap[ par.name() ] = par;
  _parOrder.push_back( par.name() );
}


// Push the two parameters in a coefficient.
void PdfModel::push( const Coef& coef )
{
  push( coef.real() );
  push( coef.imag() );
}


// Push all the parameters in a parameter expression.
void PdfModel::push( const ParameterExpr& expr )
{
  const std::map< std::string, Parameter >& pars = expr.getPars();
  typedef std::map< const std::string, Parameter >::const_iterator pIter;
  for ( pIter par = pars.begin(); par != pars.end(); ++par )
    push( par->second );
}


// Push all the parameters in a coefficient expression.
void PdfModel::push( const CoefExpr& expr )
{
  const std::map< std::string, Parameter >& pars = expr.getPars();
  typedef std::map< const std::string, Parameter >::const_iterator pIter;
  for ( pIter par = pars.begin(); par != pars.end(); ++par )
    push( par->second );
}


// Push all the parameters in a resonance.
void PdfModel::push( const Resonance& reso )
{
  const std::map< std::string, Parameter >& pars = reso.getPars();
  typedef std::map< const std::string, Parameter >::const_iterator pIter;
  for ( pIter par = pars.begin(); par != pars.end(); ++par )
    push( par->second );
}


// Push all the parameters in an amplitude.
void PdfModel::push( const Amplitude& amp )
{
  const std::map< std::string, Parameter >& pars = amp.getPars();
  typedef std::map< const std::string, Parameter >::const_iterator pIter;
  for ( pIter par = pars.begin(); par != pars.end(); ++par )
    push( par->second );
}



// Retrieve the mapped variable with same name as argument.
const Variable& PdfModel::getVar( const Variable& var ) const
{
  return _varMap.find( var.name() )->second;
}

// Retrieve the variable at specified index.
const Variable& PdfModel::getVar( const int& idx ) const
{
  return _varMap.find( _varOrder[ idx ] )->second;
}

// Retrieve the mapped parameter with same name as argument.
const Parameter& PdfModel::getPar( const Parameter& par ) const
{
  return _parMap.find( par.name() )->second;
}

// Retrieve the parameter at specified index.
const Parameter& PdfModel::getPar( const int& idx ) const
{
  return _parMap.find( _parOrder[ idx ] )->second;
}



// Set the values of the parameter map from a vector of values, sorted alphabetically by parameter name.
// It is necessary that they have the same size.
void PdfModel::setParMap( const std::vector< double >& pars )
{
  if ( _parMap.size() != pars.size() )
    throw PdfException( "PdfModel::setParMap: Number of arguments passed does not match number of required arguments." );

  typedef std::map< std::string, Parameter >::iterator pIter;
  int index = 0;
  for ( pIter par = _parMap.begin(); par != _parMap.end(); ++par )
    par->second.setValue( pars[ index++ ] );
}


// Set the values of the parameter map from another parameter map. Their sizes may differ.
void PdfModel::setParMap( const std::map< std::string, Parameter >& pars )
{
  typedef std::map< const std::string, Parameter >::iterator pIter;
  for ( pIter par = _parMap.begin(); par != _parMap.end(); ++par )
    if ( pars.count( par->first ) )
      par->second.setValue( pars.at( par->first ).value() );
}



// Set the values of the parameter map from a minuit FunctionMinimum object. Their sizes may differ.
void PdfModel::setParMap( const FunctionMinimum& min )
{
  const std::vector< MinuitParameter >& pars = min.userParameters().parameters();

  // Do not require that the number of parameters in _parMap and min be the same.
  //    The user may want to run fit 1, then set fit result as initial values for
  //    fit 2, with some other parameters.

  // if ( _parMap.size() != pars.size() )
  //   throw PdfException( "PdfModel::setPars( minimum ): Number of arguments passed does not match number of required arguments." );

  typedef std::vector< MinuitParameter >::const_iterator pIter;
  for ( pIter par = pars.begin(); par != pars.end(); ++par )
    if ( _parMap.count( par->name() ) )
      _parMap[ par->name() ].set( par->value(), par->error() );
}


// Setter for individual parameter.
void PdfModel::setPar( const std::string& name, const double& val, const double& err ) throw( PdfException )
{
  if ( ! _parMap.count( name ) )
    throw PdfException( "Cannot set unexisting parameter " + name + "." );

  _parMap[ name ].set( val, err );
}


// Set the parameters to those given as argument.
// They must be sorted alphabetically, since it's how MnUserParameters are passed
//    in the minimize function of the minimizers. It must be so, because pushing
//    two parameters with the same name would create confusion otherwise.
void PdfModel::setPars( const std::vector< double >& pars ) throw( PdfException )
{
  setParMap( pars );

  setParExpr();
}

// The function must be virtual to allow the derived decay model classes to use their
//    own setPars function.
void PdfModel::setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException )
{
  setParMap( pars );

  setParExpr();
}

void PdfModel::setPars( const FunctionMinimum& min ) throw( PdfException )
{
  setParMap( min );

  setParExpr();
}


const double PdfModel::area( const double& min, const double& max ) const throw( PdfException )
{
  throw PdfException( "You are trying to find the area of a model that does not have this property." );
}
