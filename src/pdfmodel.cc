
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
void PdfModel::push( const Variable& var ) throw( PdfException )
{
  if ( _varMap.find( var.name() ) != _varMap.end() )
    throw PdfException( "Variable \"" + var.name() + "\" does already exist in this model." );

  _varMap[ var.name() ] = var;
  _varOrder.push_back( var.name() );
}


// Add a parameter to the parameters map and add its name to the ordering vector.
void PdfModel::push( const Parameter& par ) throw( PdfException )
{
  if ( _parMap.find( par.name() ) != _parMap.end() )
    throw PdfException( "Parameter \"" + par.name() + "\" does already exist in this model." );

  _parMap[ par.name() ] = par;
  _parOrder.push_back( par.name() );
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
  if ( _parMap.size() != pars.size() )
    throw PdfException( "PdfModel::setPars( vector ): Number of arguments passed does not match number of required arguments." );

  typedef std::map< std::string, Parameter >::iterator pIter;
  int index = 0;
  for ( pIter par = _parMap.begin(); par != _parMap.end(); par++ )
    par->second.setValue( pars[ index++ ] );
}


void PdfModel::setPars( const FunctionMinimum& min ) throw( PdfException )
{
  const MnUserParameters& pars = min.userParameters();
  const std::vector< MinuitParameter >& parVec = pars.parameters();

  if ( _parMap.size() != parVec.size() )
    throw PdfException( "PdfModel::setPars( minimum ): Number of arguments passed does not match number of required arguments." );

  typedef std::vector< MinuitParameter >::const_iterator pIter;
  for ( pIter par = parVec.begin(); par != parVec.end(); ++par )
    _parMap[ par->name() ].set( par->value(), par->error() );
}


// The function must be virtual to allow the derived decay model classes to use their
//    own setPars function.
void PdfModel::setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException )
{
  typedef std::map< const std::string, Parameter >::iterator pIter;
  for ( pIter par = _parMap.begin(); par != _parMap.end(); ++par )
    par->second.setValue( pars.find( par->first )->second.value() );
}


const double PdfModel::area( const double& min, const double& max ) const throw( PdfException )
{
  throw PdfException( "You are trying to find the area of a model that does not have this property." );
}
