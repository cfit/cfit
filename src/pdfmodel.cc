
#include <vector>
#include <algorithm>
#include <functional>

#include <cfit/pdfexception.hh>
#include <cfit/pdfmodel.hh>
#include <cfit/functors.hh>

// Add a variable to the variables map.
void PdfModel::push( const Variable& var ) throw( PdfException )
{
  if ( _vars.find( var.name() ) != _vars.end() )
    throw PdfException( "Variable \"" + var.name() + "\" does already exist in this model." );

  _vars[ var.name() ] = var;
  _varOrder.push_back( var.name() );
}


// Add a parameter to the parameters map and add its name to the ordering vector.
void PdfModel::push( const Parameter& par ) throw( PdfException )
{
  if ( _pars.find( par.name() ) != _pars.end() )
    throw PdfException( "Parameter \"" + par.name() + "\" does already exist in this model." );

  _pars[ par.name() ] = par;
  _parOrder.push_back( par.name() );
}


// Retrieve the variable at specified index.
const Variable& PdfModel::getVar( int index ) const
{
  return _vars.find( _varOrder[ index ] )->second;
}

// Retrieve the parameter at specified index.
const Parameter& PdfModel::getPar( int index ) const
{
  return _pars.find( _parOrder[ index ] )->second;
}



// Set the variables to those given as argument.
// They must be sorted alphabetically.
void PdfModel::setVars( const std::vector< double >& vars ) throw( PdfException )
{
  if ( _vars.size() != vars.size() )
    throw PdfException( "Number of arguments passed does not match number of required arguments." );

  typedef std::map< std::string, Variable >::iterator vIter;
  int index = 0;
  for ( vIter var = _vars.begin(); var != _vars.end(); var++ )
    var->second.setValue( vars[ index++ ] );
}


// Set the parameters to those given as argument.
// They must be sorted alphabetically, since it's how MnUserParameters are passed
//    in the minimize function of the minimizers. It must be so, because pushing
//    two parameters with the same name would create confusion otherwise.
void PdfModel::setPars( const std::vector< double >& pars ) throw( PdfException )
{
  if ( _pars.size() != pars.size() )
    throw PdfException( "Number of arguments passed does not match number of required arguments." );

  typedef std::map< std::string, Parameter >::iterator pIter;
  int index = 0;
  for ( pIter par = _pars.begin(); par != _pars.end(); par++ )
    par->second.setValue( pars[ index++ ] );
}

