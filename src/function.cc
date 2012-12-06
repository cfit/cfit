
#include <sstream>
#include <vector>
#include <stack>
#include <cmath>

#include <algorithm>

#include <Minuit/FunctionMinimum.h>

#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/parameterexpr.hh>

#include <cfit/function.hh>
#include <cfit/operation.hh>


void Function::clear()
{
  _varMap.clear();
  _parMap.clear();

  _expression.clear();

  _opers.clear();
  _ctnts.clear();
  _varbs.clear();
  _parms.clear();
}


// Append a constant.
void Function::append( const double& ctnt )
{
  _ctnts.push_back( ctnt );

  _expression += "c"; // c = constant.
}


// Append a variable.
void Function::append( const Variable& var )
{
  _varMap[ var.name() ] = var;
  _varbs.push_back( var.name() );

  _expression += "v"; // v = variable.
}


// Append a parameter.
void Function::append( const Parameter& par )
{
  _parMap[ par.name() ] = par;
  _parms.push_back( par.name() );

  _expression += "p"; // p = parameter.
}


// Append a parameter expression.
void Function::append( const ParameterExpr& expr )
{
  for ( std::vector< Parameter >::const_iterator par = expr._parms.begin(); par != expr._parms.end(); ++par )
  {
    _parMap[ par->name() ] = *par;
    _parms.push_back( par->name() );
  }

  _opers.insert( _opers.end(), expr._opers.begin(), expr._opers.end() );
  _ctnts.insert( _ctnts.end(), expr._ctnts.begin(), expr._ctnts.end() );

  _expression += expr._expression;
}


// Append a pdf expression.
void Function::append( const Function& func )
{
  _varMap.insert(               func._varMap.begin(), func._varMap.end() );
  _parMap.insert(               func._parMap.begin(), func._parMap.end() );

  _opers .insert( _opers.end(), func._opers .begin(), func._opers .end() );
  _ctnts .insert( _ctnts.end(), func._ctnts .begin(), func._ctnts .end() );
  _varbs .insert( _varbs.end(), func._varbs .begin(), func._varbs .end() );
  _parms .insert( _parms.end(), func._parms .begin(), func._parms .end() );

  _expression += func._expression;
}



// Setter for individual variable.
void Function::setVar( const std::string& name, const double& val, const double& err ) throw( PdfException )
{
  if ( ! _varMap.count( name ) )
    throw PdfException( "Cannot set unexisting variable " + name + "." );

  _varMap[ name ].set( val, err );
}

// Setter for individual parameter.
void Function::setPar( const std::string& name, const double& val, const double& err ) throw( PdfException )
{
  if ( ! _parMap.count( name ) )
    throw PdfException( "Cannot set unexisting parameter " + name + "." );

  _parMap[ name ].set( val, err );
}



void Function::setVars( const std::map< std::string, Variable >& vars ) throw( PdfException )
{
  typedef std::map< std::string, Variable >::iterator vIter;
  for ( vIter var = _varMap.begin(); var != _varMap.end(); ++var )
    var->second.setValue( vars.find( var->first )->second.value() );
}


void Function::setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException )
{
  typedef std::map< std::string, Parameter >::iterator pIter;
  for ( pIter par = _parMap.begin(); par != _parMap.end(); ++par )
    par->second.setValue( pars.find( par->first )->second.value() );
}


void Function::setPars( const FunctionMinimum& min ) throw( PdfException )
{
  const MnUserParameters& pars = min.userParameters();
  const std::vector< MinuitParameter >& parVec = pars.parameters();

  typedef std::vector< MinuitParameter >::const_iterator pIter;
  for ( pIter par = parVec.begin(); par != parVec.end(); ++par )
    _parMap[ par->name() ].set( par->value(), par->error() );
}



// Before running this function, the Function::setVars( vars ) function must be called.
//    To avoid the risk of forgetting it, run Function::evaluate( vars ).
double Function::evaluate() const throw( PdfException )
{
  std::stack< double > values;

  double x;
  double y;

  std::vector< Operation::Op >::const_iterator ops = _opers.begin();
  std::vector< double        >::const_iterator ctt = _ctnts.begin();
  std::vector< std::string   >::const_iterator var = _varbs.begin();
  std::vector< std::string   >::const_iterator par = _parms.begin();

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
    if ( *ch == 'c' )
      values.push( *ctt++ );
    else if ( *ch == 'v' )
      values.push( _varMap.find( *var++ )->second.value() );
    else if ( *ch == 'p' )
      values.push( _parMap.find( *par++ )->second.value() );
    else
    {
      if ( *ch == 'b' )
      {
        if ( values.size() < 2 )
          throw PdfException( "Parse error: not enough values in the stack." );
        y = values.top();
        values.pop();
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, y, *ops++ ) );
      }
      else if ( *ch == 'u' )
      {
        if ( values.empty() )
          throw PdfException( "Parse error: not enough values in the stack." );
        x = values.top();
        values.pop();
        values.push( Operation::operate( x, *ops++ ) );
      }
      else
        throw PdfException( std::string( "Parse error: unknown operation " ) + *ch + "." );
    }

  if ( values.size() != 1 )
    throw PdfException( "Function parse error: too many values have been supplied." );

  return values.top();
}


const Function pow( const Function& left, const double& right )
{
  return Function( left, right, Operation::pow );
}

const Function pow( const Variable& left, const double& right )
{
  return Function( left, right, Operation::pow );
}


// Operations with variables with themselves.
const Function operator+( const Variable& left, const Variable& right )
{
  return Function( left, right, Operation::plus );
}

const Function operator-( const Variable& left, const Variable& right )
{
  return Function( left, right, Operation::minus );
}

const Function operator*( const Variable& left, const Variable& right )
{
  return Function( left, right, Operation::mult );
}

const Function operator/( const Variable& left, const Variable& right )
{
  return Function( left, right, Operation::div );
}


// Operations with variables as the left operand.
const Function operator+( const Variable& left, const double& right )
{
  return Function( left, right, Operation::plus );
}

const Function operator-( const Variable& left, const double& right )
{
  return Function( left, right, Operation::minus );
}

const Function operator*( const Variable& left, const double& right )
{
  return Function( left, right, Operation::mult );
}

const Function operator/( const Variable& left, const double& right )
{
  return Function( left, right, Operation::div );
}



const Function operator+( const Variable& left, const Parameter& right )
{
  return Function( left, right, Operation::plus );
}

const Function operator-( const Variable& left, const Parameter& right )
{
  return Function( left, right, Operation::minus );
}

const Function operator*( const Variable& left, const Parameter& right )
{
  return Function( left, right, Operation::mult );
}

const Function operator/( const Variable& left, const Parameter& right )
{
  return Function( left, right, Operation::div );
}



const Function operator+( const Variable& left, const ParameterExpr& right )
{
  return Function( left, right, Operation::plus );
}

const Function operator-( const Variable& left, const ParameterExpr& right )
{
  return Function( left, right, Operation::minus );
}

const Function operator*( const Variable& left, const ParameterExpr& right )
{
  return Function( left, right, Operation::mult );
}

const Function operator/( const Variable& left, const ParameterExpr& right )
{
  return Function( left, right, Operation::div );
}


// Operations with variables as the right operand.
const Function operator+( const double& left, const Variable& right )
{
  return Function( left, right, Operation::plus );
}

const Function operator-( const double& left, const Variable& right )
{
  return Function( left, right, Operation::minus );
}

const Function operator*( const double& left, const Variable& right )
{
  return Function( left, right, Operation::mult );
}

const Function operator/( const double& left, const Variable& right )
{
  return Function( left, right, Operation::div );
}



const Function operator+( const Parameter& left, const Variable& right )
{
  return Function( left, right, Operation::plus );
}

const Function operator-( const Parameter& left, const Variable& right )
{
  return Function( left, right, Operation::minus );
}

const Function operator*( const Parameter& left, const Variable& right )
{
  return Function( left, right, Operation::mult );
}

const Function operator/( const Parameter& left, const Variable& right )
{
  return Function( left, right, Operation::div );
}



const Function operator+( const ParameterExpr& left, const Variable& right )
{
  return Function( left, right, Operation::plus );
}

const Function operator-( const ParameterExpr& left, const Variable& right )
{
  return Function( left, right, Operation::minus );
}

const Function operator*( const ParameterExpr& left, const Variable& right )
{
  return Function( left, right, Operation::mult );
}

const Function operator/( const ParameterExpr& left, const Variable& right )
{
  return Function( left, right, Operation::div );
}



// Operations with functions with themselves.
const Function operator+( const Function& left, const Function& right )
{
  return Function( left, right, Operation::plus );
}

const Function operator-( const Function& left, const Function& right )
{
  return Function( left, right, Operation::minus );
}

const Function operator*( const Function& left, const Function& right )
{
  return Function( left, right, Operation::mult );
}

const Function operator/( const Function& left, const Function& right )
{
  return Function( left, right, Operation::div );
}


// Operations with functions as the left operand.
const Function operator+( const Function& left, const double& right )
{
  return Function( left, right, Operation::plus );
}

const Function operator-( const Function& left, const double& right )
{
  return Function( left, right, Operation::minus );
}

const Function operator*( const Function& left, const double& right )
{
  return Function( left, right, Operation::mult );
}

const Function operator/( const Function& left, const double& right )
{
  return Function( left, right, Operation::div );
}



const Function operator+( const Function& left, const Variable& right )
{
  return Function( left, right, Operation::plus );
}

const Function operator-( const Function& left, const Variable& right )
{
  return Function( left, right, Operation::minus );
}

const Function operator*( const Function& left, const Variable& right )
{
  return Function( left, right, Operation::mult );
}

const Function operator/( const Function& left, const Variable& right )
{
  return Function( left, right, Operation::div );
}



const Function operator+( const Function& left, const Parameter& right )
{
  return Function( left, right, Operation::plus );
}

const Function operator-( const Function& left, const Parameter& right )
{
  return Function( left, right, Operation::minus );
}

const Function operator*( const Function& left, const Parameter& right )
{
  return Function( left, right, Operation::mult );
}

const Function operator/( const Function& left, const Parameter& right )
{
  return Function( left, right, Operation::div );
}



const Function operator+( const Function& left, const ParameterExpr& right )
{
  return Function( left, right, Operation::plus );
}

const Function operator-( const Function& left, const ParameterExpr& right )
{
  return Function( left, right, Operation::minus );
}

const Function operator*( const Function& left, const ParameterExpr& right )
{
  return Function( left, right, Operation::mult );
}

const Function operator/( const Function& left, const ParameterExpr& right )
{
  return Function( left, right, Operation::div );
}


// Operations with functions as the right operand.
const Function operator+( const double& left, const Function& right )
{
  return Function( left, right, Operation::plus );
}

const Function operator-( const double& left, const Function& right )
{
  return Function( left, right, Operation::minus );
}

const Function operator*( const double& left, const Function& right )
{
  return Function( left, right, Operation::mult );
}

const Function operator/( const double& left, const Function& right )
{
  return Function( left, right, Operation::div );
}



const Function operator+( const Variable& left, const Function& right )
{
  return Function( left, right, Operation::plus );
}

const Function operator-( const Variable& left, const Function& right )
{
  return Function( left, right, Operation::minus );
}

const Function operator*( const Variable& left, const Function& right )
{
  return Function( left, right, Operation::mult );
}

const Function operator/( const Variable& left, const Function& right )
{
  return Function( left, right, Operation::div );
}



const Function operator+( const Parameter& left, const Function& right )
{
  return Function( left, right, Operation::plus );
}

const Function operator-( const Parameter& left, const Function& right )
{
  return Function( left, right, Operation::minus );
}

const Function operator*( const Parameter& left, const Function& right )
{
  return Function( left, right, Operation::mult );
}

const Function operator/( const Parameter& left, const Function& right )
{
  return Function( left, right, Operation::div );
}



const Function operator+( const ParameterExpr& left, const Function& right )
{
  return Function( left, right, Operation::plus );
}

const Function operator-( const ParameterExpr& left, const Function& right )
{
  return Function( left, right, Operation::minus );
}

const Function operator*( const ParameterExpr& left, const Function& right )
{
  return Function( left, right, Operation::mult );
}

const Function operator/( const ParameterExpr& left, const Function& right )
{
  return Function( left, right, Operation::div );
}

