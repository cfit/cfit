
#include <iostream> // TEMPORAL
#include <iterator> // TEMPORAL

#include <sstream>
#include <vector>
#include <stack>
#include <cmath>

#include <cfit/parameter.hh>
#include <cfit/parameterexpr.hh>
#include <cfit/pdfmodel.hh>
#include <cfit/pdf.hh>

// Append a model.
void Pdf::append( const PdfModel& model )
{
  _vars   .insert( model._vars   .begin(), model._vars   .end() );
  _pars   .insert( model._pars   .begin(), model._pars   .end() );
  _pdfVect.push_back( const_cast<PdfModel*>( &model ) );

  _expression += "m"; // m = model.
}

// Append a pdf expression.
void Pdf::append( const Pdf& pdf )
{
  _vars   .insert(                 pdf._vars   .begin(), pdf._vars   .end() );
  _pars   .insert(                 pdf._pars   .begin(), pdf._pars   .end() );
  _opers  .insert( _opers  .end(), pdf._opers  .begin(), pdf._opers  .end() );
  _consts .insert( _consts .end(), pdf._consts .begin(), pdf._consts .end() );
  _parVect.insert( _parVect.end(), pdf._parVect.begin(), pdf._parVect.end() );
  _pdfVect.insert( _pdfVect.end(), pdf._pdfVect.begin(), pdf._pdfVect.end() );

  _expression += pdf._expression;
}

// Append a parameter.
void Pdf::append( const Parameter& par )
{
  _pars[ par.name() ] = par;
  _parVect.push_back( par.name() );

  _expression += "p"; // p = parameter.
}

// Append a parameter expression.
void Pdf::append( const ParameterExpr& expr )
{
  for ( std::vector< Parameter >::const_iterator par = expr._pars.begin(); par != expr._pars.end(); ++par )
    _pars[ par->name() ] = *par;
  _opers  .insert( _opers  .end(), expr._opers .begin(), expr._opers .end() );
  _consts .insert( _consts .end(), expr._consts.begin(), expr._consts.end() );
  _parVect.insert( _parVect.end(), expr._pars  .begin(), expr._pars  .end() );

  _expression += expr._expression;
}

// Append a constant.
void Pdf::append( const double& ctnt )
{
  _consts.push_back( ctnt );

  _expression += "c"; // c = constant.
}

// Append a binary operation. No unary operation should ever be appended.
void Pdf::append( const Operation::Op& oper )
{
  _expression += "b";
  _opers.push_back( oper );
}


// Assignment operation.
const Pdf& Pdf::operator=( const PdfModel& right )
{
  append( right );
  return *this;
}

// Assignment operations with a pdf model.
const Pdf& Pdf::operator+=( const PdfModel& right ) throw( PdfException )
{
  // Cannot add two pdfs that do not depend on exactly the same variables.
  if ( this->varNames() != right.varNames() )
    throw PdfException( "Cannot add two pdfs that do not depend on the same variables." );

  append( right           );
  append( Operation::plus );
  return *this;
}

const Pdf& Pdf::operator*=( const PdfModel& right ) throw( PdfException )
{
  // Cannot multiply two pdfs that share some variable.
  std::vector< std::string > theseVars = this->varNames();
  std::vector< std::string > rightVars = right.varNames();
  std::vector< std::string > intersect;
  std::set_intersection( theseVars.begin(), theseVars.end(),
			 rightVars.begin(), rightVars.end(),
			 std::back_inserter( intersect ) );

  if ( ! intersect.empty() )
    throw PdfException( "Cannot multiply two pdfs that depend on some common variable." );

  append( right           );
  append( Operation::mult );
  return *this;
}

// Assignment operations with a pdf expression.
const Pdf& Pdf::operator+=( const Pdf& right ) throw( PdfException )
{
  // Cannot add two pdfs that do not depend on exactly the same variables.
  if ( this->varNames() != right.varNames() )
    throw PdfException( "Cannot add two pdfs that do not depend on the same variables." );

  append( right           );
  append( Operation::plus );
  return *this;
}

const Pdf& Pdf::operator*=( const Pdf& right ) throw( PdfException )
{
  // Cannot multiply two pdfs that share some variable.
  std::vector< std::string > theseVars = this->varNames();
  std::vector< std::string > rightVars = right.varNames();
  std::vector< std::string > intersect;
  std::set_intersection( theseVars.begin(), theseVars.end(),
			 rightVars.begin(), rightVars.end(),
			 std::back_inserter( intersect ) );

  if ( ! intersect.empty() )
    throw PdfException( "Cannot multiply two pdfs that depend on some common variable." );

  append( right           );
  append( Operation::mult );
  return *this;
}

// Assignment operations with a parameter.
const Pdf& Pdf::operator*=( const Parameter& right )
{
  append( right           );
  append( Operation::mult );
  return *this;
}

const Pdf& Pdf::operator/=( const Parameter& right )
{
  append( right          );
  append( Operation::div );
  return *this;
}

// Assignment operations with a parameter expression.
const Pdf& Pdf::operator*=( const ParameterExpr& right )
{
  append( right           );
  append( Operation::mult );
  return *this;
}

const Pdf& Pdf::operator/=( const ParameterExpr& right )
{
  append( right          );
  append( Operation::div );
  return *this;
}

// Assignment operations with a constant.
const Pdf& Pdf::operator*=( const double& right )
{
  append( right           );
  append( Operation::mult );
  return *this;
}

const Pdf& Pdf::operator/=( const double& right )
{
  append( right          );
  append( Operation::div );
  return *this;
}



void Pdf::setVars( const std::vector< double >& vars ) throw( PdfException )
{
  if ( _vars.size() != vars.size() )
    throw PdfException( "Number of arguments passed does not match number of required arguments." );

  // Set the local values of the variables.
  typedef std::map< std::string, Variable >::iterator vIter;
  int index = 0;
  for ( vIter var = _vars.begin(); var != _vars.end(); ++var )
    var->second.setValue( vars[ index++ ] );

  // Propagate the values to the list of pdfs.
  typedef std::vector< PdfModel* >::const_iterator pdfIter;
  for ( pdfIter pdf = _pdfVect.begin(); pdf != _pdfVect.end(); ++pdf )
    {
      std::map< std::string, Variable >& pdfVars = (*pdf)->_vars;
      for ( vIter var = pdfVars.begin(); var != pdfVars.end(); ++var )
	var->second.setValue( _vars[ var->second.name() ].value() );
    }
}


void Pdf::setPars( const std::vector< double >& pars ) throw( PdfException )
{
  if ( _pars.size() != pars.size() )
    throw PdfException( "Number of arguments passed does not match number of required arguments." );

  // Set the local values of the parameters.
  typedef std::map< std::string, Parameter >::iterator pIter;
  int index = 0;
  for ( pIter par = _pars.begin(); par != _pars.end(); par++ )
    par->second.setValue( pars[ index++ ] );

  // Propagate the values to the list of pdfs.
  typedef std::vector< PdfModel* >::const_iterator pdfIter;
  for ( pdfIter pdf = _pdfVect.begin(); pdf != _pdfVect.end(); ++pdf )
    {
      std::map< std::string, Parameter >& pdfPars = (*pdf)->_pars;
      for ( pIter par = pdfPars.begin(); par != pdfPars.end(); ++par )
	par->second.setValue( _pars[ par->second.name() ].value() );
    }
}


void Pdf::cache()
{
  typedef std::vector< PdfModel* >::const_iterator pIter;
  for ( pIter pdf = _pdfVect.begin(); pdf != _pdfVect.end(); pdf++ )
    (*pdf)->cache();

  return;
}

// Binary operations.
double Pdf::operate( const double& x, const double& y, const Operation::Op& oper ) throw( PdfException )
{
  if ( oper == Operation::plus )
    return x + y;
  if ( oper == Operation::minus )
    return x - y;
  if ( oper == Operation::mult )
    return x * y;
  if ( oper == Operation::div )
    return x / y;
  if ( oper == Operation::pow )
    return std::pow( x, y );

  throw PdfException( std::string( "Parse error: unknown binary operation " ) + Operation::tostring( oper ) + "." );
}

// Unary operations.
double Pdf::operate( const double& x, const Operation::Op& oper ) throw( PdfException )
{
  if ( oper == Operation::minus )
    return -x;
  if ( oper == Operation::exp )
    return std::exp( x );
  if ( oper == Operation::log )
    return std::log( x );
  if ( oper == Operation::sin )
    return std::sin( x );
  if ( oper == Operation::cos )
    return std::cos( x );
  if ( oper == Operation::tan )
    return std::tan( x );

  throw PdfException( std::string( "Parse error: unknown unary operation " ) + Operation::tostring( oper ) + "." );
}

// Before running this function, the Pdf::setVars( vars ) function must be called.
//    To avoid the risk of forgetting it, run Pdf::evaluate( vars ).
double Pdf::evaluate() const throw( PdfException )
{
  std::stack< double > values;

  double x;
  double y;
  std::vector< PdfModel*     >::const_iterator pdf = _pdfVect.begin();
  std::vector< Parameter     >::const_iterator par = _parVect.begin();
  std::vector< double        >::const_iterator ctt = _consts .begin();
  std::vector< Operation::Op >::const_iterator ops = _opers  .begin();

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ch++ )
    if ( *ch == 'm' )
      values.push( (*pdf++)->evaluate() );
    else if ( *ch == 'p' )
      values.push( _pars.find( par++->name() )->second.value() );
    else if ( *ch == 'c' )
      values.push( *ctt++ );
    else
      if ( values.size() < 2 )
	throw PdfException( "Parse error: not enough values in the stack." );
      else
	{
	  if ( *ch == 'b' )
	    {
	      y = values.top();
	      values.pop();
	      x = values.top();
	      values.pop();
	      values.push( operate( x, y, *ops++ ) );
	    }
	  else if ( *ch == 'u' )
	    {
	      x = values.top();
	      values.pop();
	      values.push( operate( x, *ops++ ) );
	    }
	  else
	    throw PdfException( std::string( "Parse error: unknown command " ) + *ch + "." );
	}

  if ( values.size() != 1 )
    throw PdfException( "Parse error: too many values have been supplied." );

  return values.top();
}


double Pdf::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  if ( _vars.size() != vars.size() )
    throw PdfException( "Number of arguments passed does not match number of required arguments." );

  // Dictionary of the variable names with the values passed.
  std::map< std::string, double > localVars;
  std::vector< double > modelVars;

  // Set the local values of the variables.
  int index = 0;
  typedef std::map< std::string, Variable >::const_iterator vIter;
  for ( vIter var = _vars.begin(); var != _vars.end(); var++ )
    localVars[ var->first ] = vars[ index++ ];

  std::stack< double > values;

  double x;
  double y;
  std::vector< PdfModel*     >::const_iterator pdf = _pdfVect.begin();
  std::vector< Parameter     >::const_iterator par = _parVect.begin();
  std::vector< double        >::const_iterator ctt = _consts .begin();
  std::vector< Operation::Op >::const_iterator ops = _opers  .begin();

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
    if ( *ch == 'm' )
      {
	// Determine the variables that the pdf depends on.
	modelVars.clear();
	std::map< std::string, Variable >& pdfVars = (*pdf)->_vars;
	for ( vIter var = pdfVars.begin(); var != pdfVars.end(); ++var )
	  modelVars.push_back( localVars.find( var->second.name() )->second );

	// Evaluate the function at the given point.
	values.push( (*pdf++)->evaluate( modelVars ) );
      }
    else if ( *ch == 'p' )
      values.push( _pars.find( par++->name() )->second.value() );
    else if ( *ch == 'c' )
      values.push( *ctt++ );
    else
      if ( values.size() < 2 )
	throw PdfException( "Parse error: not enough values in the stack." );
      else
	{
	  if ( *ch == 'b' )
	    {
	      y = values.top();
	      values.pop();
	      x = values.top();
	      values.pop();
	      values.push( operate( x, y, *ops++ ) );
	    }
	  else if ( *ch == 'u' )
	    {
	      x = values.top();
	      values.pop();
	      values.push( operate( x, *ops++ ) );
	    }
	  else
	    throw PdfException( std::string( "Parse error: unknown operation " ) + *ch + "." );
	}

  if ( values.size() != 1 )
    throw PdfException( "Parse error: too many values have been supplied." );

  return values.top();
}



// Get a vector with the names of the variables common in all the
//    products of pdfs. This is the list of variables that can be
//    integrated over in a convolution operation.
std::vector< std::string > Pdf::commonVars() const throw( PdfException )
{
  // Vectors that will contain the variables each model depends on.
  std::vector< std::string > x;
  std::vector< std::string > y;
  std::vector< std::string > z;

  // Empty vector.
  std::vector< std::string > voidVec;

  // Stack for partial calculations.
  std::stack< std::vector< std::string > > calcs;

  // Iterator over the models that the pdf depends on.
  std::vector< PdfModel* >::const_iterator models = _pdfVect.begin();

  typedef std::string::const_iterator eIter;
  std::vector< std::string > varNames;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
    if ( *ch == 'm' )
      {
	std::map< std::string, Variable >& vars = (*models++)->_vars;
	std::transform( vars.begin(), vars.end(), std::back_inserter( varNames ), Select1st() );
	calcs.push( varNames );
      }
    else if ( *ch == 'p' )
      calcs.push( voidVec );
    else
      if ( calcs.size() < 2 )
	throw PdfException( "Parse error computing convolution: not enough values in the stack." );
      else
	{
	  x = calcs.top();
	  calcs.pop();
	  y = calcs.top();
	  calcs.pop();
	  z.clear();

	  if ( *ch == '+' )
	    std::set_intersection( x.begin(), x.end(), y.begin(), y.end(), std::back_inserter( z ) );
	  else if ( *ch == '*' )
	    std::set_union       ( x.begin(), x.end(), y.begin(), y.end(), std::back_inserter( z ) );

	  calcs.push( z );
	}

  if ( calcs.size() != 1 )
    throw PdfException( "Parse error computing convolution: too many values have been supplied." );

  return calcs.top();
}

