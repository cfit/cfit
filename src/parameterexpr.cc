
#include <stack>

#include <cfit/parameterexpr.hh>


void ParameterExpr::append( const double& val )
{
  _ctnts.push_back( val );
  _expression += "c"; // c = constant.
}

void ParameterExpr::append( const Parameter& par )
{
  _parms.push_back( par );
  _expression += "p"; // p = parameter.
}

void ParameterExpr::append( const ParameterExpr& expr )
{
  _ctnts.insert( _ctnts.end(), expr._ctnts.begin(), expr._ctnts.end() );
  _parms.insert( _parms.end(), expr._parms.begin(), expr._parms.end() );
  _opers.insert( _opers.end(), expr._opers.begin(), expr._opers.end() );
  _expression += expr._expression;
}

void ParameterExpr::append( const Operation::Op& oper )
{
  _opers.push_back( oper );
  _expression += "b"; // b = binary operation.
}


// Return a map from parameter names to parameters for all parameters used in the expression.
const std::map< std::string, Parameter > ParameterExpr::getPars() const
{
  std::map< std::string, Parameter > parMap;

  typedef std::vector< Parameter >::const_iterator pIter;
  for ( pIter par = _parms.begin(); par != _parms.end(); ++par )
    parMap.emplace( par->name(), *par );

  return parMap;
}


void ParameterExpr::setPars( const std::map< std::string, Parameter >& pars )
{
  // Set the values of the parameters.
  typedef std::vector< Parameter >::iterator pIter;
  for ( pIter par = _parms.begin(); par != _parms.end(); ++par )
    par->setValue( pars.find( par->name() )->second.value() );
}


void ParameterExpr::setPars( const std::map< std::string, double >& pars )
{
  // Set the values of the parameters.
  typedef std::vector< Parameter >::iterator pIter;
  for ( pIter par = _parms.begin(); par != _parms.end(); ++par )
    par->setValue( pars.find( par->name() )->second );
}


const double ParameterExpr::evaluate() const
{
  std::stack< double > values;

  double x;
  double y;
  std::vector< double        >::const_iterator ctt = _ctnts.begin();
  std::vector< Parameter     >::const_iterator par = _parms.begin();
  std::vector< Operation::Op >::const_iterator ops = _opers.begin();

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
  {
    if ( *ch == 'c' )
      values.push( *ctt++ );
    else if ( *ch == 'p' )
      values.push( (par++)->value() );
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
  }

  if ( values.size() != 1 )
    throw PdfException( "PdfExpr parse error: too many values have been supplied." );

  return values.top();
}



const ParameterExpr& ParameterExpr::operator=( const Parameter& right )
{
  clear();
  append( right );

  return *this;
}



const ParameterExpr& ParameterExpr::operator+=( const ParameterExpr& right )
{
  if ( _expression.empty() )
  {
    append( right );
    return *this;
  }

  append( right           );
  append( Operation::plus );

  return *this;
}


const ParameterExpr& ParameterExpr::operator-=( const ParameterExpr& right )
{
  if ( _expression.empty() )
  {
    append( right );
    _opers.push_back( Operation::minus );
    _expression += "u"; // u = unary operation.
    return *this;
  }

  append( right            );
  append( Operation::minus );

  return *this;
}


const ParameterExpr& ParameterExpr::operator*=( const ParameterExpr& right )
{
  if ( _expression.empty() )
    return *this;

  append( right           );
  append( Operation::mult );

  return *this;
}


const ParameterExpr& ParameterExpr::operator/=( const ParameterExpr& right )
{
  if ( _expression.empty() )
    return *this;

  append( right          );
  append( Operation::div );

  return *this;
}



const ParameterExpr& ParameterExpr::operator+=( const Parameter& right )
{
  if ( _expression.empty() )
  {
    append( right );
    return *this;
  }

  append( right           );
  append( Operation::plus );

  return *this;
}


const ParameterExpr& ParameterExpr::operator-=( const Parameter& right )
{
  if ( _expression.empty() )
  {
    append( right );
    _opers.push_back( Operation::minus );
    _expression += "u"; // u = unary operation.
    return *this;
  }

  append( right            );
  append( Operation::minus );

  return *this;
}


const ParameterExpr& ParameterExpr::operator*=( const Parameter& right )
{
  if ( _expression.empty() )
    return *this;

  append( right           );
  append( Operation::mult );

  return *this;
}


const ParameterExpr& ParameterExpr::operator/=( const Parameter& right )
{
  if ( _expression.empty() )
    return *this;

  append( right          );
  append( Operation::div );

  return *this;
}



const ParameterExpr& ParameterExpr::operator+=( const double& right )
{
  if ( _expression.empty() )
  {
    append( right );
    return *this;
  }

  append( right           );
  append( Operation::plus );

  return *this;
}


const ParameterExpr& ParameterExpr::operator-=( const double& right )
{
  if ( _expression.empty() )
  {
    append( right );
    _opers.push_back( Operation::minus );
    _expression += "u"; // u = unary operation.
    return *this;
  }

  append( right            );
  append( Operation::minus );

  return *this;
}


const ParameterExpr& ParameterExpr::operator*=( const double& right )
{
  if ( _expression.empty() )
    return *this;

  append( right           );
  append( Operation::mult );

  return *this;
}


const ParameterExpr& ParameterExpr::operator/=( const double& right )
{
  if ( _expression.empty() )
    return *this;

  append( right          );
  append( Operation::div );

  return *this;
}
