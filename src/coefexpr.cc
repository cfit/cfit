
#include <vector>
#include <stack>

// #include <cfit/coef.hh>
#include <cfit/coefexpr.hh>


void CoefExpr::append( const std::complex< double >& val )
{
  _ctnts.push_back( val );
  _expression += "c"; // c = constant.
}

void CoefExpr::append( const double& val )
{
  _ctnts.push_back( val );
  _expression += "c"; // c = constant.
}

void CoefExpr::append( const Parameter& par )
{
  _parms.push_back( par );
  _expression += "p"; // p = parameter.
}

void CoefExpr::append( const ParameterExpr& expr )
{
  // Inserting real constants as complex should not be a problem.
  _ctnts.insert( _ctnts.end(), expr._ctnts.begin(), expr._ctnts.end() );
  _parms.insert( _parms.end(), expr._parms.begin(), expr._parms.end() );
  _opers.insert( _opers.end(), expr._opers.begin(), expr._opers.end() );

  _expression += expr._expression; // p = parameter.
}

void CoefExpr::append( const Coef& coef )
{
  _coefs.push_back( coef );
  _expression += "k"; // k = coefficient.
}

void CoefExpr::append( const CoefExpr& expr )
{
  _ctnts.insert( _ctnts.end(), expr._ctnts.begin(), expr._ctnts.end() );
  _parms.insert( _parms.end(), expr._parms.begin(), expr._parms.end() );
  _coefs.insert( _coefs.end(), expr._coefs.begin(), expr._coefs.end() );
  _opers.insert( _opers.end(), expr._opers.begin(), expr._opers.end() );

  _expression += expr._expression;
}



void CoefExpr::setPars( const std::map< std::string, Parameter >& pars )
{
  // Set the values of the parameters.
  typedef std::vector< Parameter >::iterator pIter;
  for ( pIter par = _parms.begin(); par != _parms.end(); ++par )
    par->setValue( pars.find( par->name() )->second.value() );

  // Set the values of the coefficients.
  typedef std::vector< Coef >::iterator cIter;
  for ( cIter coef = _coefs.begin(); coef != _coefs.end(); ++coef )
    coef->setValue( pars.find( coef->real().name() )->second.value(),
                    pars.find( coef->imag().name() )->second.value() );
}


void CoefExpr::setPars( const std::map< std::string, double >& pars )
{
  // Set the values of the parameters.
  typedef std::vector< Parameter >::iterator pIter;
  for ( pIter par = _parms.begin(); par != _parms.end(); ++par )
    par->setValue( pars.find( par->name() )->second );

  // Set the values of the coefficients.
  typedef std::vector< Coef >::iterator cIter;
  for ( cIter coef = _coefs.begin(); coef != _coefs.end(); ++coef )
    coef->setValue( pars.find( coef->real().name() )->second,
                    pars.find( coef->imag().name() )->second );
}


const std::map< std::string, Parameter > CoefExpr::getPars() const
{
  std::map< std::string, Parameter > parMap;

  typedef std::vector< Parameter >::const_iterator pIter;
  for ( pIter par = _parms.begin(); par != _parms.end(); ++par )
    parMap.emplace( par->name(), *par );

  return parMap;
}



const std::complex< double > CoefExpr::evaluate() const throw( PdfException )
{
  std::stack< std::complex< double > > values;

  std::complex< double > x;
  std::complex< double > y;
  std::vector< Coef                   >::const_iterator coe = _coefs.begin();
  std::vector< Parameter              >::const_iterator par = _parms.begin();
  std::vector< std::complex< double > >::const_iterator ctt = _ctnts.begin();
  std::vector< Operation::Op >::const_iterator ops = _opers.begin();

  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
  {
    if ( *ch == 'k' )
      values.push( (coe++)->value() );
    else if ( *ch == 'p' )
      values.push( (par++)->value() );
    else if ( *ch == 'c' )
      values.push( *ctt++ );
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
