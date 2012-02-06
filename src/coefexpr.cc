
// #include <cfit/coef.hh>
#include <cfit/coefexpr.hh>


void CoefExpr::append( const std::complex< double >& val )
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
//   _ctnts.insert( _ctnts.end(), expr._ctnts.begin(), expr._ctnts.end() );
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

