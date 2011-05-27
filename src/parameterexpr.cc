
// #include <cfit/parameter.hh>
#include <cfit/parameterexpr.hh>


void ParameterExpr::append( const double& val )
{
  _consts.push_back( val );
  _expression += "c"; // c = constant.
}

void ParameterExpr::append( const Parameter& par )
{
  _pars.push_back( par );
  _expression += "p"; // p = parameter.
}

void ParameterExpr::append( const ParameterExpr& expr )
{
  _consts.insert( _consts.end(), expr._consts.begin(), expr._consts.end() );
  _pars  .insert( _pars  .end(), expr._pars  .begin(), expr._pars  .end() );
  _opers .insert( _opers .end(), expr._opers .begin(), expr._opers .end() );
  _expression += expr._expression;
}

