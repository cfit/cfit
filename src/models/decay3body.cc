
#include <cfit/models/decay3body.hh>
#include <cfit/function.hh>

Decay3Body::Decay3Body( const Variable&   mSq12,
			const Variable&   mSq13,
			const Variable&   mSq23,
			const Amplitude&  amp  ,
			const PhaseSpace& ps     )
  : DecayModel( mSq12, mSq13, mSq23, amp, ps ), _norm( 1. )
{
  // Do calculations common to all values of variables
  //    (usually compute norm).
  cache();
}


Decay3Body* Decay3Body::copy() const
{
  return new Decay3Body( *this );
}


// Need to overwrite setters defined in PdfModel, since function variables may need to be set.
void Decay3Body::setVars( const std::vector< double >& vars ) throw( PdfException )
{
  if ( _varMap.size() != vars.size() )
    throw PdfException( "Number of arguments passed does not match number of required arguments." );

  typedef std::map< std::string, Variable >::iterator vIter;
  int index = 0;
  for ( vIter var = _varMap.begin(); var != _varMap.end(); ++var )
    var->second.setValue( vars[ index++ ] );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setVars( _varMap );
}


// Need to overwrite setters defined in PdfModel, since function variables may need to be set.
void Decay3Body::setVars( const std::map< std::string, Variable >& vars ) throw( PdfException )
{
  typedef std::map< const std::string, Variable >::iterator vIter;
  for ( vIter var = _varMap.begin(); var != _varMap.end(); ++var )
    var->second.setValue( vars.find( var->first )->second.value() );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setVars( _varMap );
}


void Decay3Body::setVars( const std::map< std::string, double >& vars ) throw( PdfException )
{
  typedef std::map< const std::string, Variable >::iterator vIter;
  for ( vIter var = _varMap.begin(); var != _varMap.end(); ++var )
    var->second.setValue( vars.find( var->first )->second );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setVars( _varMap );
}




// Set the parameters to those given as argument.
// They must be sorted alphabetically, since it's how MnUserParameters are passed
//    in the minimize function of the minimizers. It must be so, because pushing
//    two parameters with the same name would create confusion otherwise.
void Decay3Body::setPars( const std::vector< double >& pars ) throw( PdfException )
{
  if ( _parMap.size() != pars.size() )
    throw PdfException( "Number of arguments passed does not match number of required arguments." );

  typedef std::map< std::string, Parameter >::iterator pIter;
  int index = 0;
  for ( pIter par = _parMap.begin(); par != _parMap.end(); ++par )
    par->second.setValue( pars[ index++ ] );

  _amp.setPars( _parMap );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setPars( _parMap );
}


// Need to overwrite setters defined in PdfModel, since function parameters may need to be set.
void Decay3Body::setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException )
{
  typedef std::map< const std::string, Parameter >::iterator pIter;
  for ( pIter par = _parMap.begin(); par != _parMap.end(); ++par )
    par->second.setValue( pars.find( par->first )->second.value() );

  _amp.setPars( _parMap );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setPars( _parMap );
}


// Need to overwrite setters defined in PdfModel, since function parameters may need to be set.
void Decay3Body::setPars( const FunctionMinimum& min ) throw( PdfException )
{
  const MnUserParameters& pars = min.userParameters();
  typedef std::vector< MinuitParameter >::const_iterator pIter;
  for( pIter par = pars.parameters().begin(); par != pars.parameters().end(); ++par )
    _parMap[ par->name() ].set( par->value(), par->error() );

  _amp.setPars( _parMap );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setPars( _parMap );
}



const double Decay3Body::evaluateFuncs( const double& mSq12, const double& mSq13, const double& mSq23 )
{
  double value = 1.0;

  const std::string& name0 = getVar( 0 ).name(); // mSq12
  const std::string& name1 = getVar( 1 ).name(); // mSq13
  const std::string& name2 = getVar( 2 ).name(); // mSq23

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
  {
    if ( func->dependsOn( name0 ) ) func->setVar( name0, mSq12 );
    if ( func->dependsOn( name1 ) ) func->setVar( name1, mSq13 );
    if ( func->dependsOn( name2 ) ) func->setVar( name2, mSq23 );

    value *= func->evaluate();
  }

  return value;
}


// The values of the variables must be set with setVars before using this function.
const double Decay3Body::evaluateFuncs() const
{
  double value = 1.0;

  typedef std::vector< Function >::const_iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    value *= func->evaluate();

  return value;
}



void Decay3Body::cache()
{
  // Compute the value of _norm.
  _norm = 0.0;

  // Define the properties of the integration method.
  const int    nBins = 400;
  const double min   = std::pow( _ps.m1()      + _ps.m2(), 2 );
  const double max   = std::pow( _ps.mMother() - _ps.m3(), 2 );
  const double step  = ( max - min ) / double( nBins );

  const double mSqSum = _ps.mSqMother() + _ps.mSq1() + _ps.mSq2() + _ps.mSq3();

  // Define the variables at each bin.
  double mSq12;
  double mSq13;
  double mSq23;

  // Compute the integral on the grid.
  for ( int binX = 0; binX < nBins; ++binX )
    for ( int binY = 0; binY < nBins; ++binY )
    {
      mSq12 = min + step * ( binX + .5 );
      mSq13 = min + step * ( binY + .5 );
      mSq23 = mSqSum - mSq12 - mSq13;

      // Proceed only if the point lies inside the kinematically allowed Dalitz region.
      // std::norm returns the squared modulus of the complex number, not its norm.
      if ( _ps.contains( mSq12, mSq13, mSq23 ) )
        _norm += std::norm( _amp.evaluate( _ps, mSq12, mSq13, mSq23 ) ) * evaluateFuncs( mSq12, mSq13, mSq23 );
    }

  _norm *= std::pow( step, 2 );

  return;
}


const double Decay3Body::evaluate() const throw( PdfException )
{
  // Phase space amplitude of the decay of the particle.
  std::complex< double > amp = _amp.evaluate( _ps, mSq12(), mSq13(), mSq23() );

  // std::norm returns the squared modulus of the complex number, not its norm.
  return std::norm( amp ) * evaluateFuncs() / _norm;
}


const double Decay3Body::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  throw PdfException( "Do not use evaluate( vars ) in decay3body" );
}


// No need to append an operator, since it can only be multiplication.
const Decay3Body& Decay3Body::operator*=( const Function& right ) throw( PdfException )
{
  // Check that the function does not depend on any variables that the model does not.
  const std::map< std::string, Variable >& varMap = right.getVarMap();
  for ( std::map< std::string, Variable >::const_iterator var = varMap.begin(); var != varMap.end(); ++var )
    if ( ! _varMap.count( var->second.name() ) )
      throw PdfException( "Cannot multiply a Decay3Body pdf model by a function that depends on other variables." );

  // Consider the function parameters as own ones.
  const std::map< std::string, Parameter >& parMap = right.getParMap();
  _parMap.insert( parMap.begin(), parMap.end() );

  // Append the function to the functions vector.
  _funcs.push_back( right );

  // Recompute the norm, since the pdf shape has changed under this operation.
  cache();

  return *this;
}



const Decay3Body operator*( Decay3Body left, const Function& right )
{
  // Check that the function does not depend on any variables that the model does not.
  const std::map< std::string, Variable >& varMap = right.getVarMap();
  for ( std::map< std::string, Variable >::const_iterator var = varMap.begin(); var != varMap.end(); ++var )
    if ( ! left._varMap.count( var->second.name() ) )
      throw PdfException( "Cannot multiply a Decay3Body pdf model by a function that depends on other variables." );

  // Consider the function parameters as own ones.
  const std::map< std::string, Parameter >& parMap = right.getParMap();
  left._parMap.insert( parMap.begin(), parMap.end() );

  // Append the function to the functions vector.
  left._funcs.push_back( right );

  // Recompute the norm, since the pdf shape has changed under this operation.
  left.cache();

  return left;
}



const Decay3Body operator*( const Function& left, Decay3Body right )
{
  // Check that the function does not depend on any variables that the model does not.
  const std::map< std::string, Variable >& varMap = left.getVarMap();
  for ( std::map< std::string, Variable >::const_iterator var = varMap.begin(); var != varMap.end(); ++var )
    if ( ! right._varMap.count( var->second.name() ) )
      throw PdfException( "Cannot multiply a Decay3Body pdf model by a function that depends on other variables." );

  // Consider the function parameters as own ones.
  const std::map< std::string, Parameter >& parMap = left.getParMap();
  right._parMap.insert( parMap.begin(), parMap.end() );

  // Append the function to the functions vector.
  right._funcs.push_back( left );

  // Recompute the norm, since the pdf shape has changed under this operation.
  right.cache();

  return right;
}

