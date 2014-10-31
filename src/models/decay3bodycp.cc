
#include <complex>

#include <cfit/dataset.hh>
#include <cfit/function.hh>
#include <cfit/random.hh>

#include <cfit/models/decay3bodycp.hh>


Decay3BodyCP::Decay3BodyCP( const Variable&   mSq12  ,
                            const Variable&   mSq13  ,
                            const Variable&   mSq23  ,
                            const Amplitude&  amp    ,
                            const CoefExpr&   z      ,
                            const PhaseSpace& ps     ,
                            bool              docache  )
  : _amp( amp ), _ps( ps ), _hasKappa( false ), _z( z ),
    _nDir( 0.0 ), _nCnj( 0.0 ), _nXed( 0.0 ), _norm( 1.0 ), _fixedAmp( false ),
    _maxPdf( 14.0 ), _cacheAmps( false ), _ampDirCache( 0 ), _ampCnjCache( 0 )
{
  const std::map< std::string, Parameter >& pars = _amp.getPars();
  typedef std::map< const std::string, Parameter >::const_iterator pIter;
  for ( pIter par = pars.begin(); par != pars.end(); ++par )
    push( par->second );

  const std::map< std::string, Parameter >& zpars = _z.getPars();
  for ( pIter par = zpars.begin(); par != zpars.end(); ++par )
    push( par->second );

  push( mSq12 );
  push( mSq13 );
  push( mSq23 );

  // Do calculations common to all values of variables
  //    (usually compute norm).
  if ( docache )
    cache();
}




Decay3BodyCP::Decay3BodyCP( const Variable&   mSq12  ,
                            const Variable&   mSq13  ,
                            const Variable&   mSq23  ,
                            const Amplitude&  amp    ,
                            const CoefExpr&   z      ,
                            const Parameter&  kappa  ,
                            const PhaseSpace& ps     ,
                            bool              docache  )
  : _amp( amp ), _ps( ps ), _hasKappa( true ), _kappa( kappa ), _z( z ),
    _nDir( 0.0 ), _nCnj( 0.0 ), _nXed( 0.0 ), _norm( 1.0 ), _fixedAmp( false ),
    _maxPdf( 14.0 ), _cacheAmps( false ), _ampDirCache( 0 ), _ampCnjCache( 0 )
{
  const std::map< std::string, Parameter >& pars = _amp.getPars();
  typedef std::map< const std::string, Parameter >::const_iterator pIter;
  for ( pIter par = pars.begin(); par != pars.end(); ++par )
    push( par->second );

  const std::map< std::string, Parameter >& zpars = _z.getPars();
  for ( pIter par = zpars.begin(); par != zpars.end(); ++par )
    push( par->second );

  push( kappa );

  push( mSq12 );
  push( mSq13 );
  push( mSq23 );

  // Do calculations common to all values of variables
  //    (usually compute norm).
  if ( docache )
    cache();
}



Decay3BodyCP* Decay3BodyCP::copy() const
{
  return new Decay3BodyCP( *this );
}


// Need to overwrite setters defined in PdfModel, since function variables may need to be set.
void Decay3BodyCP::setVars( const std::vector< double >& vars ) throw( PdfException )
{
  if ( _varMap.size() != vars.size() )
    throw PdfException( "Decay3BodyCP::setVars: Number of arguments passed does not match number of required arguments." );

  typedef std::map< std::string, Variable >::iterator vIter;
  int index = 0;
  for ( vIter var = _varMap.begin(); var != _varMap.end(); ++var )
    var->second.setValue( vars[ index++ ] );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setVars( _varMap );
}


// Need to overwrite setters defined in PdfModel, since function variables may need to be set.
void Decay3BodyCP::setVars( const std::map< std::string, Variable >& vars ) throw( PdfException )
{
  typedef std::map< const std::string, Variable >::iterator vIter;
  for ( vIter var = _varMap.begin(); var != _varMap.end(); ++var )
    var->second.setValue( vars.find( var->first )->second.value() );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setVars( _varMap );
}


void Decay3BodyCP::setVars( const std::map< std::string, double >& vars ) throw( PdfException )
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
void Decay3BodyCP::setPars( const std::vector< double >& pars ) throw( PdfException )
{
  if ( _parMap.size() != pars.size() )
    throw PdfException( "Decay3BodyCP::setPars: Number of arguments passed does not match number of required arguments." );

  typedef std::map< std::string, Parameter >::iterator pIter;
  int index = 0;
  for ( pIter par = _parMap.begin(); par != _parMap.end(); ++par )
    par->second.setValue( pars[ index++ ] );

  _amp.setPars( _parMap );
  _z  .setPars( _parMap );

  if ( _hasKappa )
    _kappa.setValue( _parMap[ _kappa.name() ].value() );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setPars( _parMap );
}


// Need to overwrite setters defined in PdfModel, since function parameters may need to be set.
void Decay3BodyCP::setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException )
{
  typedef std::map< const std::string, Parameter >::iterator pIter;
  for ( pIter par = _parMap.begin(); par != _parMap.end(); ++par )
    par->second.setValue( pars.find( par->first )->second.value() );

  _amp.setPars( _parMap );
  _z  .setPars( _parMap );

  if ( _hasKappa )
    _kappa.setValue( pars.find( _kappa.name() )->second.value() );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setPars( _parMap );
}


// Need to overwrite setters defined in PdfModel, since function parameters may need to be set.
void Decay3BodyCP::setPars( const FunctionMinimum& min ) throw( PdfException )
{
  const std::vector< MinuitParameter >& pars = min.userParameters().parameters();
  typedef std::map< const std::string, Parameter >::iterator pIter;
  std::vector< MinuitParameter >::const_iterator mpar;
  for ( pIter par = _parMap.begin(); par != _parMap.end(); ++par )
  {
    mpar = std::find_if( pars.begin(), pars.end(),
                         [par]( const MinuitParameter& mpar )
                         { return mpar.name() == par->first; } );

    if ( mpar != pars.end() )
      par->second.setValue( mpar->value() );
  }

  _amp.setPars( _parMap );
  _z  .setPars( _parMap );

  if ( _hasKappa )
    _kappa.setValue( _parMap.find( _kappa.name() )->second.value() );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setPars( _parMap );
}



const double Decay3BodyCP::evaluateFuncs( const double& mSq12, const double& mSq13, const double& mSq23 ) const
{
  double value = 1.0;

  const std::string& name12 = getVar( 0 ).name(); // mSq12
  const std::string& name13 = getVar( 1 ).name(); // mSq13
  const std::string& name23 = getVar( 2 ).name(); // mSq23

  typedef std::vector< Function >::const_iterator fIter;

  std::map< std::string, double > varMap;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
  {
    if ( func->dependsOn( name12 ) ) varMap[ name12 ] = mSq12;
    if ( func->dependsOn( name13 ) ) varMap[ name13 ] = mSq13;
    if ( func->dependsOn( name23 ) ) varMap[ name23 ] = mSq23;

    value *= func->evaluate( varMap );
  }

  // Always return a non-negative value. Default to zero.
  return std::max( value, 0.0 );
}


const double Decay3BodyCP::evaluateFuncs( const double& mSq12, const double& mSq13 ) const
{
  const double& mSq23 = _ps.mSqSum() - mSq12 - mSq13;

  return evaluateFuncs( mSq12, mSq13, mSq23 );
};


// The values of the variables must be set with setVars before using this function.
const double Decay3BodyCP::evaluateFuncs() const
{
  double value = 1.0;

  typedef std::vector< Function >::const_iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    value *= func->evaluate();

  // Always return a non-negative value. Default to zero.
  return std::max( value, 0.0 );
}


void Decay3BodyCP::cache()
{
  const std::complex< double >& vz     = _z.evaluate();
  const double&                 vKappa = _hasKappa ? _kappa.value() : 1.0;

  if ( _fixedAmp )
  {
    _norm = _nDir + std::norm( vz ) * _nCnj + 2.0 * vKappa * real( vz * _nXed );
    return;
  }

  // Compute the value of _norm.
  _nDir = 0.0;
  _nCnj = 0.0;
  _nXed = 0.0;
  _norm = 0.0;

  // Define the properties of the integration method.
  const int    nBins = 400;
  const double min   = _ps.mSq12min();
  const double max   = _ps.mSq12max();
  const double step  = ( max - min ) / double( nBins );

  const double mSqSum = _ps.mSqSum();

  // Define the variables at each bin.
  double mSq12;
  double mSq13;
  double mSq23;

  double funcs;
  std::complex< double > ampDir;
  std::complex< double > ampCnj;

  // Compute the integral on the grid.
  for ( int binX = 0; binX < nBins; ++binX )
    for ( int binY = 0; binY < nBins; ++binY )
    {
      mSq12 = min + step * ( binX + 0.5 );
      mSq13 = min + step * ( binY + 0.5 );
      mSq23 = mSqSum - mSq12 - mSq13;

      // Proceed only if the point lies inside the kinematically allowed Dalitz region.
      // std::norm returns the squared modulus of the complex number, not its norm.
      if ( _ps.contains( mSq12, mSq13, mSq23 ) )
      {
        funcs = evaluateFuncs( mSq12, mSq13, mSq23 );
        ampDir = _amp.evaluate( _ps, mSq12, mSq13, mSq23 );
        ampCnj = _amp.evaluate( _ps, mSq13, mSq12, mSq23 );

        _nDir += std::norm( ampDir ) * funcs;
        _nCnj += std::norm( ampCnj ) * funcs;
        _nXed += conj( ampDir ) * ampCnj * funcs;
      }
    }

  const double& stepSq = std::pow( step, 2 );
  _nDir *= stepSq;
  _nCnj *= stepSq;
  _nXed *= stepSq;

  _fixedAmp = _amp.isFixed();

  _norm  = _nDir + std::norm( vz ) * _nCnj + 2.0 * vKappa * std::real( vz * _nXed );

  return;
}


const std::map< unsigned, std::vector< std::complex< double > > > Decay3BodyCP::cacheComplex( const Dataset& data )
{
  // Determine whether the amplitudes should be cached, i.e. only if all their parameters are fixed.
  _cacheAmps = true;

  std::map< unsigned, std::vector< std::complex< double > > > cached;

  if ( ! _cacheAmps )
    return cached;

  // Get an index for the cached complex amplitudes.
  _ampDirCache = assignCacheIdxComplex();
  _ampCnjCache = assignCacheIdxComplex();

  const std::string& mSq12name = getVar( 0 ).name();
  const std::string& mSq13name = getVar( 1 ).name();
  const std::string& mSq23name = getVar( 2 ).name();

  double mSq12;
  double mSq13;
  double mSq23;

  const std::size_t& size = data.size();
  for ( std::size_t entry = 0; entry < size; ++entry )
  {
    mSq12 = data.value( mSq12name, entry );
    mSq13 = data.value( mSq13name, entry );
    mSq23 = data.value( mSq23name, entry );

    cached[ _ampDirCache ].push_back( _amp.evaluate( _ps, mSq12, mSq13, mSq23 ) );
    cached[ _ampCnjCache ].push_back( _amp.evaluate( _ps, mSq13, mSq12, mSq23 ) );
  }

  return cached;
}


// Unnormalized evaluation.
const double Decay3BodyCP::evaluateUnnorm( const double& mSq12, const double& mSq13, const double& mSq23 ) const throw( PdfException )
{
  // Phase space amplitude of the decay of the particle.
  std::complex< double > ampDir = _amp.evaluate( _ps, mSq12, mSq13, mSq23 );
  std::complex< double > ampCnj = _amp.evaluate( _ps, mSq13, mSq12, mSq23 );

  const std::complex< double >& vz = _z.evaluate();

  if ( ! _hasKappa )
    // std::norm returns the squared modulus of the complex number, not its norm.
    return std::norm( ampDir + vz * ampCnj ); // * evaluateFuncs( mSq12, mSq13, mSq23 );

  const std::complex< double >& interf = conj( ampDir ) * ampCnj;

  double ampSq = std::norm( ampDir ) + std::norm( vz ) * std::norm( ampCnj );
  ampSq += 2.0 * _kappa.value() * std::real( vz * interf );

  return ampSq; // * evaluateFuncs( mSq12, mSq13, mSq23 );
}


const double Decay3BodyCP::evaluate( const double& mSq12, const double& mSq13, const double& mSq23 ) const throw( PdfException )
{
  return evaluateUnnorm( mSq12, mSq13, mSq23 ) / _norm;
}


const double Decay3BodyCP::evaluate( const double& mSq12, const double& mSq13 ) const throw( PdfException )
{
  const double& mSq23 = _ps.mSqSum() - mSq12 - mSq13;

  return evaluateUnnorm( mSq12, mSq13, mSq23 ) / _norm;
}


const double Decay3BodyCP::evaluate() const throw( PdfException )
{
  throw PdfException( "Decay3BodyCP: do not use evaluate() without arguments." );
}


const double Decay3BodyCP::project( const std::string& varName, const double& x ) const throw( PdfException )
{
  // Find the index of the variable to be projected.
  int index = -1;
  for ( unsigned var = 0; var < 3; ++var )
    if ( varName == getVar( var ).name() )
      index = var;

  // If the pdf does not depend on the passed variable name, the projection is 1.
  if ( index == -1 )
    return 1.0;

  // Minimum and maximum values of the variable wrt which the integration
  //    is done i.e. the next to the projected variable.
  const double& min = _ps.mSqMin( ( index + 1 ) % 3 );
  const double& max = _ps.mSqMax( ( index + 1 ) % 3 );

  // Integrate the model.
  const int& nbins = 400;
  double proj = 0.0;
  for ( int yBin = 0; yBin < nbins; ++yBin )
  {
    double y = binCenter( yBin, nbins, min, max );
    double z = _ps.mSqSum() - x - y;

    if ( index == 0 ) proj += evaluate( x, y, z );
    if ( index == 1 ) proj += evaluate( z, x, y );
    if ( index == 2 ) proj += evaluate( y, z, x );
  }

  proj *= ( max - min ) / double( nbins );

  return proj;
}


const double Decay3BodyCP::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  const std::size_t& size = vars.size();

  if ( size == 2 )
    return evaluate( vars[ 0 ], vars[ 1 ] );

  if ( size == 3 )
    return evaluate( vars[ 0 ], vars[ 1 ], vars[ 2 ] );

  throw PdfException( "Decay3BodyCP can only take either 2 or 3 arguments." );
}



const double Decay3BodyCP::evaluate( const std::vector< double >&                 vars  ,
                                     const std::vector< double >&                 cacheR,
                                     const std::vector< std::complex< double > >& cacheC ) const throw( PdfException )
{
  if ( ! _cacheAmps )
    return evaluate( vars );

  const std::size_t& size = vars.size();

  std::complex< double > ampDir = cacheC[ _ampDirCache ];
  std::complex< double > ampCnj = cacheC[ _ampCnjCache ];

  const std::complex< double >& vz = _z.evaluate();

  double funcs;
  if ( size == 2 )
    funcs = evaluateFuncs( vars[ 0 ], vars[ 1 ] );
  else if ( size == 3 )
    funcs = evaluateFuncs( vars[ 0 ], vars[ 1 ], vars[ 2 ] );
  else
    throw PdfException( "Decay3BodyCP can only take either 2 or 3 arguments." );

  if ( ! _hasKappa )
    return std::norm( ampDir + vz * ampCnj ) * funcs / _norm;

  const std::complex< double >& interf = std::conj( ampDir ) * ampCnj;

  double ampSq = std::norm( ampDir ) + std::norm( vz ) * std::norm( ampCnj );
  ampSq += 2.0 * _kappa.value() * std::real( vz * interf );

  return ampSq * funcs / _norm;
}




// No need to append an operator, since it can only be multiplication.
const Decay3BodyCP& Decay3BodyCP::operator*=( const Function& right ) throw( PdfException )
{
  // Check that the function does not depend on any variables that the model does not.
  const std::map< std::string, Variable >& varMap = right.getVarMap();
  for ( std::map< std::string, Variable >::const_iterator var = varMap.begin(); var != varMap.end(); ++var )
    if ( ! _varMap.count( var->second.name() ) )
      throw PdfException( "Cannot multiply a Decay3BodyCP pdf model by a function that depends on other variables." );

  // Consider the function parameters as own ones.
  const std::map< std::string, Parameter >& parMap = right.getParMap();
  _parMap.insert( parMap.begin(), parMap.end() );

  // Append the function to the functions vector.
  _funcs.push_back( right );

  // Recompute the norm, since the pdf shape has changed under this operation.
  cache();

  return *this;
}



const Decay3BodyCP operator*( Decay3BodyCP left, const Function& right )
{
  // Check that the function does not depend on any variables that the model does not.
  const std::map< std::string, Variable >& varMap = right.getVarMap();
  for ( std::map< std::string, Variable >::const_iterator var = varMap.begin(); var != varMap.end(); ++var )
    if ( ! left._varMap.count( var->second.name() ) )
      throw PdfException( "Cannot multiply a Decay3BodyCP pdf model by a function that depends on other variables." );

  // Consider the function parameters as own ones.
  const std::map< std::string, Parameter >& parMap = right.getParMap();
  left._parMap.insert( parMap.begin(), parMap.end() );

  // Append the function to the functions vector.
  left._funcs.push_back( right );

  // Recompute the norm, since the pdf shape has changed under this operation.
  left.cache();

  return left;
}



const Decay3BodyCP operator*( const Function& left, Decay3BodyCP right )
{
  // Check that the function does not depend on any variables that the model does not.
  const std::map< std::string, Variable >& varMap = left.getVarMap();
  for ( std::map< std::string, Variable >::const_iterator var = varMap.begin(); var != varMap.end(); ++var )
    if ( ! right._varMap.count( var->second.name() ) )
      throw PdfException( "Cannot multiply a Decay3BodyCP pdf model by a function that depends on other variables." );

  // Consider the function parameters as own ones.
  const std::map< std::string, Parameter >& parMap = left.getParMap();
  right._parMap.insert( parMap.begin(), parMap.end() );

  // Append the function to the functions vector.
  right._funcs.push_back( left );

  // Recompute the norm, since the pdf shape has changed under this operation.
  right.cache();

  return right;
}


const std::map< std::string, double > Decay3BodyCP::generate() const throw( PdfException )
{
  // Generate mSq12 and mSq13, and compute mSq23 from these.
  const double& min12 = std::pow( _ps.m1()      + _ps.m2(), 2 );
  const double& min13 = std::pow( _ps.m1()      + _ps.m3(), 2 );
  const double& max12 = std::pow( _ps.mMother() - _ps.m3(), 2 );
  const double& max13 = std::pow( _ps.mMother() - _ps.m2(), 2 );

  const std::string& mSq12name = getVar( 0 ).name();
  const std::string& mSq13name = getVar( 1 ).name();
  const std::string& mSq23name = getVar( 2 ).name();

  // Sum of squared invariant masses of all particles (mother and daughters).
  const double& mSqSum = _ps.mSqMother() + _ps.mSq1() + _ps.mSq2() + _ps.mSq3();

  double mSq12 = 0.0;
  double mSq13 = 0.0;
  double mSq23 = 0.0;

  double pdfVal  = 0.0;

  // Attempts to generate an event.
  int count = 10000;

  std::map< std::string, double > values;

  while ( count-- )
  {
    // Generate uniform mSq12 and mSq13.
    mSq12 = Random::flat( min12, max12 );
    mSq13 = Random::flat( min13, max13 );
    mSq23 = mSqSum - mSq12 - mSq13;

    values[ mSq12name ] = mSq12;
    values[ mSq13name ] = mSq13;
    values[ mSq23name ] = mSq23;

    pdfVal = this->evaluate( mSq12, mSq13, mSq23 );

    if ( pdfVal > _maxPdf )
      std::cout << "Problem: " << pdfVal << " > " << _maxPdf
                << " for variables " << mSq12 << " " << mSq13 << " " << mSq23 << std::endl;

    // Apply the accept-reject decision.
    if ( Random::flat( 0.0, _maxPdf ) < pdfVal )
      return values;
  }

  values[ mSq12name ] = 0.0;
  values[ mSq13name ] = 0.0;
  values[ mSq23name ] = 0.0;

  return values;
}
