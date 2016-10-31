
#include <complex>

#include <cfit/dataset.hh>
#include <cfit/function.hh>
#include <cfit/random.hh>

#include <cfit/models/decay3bodymix.hh>

// Main constructor.
Decay3BodyMix::Decay3BodyMix( const Variable&      mSq12  ,
                              const Variable&      mSq13  ,
                              const Variable&      mSq23  ,
                              const Variable&      t      ,
                              const ParameterExpr& width  ,
                              const Amplitude&     amp    ,
                              const CoefExpr&      z      ,
                              const PhaseSpace&    ps     ,
                              bool                 docache  )
  : DecayModel( mSq12, mSq13, mSq23, amp, ps ),
    _mSq12( mSq12.name() ), _mSq13( mSq13.name() ), _mSq23( mSq23.name() ), _t( t.name() ),
    _width( width ),
    _z( z ),
    _qoverp( 1.0 ),
    _hasMixing( true  ),
    _hasCPV   ( false ),
    _nDir( 0.0 ), _nCnj( 0.0 ), _nXed( 0.0 ), _norm( 1.0 ), _fixedAmp( false ),
    _maxPdf( 54.0 ), _cacheAmps( false ), _ampDirCache( 0 ), _ampCnjCache( 0 )
{
  // Make the variables available to cfit.
  // The squared invariant masses are already made available by the DecayModel constructor.
  push( t     );

  // Make all the parameters in width and z available to cfit.
  push( width );
  push( z     );

  // Do calculations common to all values of variables (compute norm).
  if ( docache )
    cache();
}


// Main constructor.
Decay3BodyMix::Decay3BodyMix( const Variable&      mSq12  ,
                              const Variable&      mSq13  ,
                              const Variable&      mSq23  ,
                              const Variable&      t      ,
                              const ParameterExpr& width  ,
                              const Amplitude&     amp    ,
                              const CoefExpr&      z      ,
                              const CoefExpr&      qoverp ,
                              const PhaseSpace&    ps     ,
                              bool                 docache  )
  : DecayModel( mSq12, mSq13, mSq23, amp, ps ),
    _mSq12( mSq12.name() ), _mSq13( mSq13.name() ), _mSq23( mSq23.name() ), _t( t.name() ),
    _width( width ),
    _z( z ),
    _qoverp( qoverp ),
    _hasMixing( true ),
    _hasCPV   ( true ),
    _nDir( 0.0 ), _nCnj( 0.0 ), _nXed( 0.0 ), _norm( 1.0 ), _fixedAmp( false ),
    _maxPdf( 54.0 ), _cacheAmps( false ), _ampDirCache( 0 ), _ampCnjCache( 0 )
{
  // Make the variables available to cfit.
  // The squared invariant masses are already made available by the DecayModel constructor.
  push( t     );

  // Make all the parameters in width, z and qoverp available to cfit.
  push( width  );
  push( z      );
  push( qoverp );

  // Do calculations common to all values of variables (compute norm).
  if ( docache )
    cache();
}




// Copy function using default copy constructor.
Decay3BodyMix* Decay3BodyMix::copy() const
{
  return new Decay3BodyMix( *this );
}


// Time evolution functions.
const double Decay3BodyMix::psip( const double& t ) const
{
  return std::exp( - ( 1.0 - x() ) * gamma() * t );
}


const double Decay3BodyMix::psim( const double& t ) const
{
  return std::exp( - ( 1.0 + x() ) * gamma() * t );
}


const std::complex< double > Decay3BodyMix::psii( const double& t ) const
{
  return std::exp( - std::complex< double >( 1.0, - y() ) * gamma() * t );
}



// Setters for the mixing and CP violation parameters.
void Decay3BodyMix::setMixing( const CoefExpr& z )
{
  if ( _hasMixing )
    throw PdfException( "Decay3BodyMix: mixing coefficient z has already been set." );

  _hasMixing = true;
  push( _z = z );
}

void Decay3BodyMix::setCPV( const CoefExpr& qoverp )
{
  if ( _hasCPV )
    throw PdfException( "Decay3BodyMix: mixing coefficient z has already been set." );

  _hasCPV = true;
  push( _qoverp = qoverp );
}



void Decay3BodyMix::setParExpr()
{
  _width.setPars( _parMap );

  if ( _hasMixing ) _z     .setPars( _parMap );
  if ( _hasCPV    ) _qoverp.setPars( _parMap );
}



const double Decay3BodyMix::evaluateFuncs( const double& mSq12, const double& mSq13, const double& mSq23 ) const
{
  double value = 1.0;

  typedef std::vector< Function >::const_iterator fIter;

  std::map< std::string, double > varMap;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
  {
    if ( func->dependsOn( _mSq12 ) ) varMap[ _mSq12 ] = mSq12;
    if ( func->dependsOn( _mSq13 ) ) varMap[ _mSq13 ] = mSq13;
    if ( func->dependsOn( _mSq23 ) ) varMap[ _mSq23 ] = mSq23;

    value *= func->evaluate( varMap );
  }

  // Always return a non-negative value. Default to zero.
  return std::max( value, 0.0 );
}


const double Decay3BodyMix::evaluateFuncs( const double& mSq12, const double& mSq13 ) const
{
  const double& mSq23 = _ps.mSqSum() - mSq12 - mSq13;

  return evaluateFuncs( mSq12, mSq13, mSq23 );
}




void Decay3BodyMix::cacheNormComponents()
{
  // If the amplitude is fixed and the components have already
  //    been computed, just return without recomputing anything.
  if ( _fixedAmp )
    return;

  // Compute the value of _norm.
  _nDir = 0.0;
  _nCnj = 0.0;
  _nXed = 0.0;

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
}



void Decay3BodyMix::cache()
{
  // Compute the norm components, only if the amplitude
  //    is not fixed or they have not yet been computed.
  cacheNormComponents();

  _fixedAmp = _amp.isFixed();

  // Calculate the norm.
  const double&& xval = x();
  const double&& yval = y();
  const double&& xSq  = std::pow( xval, 2 );
  const double&& ySq  = std::pow( yval, 2 );

  const std::complex< double >&& qoverp = _qoverp.evaluate();
  const double&& cpvp                   = ( 1.0 + std::norm( qoverp ) ) / 2.0;
  const double&& cpvm                   = ( 1.0 - std::norm( qoverp ) ) / 2.0;

  _norm  = ( _nDir * cpvp + xval * std::real( _nXed * qoverp ) ) / ( 1.0 - xSq );
  _norm += ( _nDir * cpvm - yval * std::imag( _nXed * qoverp ) ) / ( 1.0 + ySq );
  _norm /= gamma();

  return;
}



const std::map< unsigned, std::vector< std::complex< double > > > Decay3BodyMix::cacheComplex( const Dataset& data )
{
  // Determine whether the amplitudes should be cached, i.e. only if all their parameters are fixed.
  _cacheAmps = _amp.isFixed();

  std::map< unsigned, std::vector< std::complex< double > > > cached;

  if ( ! _cacheAmps )
    return cached;

  // Get an index for the cached complex amplitudes.
  _ampDirCache = _cacheIdxComplex++;
  _ampCnjCache = _cacheIdxComplex++;

  double mSq12;
  double mSq13;
  double mSq23;

  // Cache the direct and conjugated amplitudes for every point in the given dataset.
  const std::size_t& size = data.size();
  for ( std::size_t entry = 0; entry < size; ++entry )
  {
    mSq12 = data.value( _mSq12, entry );
    mSq13 = data.value( _mSq13, entry );
    mSq23 = data.value( _mSq23, entry );

    cached[ _ampDirCache ].push_back( _amp.evaluate( _ps, mSq12, mSq13, mSq23 ) );
    cached[ _ampCnjCache ].push_back( _amp.evaluate( _ps, mSq13, mSq12, mSq23 ) );
  }

  return cached;
}


// Unnormalized evaluation. Calculate:
// | A + Abar |^2            | A* - Abar* |^2                [ A + Abar   A* - Abar*            ]
// | -------- |   psi_+(t) + | ---------- |   psi_-(t) + 2 Re[ -------- * ---------- * psi_i(t) ]
// |    2     |              |     2      |                  [    2           2                 ]
const double Decay3BodyMix::evaluateUnnorm( const double& mSq12, const double& mSq13, const double& mSq23, const double& t ) const
{
  // Particle decay amplitude.
  std::complex< double > ampDir = _amp.evaluate( _ps, mSq12, mSq13, mSq23 );
  std::complex< double > ampCnj = _amp.evaluate( _ps, mSq13, mSq12, mSq23 );
  if ( _hasCPV )
    ampCnj *= _qoverp.evaluate();

  const std::complex< double >&& apb2 =          ( ampDir + ampCnj ) / 2.0;
  const std::complex< double >&& amb2 = std::conj( ampDir - ampCnj ) / 2.0;

  double ampSq = 0.0;
  ampSq += std::norm( apb2 ) * psip( t );
  ampSq += std::norm( amb2 ) * psim( t );
  ampSq += 2.0 * std::real( apb2 * amb2 * psii( t ) );

  return ampSq * evaluateFuncs( mSq12, mSq13, mSq23 );
}


const double Decay3BodyMix::evaluate( const double& mSq12, const double& mSq13, const double& mSq23, const double& t ) const
{
  return evaluateUnnorm( mSq12, mSq13, mSq23, t ) / _norm;
}


const double Decay3BodyMix::evaluate( const double& mSq12, const double& mSq13, const double& t ) const
{
  const double& mSq23 = _ps.mSqSum() - mSq12 - mSq13;

  return evaluateUnnorm( mSq12, mSq13, mSq23, t ) / _norm;
}


const double Decay3BodyMix::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  const std::size_t& size = vars.size();

  if ( size == 3 )
    return evaluate( vars[ 0 ], vars[ 1 ], vars[ 2 ] );

  if ( size == 4 )
    return evaluate( vars[ 0 ], vars[ 1 ], vars[ 2 ], vars[ 3 ] );

  throw PdfException( "Decay3BodyMix can only take either 3 or 4 arguments." );
}



const double Decay3BodyMix::evaluate( const std::vector< double >&                 vars  ,
                                      const std::vector< double >&                 cacheR,
                                      const std::vector< std::complex< double > >& cacheC ) const throw( PdfException )
{
  if ( ! _cacheAmps )
    return evaluate( vars );

  std::map< std::string, Variable >::const_iterator&& tpos = _varMap.find( _t );
  if ( tpos == _varMap.end() )
    throw PdfException( "Decay3BodyMix: model does not depend on required variable. This is a bug." );

  const double& t = vars[ std::distance( _varMap.begin(), tpos ) ];

  // Particle decay amplitude.
  std::complex< double > ampDir = cacheC[ _ampDirCache ];
  std::complex< double > ampCnj = cacheC[ _ampCnjCache ];
  if ( _hasCPV )
    ampCnj *= _qoverp.evaluate();

  const std::complex< double >&& apb2 =          ( ampDir + ampCnj ) / 2.0;
  const std::complex< double >&& amb2 = std::conj( ampDir - ampCnj ) / 2.0;

  // Calculate the squared amplitude.
  double ampSq = 0.0;
  ampSq += std::norm( apb2 ) * psip( t );
  ampSq += std::norm( amb2 ) * psim( t );
  ampSq += 2.0 * std::real( apb2 * amb2 * psii( t ) );

  // Evaluate the efficiency functions. Could be cached if necessary.
  double funcs;
  const std::size_t& size = vars.size();
  if ( size == 3 )
    funcs = evaluateFuncs( vars[ 0 ], vars[ 1 ] );
  else if ( size == 4 )
    funcs = evaluateFuncs( vars[ 0 ], vars[ 1 ], vars[ 2 ] );
  else
    throw PdfException( "Decay3BodyMix can only take either 3 or 4 arguments." );

  return ampSq * funcs / _norm;
}



const double Decay3BodyMix::project( const std::string& varName, const double& x ) const throw( PdfException )
{
  throw PdfException( "Decay3BodyMix::project is not implemented yet" );

  return 0.0;
}


// No need to append an operator, since it can only be multiplication.
const Decay3BodyMix& Decay3BodyMix::operator*=( const Function& right )
{
  // Check that the function does not depend on any variables that the model does not.
  const std::map< std::string, Variable >& varMap = right.getVarMap();
  for ( std::map< std::string, Variable >::const_iterator var = varMap.begin(); var != varMap.end(); ++var )
    if ( ! _varMap.count( var->second.name() ) )
      throw PdfException( "Cannot multiply a Decay3BodyMix pdf model by a function that depends on other variables." );

  // Consider the function parameters as own ones.
  const std::map< std::string, Parameter >& parMap = right.getParMap();
  _parMap.insert( parMap.begin(), parMap.end() );

  // Append the function to the functions vector.
  _funcs.push_back( right );

  // Recompute the norm, since the pdf shape has changed under this operation.
  cache();

  return *this;
}



const Decay3BodyMix operator*( Decay3BodyMix left, const Function& right )
{
  // Check that the function does not depend on any variables that the model does not.
  const std::map< std::string, Variable >& varMap = right.getVarMap();
  for ( std::map< std::string, Variable >::const_iterator var = varMap.begin(); var != varMap.end(); ++var )
    if ( ! left._varMap.count( var->second.name() ) )
      throw PdfException( "Cannot multiply a Decay3BodyMix pdf model by a function that depends on other variables." );

  // Consider the function parameters as own ones.
  const std::map< std::string, Parameter >& parMap = right.getParMap();
  left._parMap.insert( parMap.begin(), parMap.end() );

  // Append the function to the functions vector.
  left._funcs.push_back( right );

  // Recompute the norm, since the pdf shape has changed under this operation.
  left.cache();

  return left;
}



const Decay3BodyMix operator*( const Function& left, Decay3BodyMix right )
{
  // Check that the function does not depend on any variables that the model does not.
  const std::map< std::string, Variable >& varMap = left.getVarMap();
  for ( std::map< std::string, Variable >::const_iterator var = varMap.begin(); var != varMap.end(); ++var )
    if ( ! right._varMap.count( var->second.name() ) )
      throw PdfException( "Cannot multiply a Decay3BodyMix pdf model by a function that depends on other variables." );

  // Consider the function parameters as own ones.
  const std::map< std::string, Parameter >& parMap = left.getParMap();
  right._parMap.insert( parMap.begin(), parMap.end() );

  // Append the function to the functions vector.
  right._funcs.push_back( left );

  // Recompute the norm, since the pdf shape has changed under this operation.
  right.cache();

  return right;
}



const std::map< std::string, double > Decay3BodyMix::generate() const throw( PdfException )
{
  // Generate mSq12 and mSq13, and compute mSq23 from these.
  const double& min12 = std::pow( _ps.m1()      + _ps.m2(), 2 );
  const double& min13 = std::pow( _ps.m1()      + _ps.m3(), 2 );
  const double& max12 = std::pow( _ps.mMother() - _ps.m3(), 2 );
  const double& max13 = std::pow( _ps.mMother() - _ps.m2(), 2 );

  // Maximum value of the pdf.
  const double& max = _maxPdf * tau();

  // Sum of squared invariant masses of all particles (mother and daughters).
  const double& mSqSum = _ps.mSqSum();

  double mSq12  = 0.0;
  double mSq13  = 0.0;
  double mSq23  = 0.0;
  double gammat = 0.0;
  double time   = 0.0;

  double pdfVal = 0.0;

  // Attempts to generate an event.
  unsigned count = 10000;

  std::map< std::string, double > values;

  while ( count-- )
  {
    // Generate uniform mSq12 and mSq13.
    mSq12 = Random::flat( min12, max12 );
    mSq13 = Random::flat( min13, max13 );
    mSq23 = mSqSum - mSq12 - mSq13;

    // Generate time according to q(t) = e^[ - ( 1 - |x| ) Gamma t ]
    gammat = - std::log( Random::flat() ) / ( 1.0 - std::abs( x() ) );
    time   = gammat / gamma();

    // Prepare to run accept-reject on p(t) / q(t).
    pdfVal = this->evaluate( mSq12, mSq13, mSq23, time ) * std::exp( ( 1.0 - std::abs( x() ) ) * gammat );

    // Complain if the maximum value of the pdf was not properly estimated.
    if ( pdfVal > max )
      std::cout << "Problem: " << pdfVal << " > " << max
                << " for ( " << mSq12 << ", " << mSq13 << ", " << mSq23 << ", " << time << " )" << std::endl;

    // Apply the accept-reject decision.
    if ( Random::flat( 0.0, max ) < pdfVal )
    {
      values[ _mSq12 ] = mSq12;
      values[ _mSq13 ] = mSq13;
      values[ _mSq23 ] = mSq23;
      values[ _t     ] = time;
      return values;
    }
  }

  std::cout << "Problem: too many attempts" << std::endl;
  values[ _mSq12 ] = 0.0;
  values[ _mSq13 ] = 0.0;
  values[ _mSq23 ] = 0.0;
  values[ _t     ] = 0.0;

  return values;
}

