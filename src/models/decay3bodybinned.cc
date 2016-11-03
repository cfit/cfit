
#include <complex>

#include <cfit/dataset.hh>
#include <cfit/function.hh>
#include <cfit/random.hh>

#include <cfit/models/decay3bodybinned.hh>


Decay3BodyBinned::Decay3BodyBinned( const Variable&        bin    ,
                                    const BinnedAmplitude& amp    ,
                                    const CoefExpr&        phi    ,
                                    bool                   docache  )
  : _amp( amp ), _hasKappa( false ), _phi( phi ),
    _nDir( 0.0 ), _nXed( 0.0 ), _norm( 1.0 ), _fixedAmp( false ),
    _cacheAmps( false ), _ampDirCache( 0 ), _ampCnjCache( 0 )
{
  push( bin );

  push( amp );
  push( phi );

  // Do calculations common to all values of variables
  //    (usually compute norm).
  if ( docache )
    cache();
}



Decay3BodyBinned::Decay3BodyBinned( const Variable&         bin    ,
                                    const BinnedAmplitude&  amp    ,
                                    const CoefExpr&         phi    ,
                                    const Parameter&        kappa  ,
                                    bool                    docache  )
  : _amp( amp ), _hasKappa( true ), _kappa( kappa ), _phi( phi ),
    _nDir( 0.0 ), _nXed( 0.0 ), _norm( 1.0 ), _fixedAmp( false ),
    _cacheAmps( false ), _ampDirCache( 0 ), _ampCnjCache( 0 )
{
  push( bin );

  push( amp   );
  push( phi   );
  push( kappa );

  // Do calculations common to all values of variables
  //    (usually compute norm).
  if ( docache )
    cache();
}




Decay3BodyBinned::Decay3BodyBinned( const Variable&         bin    ,
                                    const BinnedAmplitude&  amp    ,
                                    const CoefExpr&         phi    ,
                                    const ParameterExpr&    kappa  ,
                                    bool                    docache  )
  : _amp( amp ), _hasKappa( true ), _kappa( kappa ), _phi( phi ),
    _nDir( 0.0 ), _nXed( 0.0 ), _norm( 1.0 ), _fixedAmp( false ),
    _cacheAmps( false ), _ampDirCache( 0 ), _ampCnjCache( 0 )
{
  push( bin );

  push( amp   );
  push( phi   );
  push( kappa );

  // Do calculations common to all values of variables
  //    (usually compute norm).
  if ( docache )
    cache();
}





Decay3BodyBinned* Decay3BodyBinned::copy() const
{
  return new Decay3BodyBinned( *this );
}



void Decay3BodyBinned::setParExpr()
{
  _phi.setPars( _parMap );

  // Propagate the kappa parameter value if necessary.
  if ( _hasKappa )
    _kappa.setPars( _parMap );
}



// const double Decay3BodyBinned::evaluateFuncs( const double& mSq12, const double& mSq13, const double& mSq23 ) const
// {
//   double value = 1.0;

//   const std::string& name12 = getVar( 0 ).name(); // mSq12
//   const std::string& name13 = getVar( 1 ).name(); // mSq13
//   const std::string& name23 = getVar( 2 ).name(); // mSq23

//   typedef std::vector< Function >::const_iterator fIter;

//   std::map< std::string, double > varMap;
//   for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
//   {
//     if ( func->dependsOn( name12 ) ) varMap[ name12 ] = mSq12;
//     if ( func->dependsOn( name13 ) ) varMap[ name13 ] = mSq13;
//     if ( func->dependsOn( name23 ) ) varMap[ name23 ] = mSq23;

//     value *= func->evaluate( varMap );
//   }

//   // Always return a non-negative value. Default to zero.
//   return std::max( value, 0.0 );
// }


// const double Decay3BodyBinned::evaluateFuncs( const double& mSq12, const double& mSq13 ) const
// {
//   const double& mSq23 = _ps.mSqSum() - mSq12 - mSq13;

//   return evaluateFuncs( mSq12, mSq13, mSq23 );
// };


// // The values of the variables must be set with setVars before using this function.
// const double Decay3BodyBinned::evaluateFuncs() const
// {
//   double value = 1.0;

//   typedef std::vector< Function >::const_iterator fIter;
//   for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
//     value *= func->evaluate();

//   // Always return a non-negative value. Default to zero.
//   return std::max( value, 0.0 );
// }



void Decay3BodyBinned::cacheNormComponents()
{
  // If the amplitude is fixed and the components have already
  //    been computed, just return without recomputing anything.
  if ( _fixedAmp )
    return;

  // Initialize the value of the norm components and the norm.
  _nDir = 0.0;
  _nXed = 0.0;
  _norm = 0.0;

  std::vector< Parameter >::const_iterator tpb = _amp.tpb().begin();
  std::vector< Parameter >::const_iterator tmb = _amp.tmb().begin();
  std::vector< CoefExpr  >::const_iterator xb  = _amp.xb() .begin();

  double                 tp;
  double                 tm;
  std::complex< double > tx;

  // Calculate the norm components.
  while ( tpb != _amp.tpb().end() )
  {
    tp = (tpb++)->value();
    tm = (tmb++)->value();
    tx = std::sqrt( tp * tm ) * (xb++)->evaluate();

    _nDir += tp + tm;
    _nXed += 2.0 * std::real( tx );
  }

  return;
}


void Decay3BodyBinned::cache()
{
  // Compute the norm components, only if the amplitude
  //    is not fixed or they have not yet been computed.
  cacheNormComponents();

  _fixedAmp = true; // _amp.isFixed();

  // Calculate the norm.
  const std::complex< double >& vz     = z();
  const double&                 vKappa = kappa();
  _norm = _nDir * ( 1.0 + std::norm( vz ) ) + 2.0 * vKappa * std::real( vz * _nXed );

  return;
}



// Unnormalized evaluation.
const double Decay3BodyBinned::evaluateUnnorm( const int& bin ) const throw( PdfException )
{
  const std::tuple< double, double, std::complex< double > >&& tx = _amp.evaluate( bin );

  const std::complex< double >&& vz     = z();
  const double&&                 vKappa = kappa();

  return std::get< 0 >( tx )                   +
         std::get< 1 >( tx ) * std::norm( vz ) +
         2.0 * vKappa * std::real( vz * std::get< 2 >( tx ) );
}


const double Decay3BodyBinned::evaluate( const int& bin ) const throw( PdfException )
{
  if ( bin )
    return evaluateUnnorm( bin ) / _norm;

  std::cout << "problema" << std::endl;

  return 0.0;
}




// const double Decay3BodyBinned::project( const std::string& varName, const double& x ) const throw( PdfException )
// {
//   // Find the index of the variable to be projected.
//   int index = -1;
//   for ( unsigned var = 0; var < 3; ++var )
//     if ( varName == getVar( var ).name() )
//       index = var;

//   // If the pdf does not depend on the passed variable name, the projection is 1.
//   if ( index == -1 )
//     return 1.0;

//   // Minimum and maximum values of the variable wrt which the integration
//   //    is done i.e. the next to the projected variable.
//   const double& min = _ps.mSqMin( ( index + 1 ) % 3 );
//   const double& max = _ps.mSqMax( ( index + 1 ) % 3 );

//   // Integrate the model.
//   const int& nbins = 400;
//   double proj = 0.0;
//   for ( int yBin = 0; yBin < nbins; ++yBin )
//   {
//     double y = binCenter( yBin, nbins, min, max );
//     double z = _ps.mSqSum() - x - y;

//     if ( index == 0 ) proj += evaluate( x, y, z );
//     if ( index == 1 ) proj += evaluate( z, x, y );
//     if ( index == 2 ) proj += evaluate( y, z, x );
//   }

//   proj *= ( max - min ) / double( nbins );

//   return proj;
// }


const double Decay3BodyBinned::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  if ( vars.size() != 1 )
    throw PdfException( "Decay3BodyBinned can only take 1 argument." );

  return evaluate( vars[ 0 ] );
}



const double Decay3BodyBinned::evaluate( const std::vector< double >&                 vars  ,
                                         const std::vector< double >&                 cacheR,
                                         const std::vector< std::complex< double > >& cacheC ) const throw( PdfException )
{
  return evaluate( vars );
}




// // No need to append an operator, since it can only be multiplication.
// const Decay3BodyBinned& Decay3BodyBinned::operator*=( const Function& right ) throw( PdfException )
// {
//   // Check that the function does not depend on any variables that the model does not.
//   const std::map< std::string, Variable >& varMap = right.getVarMap();
//   for ( std::map< std::string, Variable >::const_iterator var = varMap.begin(); var != varMap.end(); ++var )
//     if ( ! _varMap.count( var->second.name() ) )
//       throw PdfException( "Cannot multiply a Decay3BodyBinned pdf model by a function that depends on other variables." );

//   // Consider the function parameters as own ones.
//   const std::map< std::string, Parameter >& parMap = right.getParMap();
//   _parMap.insert( parMap.begin(), parMap.end() );

//   // Append the function to the functions vector.
//   _funcs.push_back( right );

//   // Recompute the norm, since the pdf shape has changed under this operation.
//   cache();

//   return *this;
// }



// const Decay3BodyBinned operator*( Decay3BodyBinned left, const Function& right )
// {
//   // Check that the function does not depend on any variables that the model does not.
//   const std::map< std::string, Variable >& varMap = right.getVarMap();
//   for ( std::map< std::string, Variable >::const_iterator var = varMap.begin(); var != varMap.end(); ++var )
//     if ( ! left._varMap.count( var->second.name() ) )
//       throw PdfException( "Cannot multiply a Decay3BodyBinned pdf model by a function that depends on other variables." );

//   // Consider the function parameters as own ones.
//   const std::map< std::string, Parameter >& parMap = right.getParMap();
//   left._parMap.insert( parMap.begin(), parMap.end() );

//   // Append the function to the functions vector.
//   left._funcs.push_back( right );

//   // Recompute the norm, since the pdf shape has changed under this operation.
//   left.cache();

//   return left;
// }



// const Decay3BodyBinned operator*( const Function& left, Decay3BodyBinned right )
// {
//   // Check that the function does not depend on any variables that the model does not.
//   const std::map< std::string, Variable >& varMap = left.getVarMap();
//   for ( std::map< std::string, Variable >::const_iterator var = varMap.begin(); var != varMap.end(); ++var )
//     if ( ! right._varMap.count( var->second.name() ) )
//       throw PdfException( "Cannot multiply a Decay3BodyBinned pdf model by a function that depends on other variables." );

//   // Consider the function parameters as own ones.
//   const std::map< std::string, Parameter >& parMap = left.getParMap();
//   right._parMap.insert( parMap.begin(), parMap.end() );

//   // Append the function to the functions vector.
//   right._funcs.push_back( left );

//   // Recompute the norm, since the pdf shape has changed under this operation.
//   right.cache();

//   return right;
// }


// const std::map< std::string, double > Decay3BodyBinned::generate() const throw( PdfException )
// {
//   // Generate mSq12 and mSq13, and compute mSq23 from these.
//   const double& min12 = std::pow( _ps.m1()      + _ps.m2(), 2 );
//   const double& min13 = std::pow( _ps.m1()      + _ps.m3(), 2 );
//   const double& max12 = std::pow( _ps.mMother() - _ps.m3(), 2 );
//   const double& max13 = std::pow( _ps.mMother() - _ps.m2(), 2 );

//   const std::string& mSq12name = getVar( 0 ).name();
//   const std::string& mSq13name = getVar( 1 ).name();
//   const std::string& mSq23name = getVar( 2 ).name();

//   // Sum of squared invariant masses of all particles (mother and daughters).
//   const double& mSqSum = _ps.mSqMother() + _ps.mSq1() + _ps.mSq2() + _ps.mSq3();

//   double mSq12 = 0.0;
//   double mSq13 = 0.0;
//   double mSq23 = 0.0;

//   double pdfVal  = 0.0;

//   // Attempts to generate an event.
//   int count = 10000;

//   std::map< std::string, double > values;

//   while ( count-- )
//   {
//     // Generate uniform mSq12 and mSq13.
//     mSq12 = Random::flat( min12, max12 );
//     mSq13 = Random::flat( min13, max13 );
//     mSq23 = mSqSum - mSq12 - mSq13;

//     values[ mSq12name ] = mSq12;
//     values[ mSq13name ] = mSq13;
//     values[ mSq23name ] = mSq23;

//     pdfVal = this->evaluate( mSq12, mSq13, mSq23 );

//     if ( pdfVal > _maxPdf )
//       std::cout << "Problem: " << pdfVal << " > " << _maxPdf
//                 << " for variables " << mSq12 << " " << mSq13 << " " << mSq23 << std::endl;

//     // Apply the accept-reject decision.
//     if ( Random::flat( 0.0, _maxPdf ) < pdfVal )
//       return values;
//   }

//   values[ mSq12name ] = 0.0;
//   values[ mSq13name ] = 0.0;
//   values[ mSq23name ] = 0.0;

//   return values;
// }
