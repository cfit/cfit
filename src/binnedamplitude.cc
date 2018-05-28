
#include <stack>
#include <algorithm>
#include <functional>

#include <cfit/phasespace.hh>
#include <cfit/binnedamplitude.hh>
#include <cfit/coef.hh>
#include <cfit/resonance.hh>
#include <cfit/operation.hh>



BinnedAmplitude::BinnedAmplitude( const std::vector< Parameter >& npb,
                                  const std::vector< Parameter >& nmb,
                                  const std::vector< CoefExpr  >& xb  )
  : _nbins( npb.size() ),
    _xb( xb )
{
  for ( std::vector< Parameter >::const_iterator par = npb.begin(); par != npb.end(); ++par )
    _npb.push_back( *par );

  for ( std::vector< Parameter >::const_iterator par = nmb.begin(); par != nmb.end(); ++par )
    _nmb.push_back( *par );

  if ( _nmb.size() != _nbins )
    throw PdfException( "BinnedAmplitude: squared amplitude vectors must have the same size." );

  if ( _xb.size() != _nbins )
    throw PdfException( "BinnedAmplitude: squared amplitude and crossed amplitude vectors must have the same size." );

  for ( std::vector< Parameter >::const_iterator np = npb.begin(); np != npb.end(); ++np )
    _parMap.emplace( np->name(), *np );

  for ( std::vector< Parameter >::const_iterator nm = nmb.begin(); nm != nmb.end(); ++nm )
    _parMap.emplace( nm->name(), *nm );

  for ( std::vector< CoefExpr >::const_iterator x = xb.begin(); x != xb.end(); ++x )
  {
    const std::map< std::string, Parameter >& pars = x->getPars();
    _parMap.insert( pars.begin(), pars.end() );

    // _parMap.emplace( x->real().name(), x->real() );
    // _parMap.emplace( x->imag().name(), x->imag() );
  }
}




BinnedAmplitude::BinnedAmplitude( const std::vector< ParameterExpr >& npb,
                                  const std::vector< ParameterExpr >& nmb,
                                  const std::vector< CoefExpr      >& xb  )
  : _nbins( npb.size() ),
    _npb( npb ), _nmb( nmb ), _xb( xb )
{
  if ( _nmb.size() != _nbins )
    throw PdfException( "BinnedAmplitude: squared amplitude vectors must have the same size." );

  if ( _xb.size() != _nbins )
    throw PdfException( "BinnedAmplitude: squared amplitude and crossed amplitude vectors must have the same size." );

  for ( std::vector< ParameterExpr >::const_iterator np = npb.begin(); np != npb.end(); ++np )
  {
    const std::map< std::string, Parameter >& pars = np->getPars();
    _parMap.insert( pars.begin(), pars.end() );
  }

  for ( std::vector< ParameterExpr >::const_iterator nm = nmb.begin(); nm != nmb.end(); ++nm )
  {
    const std::map< std::string, Parameter >& pars = nm->getPars();
    _parMap.insert( pars.begin(), pars.end() );
  }

  for ( std::vector< CoefExpr >::const_iterator x = xb.begin(); x != xb.end(); ++x )
  {
    const std::map< std::string, Parameter >& pars = x->getPars();
    _parMap.insert( pars.begin(), pars.end() );
  }
}



void BinnedAmplitude::setPars( const std::map< std::string, Parameter >& pars )
{
  typedef std::map< std::string, Parameter >::iterator mIter;
  for ( mIter par = _parMap.begin(); par != _parMap.end(); ++par )
    par->second.setValue( pars.find( par->first )->second.value() );

  // // Set the values of the parameters.
  // typedef std::vector< Parameter >::iterator pIter;
  // for ( pIter par = _parms.begin(); par != _parms.end(); ++par )
  //   par->setValue( pars.find( par->name() )->second.value() );

  // // Set the values of the coefficients.
  // typedef std::vector< Coef >::iterator cIter;
  // for ( cIter coef = _coefs.begin(); coef != _coefs.end(); ++coef )
  //   coef->setValue( pars.find( coef->real().name() )->second.value(),
  //                   pars.find( coef->imag().name() )->second.value() );
}

