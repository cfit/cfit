
#include <stack>
#include <algorithm>
#include <functional>

#include <cfit/phasespace.hh>
#include <cfit/binnedamplitude.hh>
#include <cfit/coef.hh>
#include <cfit/resonance.hh>
#include <cfit/operation.hh>



BinnedAmplitude::BinnedAmplitude( const std::vector< Parameter >& tpb,
                                  const std::vector< Parameter >& tmb,
                                  const std::vector< CoefExpr  >& xb  )
  : _nbins( tpb.size() ),
    _tpb( tpb ), _tmb( tmb ), _xb( xb )
{
  if ( _tmb.size() != _nbins )
    throw PdfException( "BinnedAmplitude: squared amplitude vectors must have the same size." );

  if ( _xb.size() != _nbins )
    throw PdfException( "BinnedAmplitude: squared amplitude and crossed amplitude vectors must have the same size." );

  for ( std::vector< Parameter >::const_iterator tp = tpb.begin(); tp != tpb.end(); ++tp )
    _parMap.emplace( tp->name(), *tp );

  for ( std::vector< Parameter >::const_iterator tm = tmb.begin(); tm != tmb.end(); ++tm )
    _parMap.emplace( tm->name(), *tm );

  for ( std::vector< CoefExpr >::const_iterator x = xb.begin(); x != xb.end(); ++x )
  {
    const std::map< std::string, Parameter >& pars = x->getPars();
    _parMap.insert( pars.begin(), pars.end() );

    // _parMap.emplace( x->real().name(), x->real() );
    // _parMap.emplace( x->imag().name(), x->imag() );
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

