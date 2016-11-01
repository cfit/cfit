
#include <vector>

#include <Minuit/MnMigrad.h>

#include <cfit/minimizer.hh>
#include <cfit/variable.hh>
#include <cfit/parameter.hh>


void Minimizer::cache()
{
  const std::map< unsigned, std::vector< double >                 >& cachedR = _pdf->cacheReal   ( _data );
  const std::map< unsigned, std::vector< std::complex< double > > >& cachedC = _pdf->cacheComplex( _data );

  _cacheR.insert( cachedR.begin(), cachedR.end() );
  _cacheC.insert( cachedC.begin(), cachedC.end() );
}



FunctionMinimum Minimizer::minimize() const
{
  // Work with Minuit user defined parameters.
  MnUserParameters upar;

  typedef std::map< std::string, Parameter >::const_iterator pIter;

  // Set the Minuit parameters' name, value and uncertainty.
  const std::map< std::string, Parameter >& pars = _pdf->getPars();

  for ( pIter par = pars.begin(); par != pars.end(); ++par )
  {
    upar.add( par->first.c_str(), par->second.value(), par->second.error() );

    // Fix the parameters that are set to be fixed.
    if ( par->second.isFixed() )
      upar.fix( par->first.c_str() );

    // Set the blinding if requested.
    if ( par->second.isBlind() )
      upar.blind( par->first.c_str() );

    // Set the limits for those parameters that have some.
    if ( par->second.hasLimits() )
      upar.setLimits( par->first.c_str(), par->second.lower(), par->second.upper() );
  }

  MnMigrad migrad( *this, upar );

  return migrad();
}


double Minimizer::up() const throw( MinimizerException )
{
  if ( _up < 0.0 )
    throw MinimizerException( "The minimizer variation that specifies the sigma level of uncertainties must be positive." );

  return _up;
}

