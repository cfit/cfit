
#include <iostream>
#include <vector>

#include <Minuit/MnMigrad.h>

#include <cfit/minimizer.hh>
#include <cfit/variable.hh>
#include <cfit/parameter.hh>


FunctionMinimum Minimizer::minimize() const
{
  // Work with Minuit user defined parameters.
  MnUserParameters upar;

  typedef std::map< std::string, Parameter >::const_iterator pIter;

  // Set the Minuit parameters' name, value and uncertainty.
  const std::map< std::string, Parameter >& pars = _pdf.getParameters();
  for ( pIter par = pars.begin(); par != pars.end(); ++par )
    upar.add( par->first.c_str(), par->second.value(), par->second.error() );

  // Fix the parameters that are set to be fixed.
  for ( pIter par = pars.begin(); par != pars.end(); ++par )
    if ( par->second.isFixed() )
      upar.fix( par->first.c_str() );

  // Set the limits for those parameters that have some.
  for ( pIter par = pars.begin(); par != pars.end(); ++par )
    if ( par->second.hasLimits() )
      upar.setLimits( par->first.c_str(), par->second.lower(), par->second.upper() );

  MnMigrad migrad( *this, upar );

  return migrad();
}
