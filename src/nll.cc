
#include <iostream>

#include <vector>
#include <string>

#ifdef MPI_ON
#include <mpi.h>
#endif

#include <cfit/variable.hh>
#include <cfit/dataset.hh>
#include <cfit/pdfmodel.hh>
#include <cfit/nll.hh>


Nll::Nll( PdfBase& pdf, const Dataset& data )
  : Minimizer( pdf, data )
{
  _up = 1.0;
}


double Nll::operator()( const std::vector<double>& pars ) const throw( PdfException )
{
  if ( pars.size() != _pdf.nPars() )
    throw PdfException( "Number of parameters passed does not match number of required arguments." );

  _pdf.setPars( pars );

  // Before evaluating the pdf at all data points, cache anything common to
  //    all points (usually compute the norm).
  _pdf.cache();

  // Get the vector of variable names that the pdf depends on.
  std::vector< std::string > varNames = _pdf.varNames();

  typedef std::vector< std::string                                     >::const_iterator vIter;
  typedef std::map   < unsigned, std::vector< double >                 >::const_iterator mrIter;
  typedef std::map   < unsigned, std::vector< std::complex< double > > >::const_iterator mcIter;

  // Vector of values of the variables that the pdf must be evaluated at, and vectors of cached values.
  std::vector< double                 > vars;
  std::vector< double                 > cacheR;
  std::vector< std::complex< double > > cacheC;

  // Allocate memory for the vectors of cached variables.
  cacheR.reserve( _pdf.nCachedReal()    );
  cacheC.reserve( _pdf.nCachedComplex() );

  // Initialize the value of the nll.
  double nll = 0.;

  double value = 0.;

  // Sum of the terms of the nll.
  for ( std::size_t n = 0; n < _data.size(); ++n )
  {
    // Reset the vector of values of the variables and the previously cached values.
    vars  .clear();
    cacheR.clear();
    cacheC.clear();

    // Fill the vector of values and sum the terms of the variance.
    for ( vIter var = varNames.begin(); var != varNames.end(); ++var )
      vars.push_back( _data.value( *var, n ) );

    for ( mrIter cached = _cacheR.begin(); cached != _cacheR.end(); ++cached )
      cacheR[ cached->first ] = cached->second[ n ];

    for ( mcIter cached = _cacheC.begin(); cached != _cacheC.end(); ++cached )
      cacheC[ cached->first ] = cached->second[ n ];

    // _pdf.setVars( vars );

    // Add the term to the nll.
    value = _pdf.evaluate( vars, cacheR, cacheC );
    if ( value )
      nll += - 2. * log( value );
//       else
// 	std::cout << "Warning: pdf evaluates to zero for entry " << n
// 		  << ". Not taking this entry into account for the nll." << std::endl;
  }

#ifdef MPI_ON
  // If running with MPI, each process has only computed a piece of the chi2.
  //    Add all the pieces up and broadcast them to all the processes.
  double result = 0.;
  MPI::Comm& world = MPI::COMM_WORLD;
  world.Allreduce( &nll, &result, 1, MPI::DOUBLE, MPI::SUM );

  return result;
#else
  std::cout << "nll = " << nll << std::endl;

  return nll;
#endif
}

