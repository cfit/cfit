
#include <iostream> // TEMPORAL

#include <vector>
#include <string>

#ifdef MPI_ON
#include <mpi.h>
#endif

#include <cfit/variable.hh>
#include <cfit/dataset.hh>
#include <cfit/pdfmodel.hh>
#include <cfit/chi2.hh>


Chi2::Chi2( const PdfModel& pdf, const Variable& y, const Dataset& data )
  : Minimizer( pdf, data ), _y( y )
{
  _up = 1.0;
}


Chi2::Chi2( const PdfExpr& pdf, const Variable& y, const Dataset& data )
  : Minimizer( pdf, data ), _y( y )
{
  _up = 1.0;
}


Chi2::Chi2( const Chi2& chi2 )
  : Minimizer( chi2 ), _y( chi2._y )
{}


double Chi2::operator()( const std::vector<double>& pars ) const throw( PdfException )
{
  if ( pars.size() != _pdf->nPars() )
    throw PdfException( "Number of parameters passed does not match number of required arguments." );

  _pdf->setPars( pars );

  // Before evaluating the pdf at all data points, cache anything common to
  //    all points (usually compute the norm).
  _pdf->cache();

  // Get the vector of variable names that the pdf depends on.
  std::vector< std::string > varNames = _pdf->varNames();
  typedef std::vector< std::string >::const_iterator vIter;

  // Vector of values of the variables that the pdf must be evaluated at.
  std::vector< double > vars;

  // Initialize the value of the chi^2.
  double chi2 = 0.;

  // Sum of the terms of the chi^2.
  for ( std::size_t n = 0; n < _data.size(); ++n )
    {
      // Initialize the value of the variance for the current entry.
      //    It must be s_y^2 + Sum( s_x^2 ).
      double variance = 0.;

      // Reset the vector of values of the variables.
      vars.clear();

      // Fill the vector of values and sum the terms of the variance.
      for ( vIter var = varNames.begin(); var != varNames.end(); ++var )
	{
	  vars.push_back( _data.value( *var, n ) );
	  variance += pow( _data.error( *var, n ), 2 );
	}

      // Compute the numerator of the chi^2 term and finish computing the variance.
      double diff = _pdf->evaluate( vars ) - _data.value( _y.name(), n );
      variance += pow( _data.error( _y.name(), n ), 2 );

      // Add the term to the chi^2.
      chi2 += pow( diff, 2 ) / variance;
    }

#ifdef MPI_ON
  // If running with MPI, each process has only computed a piece of the chi2.
  //    Add all the pieces up and broadcast them to all the processes.
  double result = 0.;
  MPI::Comm& world = MPI::COMM_WORLD;
  world.Barrier();
  world.Allreduce( &chi2, &result, 1, MPI::DOUBLE, MPI::SUM );

  return result;
#else
  std::cout << "chi2 = " << chi2 << std::endl;

  return chi2;
#endif
}

