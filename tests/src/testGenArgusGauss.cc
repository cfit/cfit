#include <iostream>
#include <fstream>

#include <cfit/parameter.hh>
#include <cfit/variable.hh>
#include <cfit/dataset.hh>
#include <cfit/chi2.hh>
#include <cfit/pdf.hh>

#include <cfit/models/genargusgauss.hh>

#include <Minuit/MnPrint.h>

#define MIN   (  -1.0   )
#define MAX   (   4.5   )
#define NBIN  ( 200     )
#define NEVT  (   1.0e6 )

#define RANGE ( MAX - MIN )
#define STEP  ( RANGE / double( NBIN ) )
#define AREA  ( RANGE * NEVT / double( NBIN ) ) // Expected value of the area.


// Auxiliary functions to compute the bin number from a point, the center of
//    a bin, the (approximate) Poisson error, and to read a data file.
int    bin      ( double x ) { return int( NBIN / RANGE * ( x - MIN ) ); }
double binCenter( int    b ) { return RANGE / NBIN * ( b + .5 ) + MIN;   }
double error    ( double y ) { return y == 0 ? 1. : sqrt( y ); }
void   readData ( std::string fileName, Dataset& data )
{
  double x;
  double y[ NBIN ];

  // Initialize the bin contents to zero.
  for ( int b = 0; b < NBIN; b++ )
    y[ b ] = 0;

  // Open the input file.
  std::ifstream file( fileName.c_str() );
  if ( ! file )
    throw std::exception();

  // Read the input file.
  while ( file >> x )
    if ( ( x > MIN ) && ( x < MAX ) )
      y[ bin( x ) ]++;

  // Close the input file.
  file.close();

  // Fill the dataset with the histogram.
  for ( int b = 0; b < NBIN; b++ )
  {
    data.push( "x", binCenter( b ) );
    data.push( "y", y[ b ], error( y[ b ] ) );
  }

  return;
}



int main( int argc, char** argv )
{
  if ( argc < 2 )
  {
    std::cerr << "\nUsage: " << argv[ 0 ] << " [gen|fit]\n" << std::endl;
    return 1;
  }

  // Variables the model depends on.
  Variable x( "x" );
  Variable y( "y" );

  // Parameters of the model.
  Parameter area ( "area",  2.0e5, 300.0   );
  Parameter c    ( "c"    , 3.0, 0.003 );
  Parameter chi  ( "chi"  , 2.0, 0.003 );
  Parameter p    ( "p"    , 1.4, 0.003 );
  Parameter mu   ( "mu"   , 0.0, 0.001 );
  Parameter sigma( "sigma", 0.5, 0.001 );

  p    .fix();
  mu   .fix();
  sigma.fix();

  // Definition of the pdf.
  GenArgusGauss argus( x, c, chi, p, mu, sigma );

  if ( std::string( argv[ 1 ] ) == "gen" )
  {
    // Open the output file.
    std::ofstream output( "data/genargusgauss.dat" );

    // Fill the output file with randomly generated values.
    std::map< std::string, double > entry;
    for ( unsigned k = 0; k < NEVT; ++k )
    {
      entry = argus.generate();
      for ( std::map< std::string, double >::const_iterator key = entry.begin(); key != entry.end(); ++key )
        output << key->second << std::endl;
    }

    // Close the output file.
    output.close();
  }

  if ( std::string( argv[ 1 ] ) == "fit" )
  {
    // Data container.
    Dataset data;

    // Read the data.
    readData( "data/genargusgauss.dat", data );

    // Define the pdf to be fitted.
    Pdf sum = area * argus;

    // Definition of the minimizer from the pdf.
    Chi2 chi2( sum, y, data );

    // Compute the minimum.
    const FunctionMinimum& min = chi2.minimize();

    // Output the result.
    std::cout << min << std::endl;
  }

  return 0;
}
