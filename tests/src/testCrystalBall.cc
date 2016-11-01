#include <iostream>
#include <fstream>

#include <cfit/parameter.hh>
#include <cfit/variable.hh>
#include <cfit/dataset.hh>
#include <cfit/chi2.hh>
#include <cfit/pdfexpr.hh>

#include <cfit/models/crystalball.hh>

#include <Minuit/MnPrint.h>

#define MIN   ( -10.0   )
#define MAX   (  10.0   )
#define NBIN  ( 100     )
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
  Parameter area ( "area" ,  2.0e5, 300.0   );
  Parameter mu   ( "mu"   ,  0.0  ,   0.003 );
  Parameter sigma( "sigma",  2.0  ,   0.003 );
  Parameter alpha( "alpha", -0.7  ,   0.003 );
  Parameter n    ( "n"    ,  1.7  ,   0.003 );

  // Definition of the pdf.
  CrystalBall cb( x, mu, sigma, alpha, n );

  cb.setLimits( -8.0, 8.0 );

  if ( std::string( argv[ 1 ] ) == "gen" )
  {
    // Open the output file.
    std::ofstream output( "data/crystalball.dat" );

    // Fill the output file with randomly generated values.
    std::map< std::string, double > entry;
    for ( unsigned k = 0; k < NEVT; ++k )
    {
      entry = cb.generate();
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
    readData( "data/crystalball.dat", data );

    // Define the pdf to be fitted.
    PdfExpr sum = area * cb;

    // Definition of the minimizer from the pdf.
    Chi2 chi2( sum, y, data );

    // Compute the minimum.
    const FunctionMinimum& min = chi2.minimize();

    // Output the result.
    std::cout << min << std::endl;
  }

  return 0;
}
