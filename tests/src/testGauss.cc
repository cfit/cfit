#include <iostream>
#include <fstream>

#include <cfit/parameter.hh>
#include <cfit/variable.hh>
#include <cfit/dataset.hh>
#include <cfit/chi2.hh>
#include <cfit/pdf.hh>

#include <cfit/models/gauss.hh>

#include <Minuit/MnPrint.h>

#define MIN   ( -10.0   )
#define MAX   (  10.0   )
#define NBIN  ( 100.0   )
#define NEVT  (   1.0e6 )

#define RANGE ( MAX - MIN )
#define AREA  ( RANGE * NEVT / NBIN ) // Expected value of the area.


// Auxiliary functions to compute the bin number from a point, the center of
//    a bin, the (approximate) Poisson error, and to read a data file.
int    bin      ( double x ) { return int( NBIN / RANGE * ( x - MIN ) ); }
double binCenter( int    b ) { return RANGE / NBIN * ( b + .5 ) + MIN;   }
double error    ( double y ) { return y == 0 ? 1. : sqrt( y ); }
void   readData ( std::string fileName, Dataset& data )
{
  double x;
  double y[ 100 ];

  // Initialize the bin contents to zero.
  for ( int b = 0; b < 100; b++ )
    y[ b ] = 0;

  // Open the input file.
  std::ifstream file( fileName.c_str() );
  if ( ! file )
    throw std::exception();

  // Read the input file.
  while ( file >> x )
    y[ bin( x ) ]++;

  // Close the input file.
  file.close();

  // Fill the dataset with the histogram.
  for ( int b = 0; b < 100; b++ )
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
  Parameter area  ( "area"  ,  2.0e5, 300.0   );
  Parameter phi   ( "phi"   ,  0.5  ,   0.1   );
  Parameter mean1 ( "mean1" , -1.0  ,   0.003 );
  Parameter mean2 ( "mean2" ,  2.0  ,   0.003 );
  Parameter sigma1( "sigma1",  2.0  ,   0.003 );
  Parameter sigma2( "sigma2",  1.0  ,   0.003 );

  // Definition of the pdf.
  Gauss g1( x, mean1, sigma1 );
  Gauss g2( x, mean2, sigma2 );

  if ( std::string( argv[ 1 ] ) == "gen" )
  {
    Pdf sum = g1 + g2;

    // Open the output file.
    std::ofstream output( "data/gauss.dat" );

    // Fill the output file with randomly generated values.
    std::map< std::string, double > entry;
    for ( unsigned k = 0; k < 1000000; ++k )
    {
      entry = sum.generate();
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
    readData( "data/gauss.dat", data );

    // Define the pdf to be fitted.
    Pdf sum = area * ( pow( sin( phi ), 2 ) * g1 + pow( cos( phi ), 2 ) * g2 );

    // Definition of the minimizer from the pdf.
    Chi2 chi2( sum, y, data );

    // Compute the minimum.
    const FunctionMinimum& min = chi2.minimize();

    // Output the result.
    std::cout << min << std::endl;
  }

  return 0;
}
