
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <algorithm>

#include <cfit/binning.hh>
#include <cfit/exceptions.hh>


float Binning::distance( const std::pair< float, float >& a, const std::pair< float, float >& b )
{
  return std::sqrt( std::pow( a.first - b.first, 2 ) + std::pow( a.second - b.second, 2 ) );
}


float Binning::distance( const std::pair< float, float >& a, const float& bx, const float& by )
{
  return std::sqrt( std::pow( a.first - bx, 2 ) + std::pow( a.second - by, 2 ) );
}


Binning::Binning( const Data& xybin )
  : _xybinx( xybin ), _xybiny( xybin )
{
  std::sort( _xybinx.begin(), _xybinx.end(), ByX() );
  std::sort( _xybiny.begin(), _xybiny.end(), ByY() );
}


// Find the bin corresponding to position (x,y) over the phase space.
int Binning::bin( const float& x, const float& y ) const
{
  if ( _xybinx.empty() || _xybiny.empty() )
    throw PdfException( "Binning: cannot identify bin with an empty binning scheme." );

  if ( x < y )
    return - bin( y, x );

  // Find the first elements with x-coordinate >= x. Can be more than 1 if
  //    some points have the same x-coordinate and different y-coordinate.
  Data::const_iterator pxmin = std::lower_bound( _xybinx.begin(), _xybinx.end(), x                 , PosX() );
  Data::const_iterator pxmax = std::upper_bound( _xybinx.begin(), _xybinx.end(), pxmin->first.first, PosX() );

  // Out of these, find the first point with y-coordinate >= y.
  Data::const_iterator px = std::lower_bound( pxmin, pxmax, y, PosY() );

  // This point is likely to be close to (x,y). Find an upper limit to the closest point distance.
  float delta = distance( px->first, x, y );
  if ( px != _xybinx.begin() )
  {
    delta = std::min( delta, distance( ( px - 1 )->first, x, y ) );
    pxmin = std::lower_bound( _xybinx.begin(), _xybinx.end(), ( pxmin - 1 )->first.first, PosX() );
    pxmax = std::upper_bound( _xybinx.begin(), _xybinx.end(),   pxmin      ->first.first, PosX() );
    px = std::lower_bound( pxmin, pxmax, y, PosY() );
    delta = std::min( delta, distance( px->first, x, y ) );
    if ( px != _xybinx.begin() )
      delta = std::min( delta, distance( ( px - 1 )->first, x, y ) );
  }

  // Find all points with y-coordinate not further away than delta from the given point.
  Data::const_iterator pymin = std::lower_bound( _xybiny.begin(), _xybiny.end(), y - delta, PosY() );
  Data::const_iterator pymax = std::upper_bound( _xybiny.begin(), _xybiny.end(), y + delta, PosY() );

  Data::const_iterator py;

  float dist;
  Data::const_iterator point;
  while( pymin != pymax )
  {
    // Build (pymin, py) as a range of points with equal y, sorted by x.
    //    Find the two points closest to (x,y).
    py    = std::upper_bound( pymin, pymax, pymin->first.second, PosY() );
    pymin = std::lower_bound( pymin, py   , x                  , PosX() );
    if ( std::fabs( ( pymin - 1 )->first.first - x ) <= delta )
    {
      dist = distance( ( pymin - 1 )->first, x, y );
      if ( dist <= delta )
      {
        delta = dist;
        point = pymin - 1;
      }
    }
    if ( std::fabs( pymin->first.first - x ) <= delta )
    {
      dist = distance( pymin->first, x, y );
      if ( dist <= delta )
      {
        delta = dist;
        point = pymin;
      }
    }
    pymin = py;
  }

  return point->second;
}


// // Find the bin corresponding to position (x,y) over the phase space.
// int Binning::bin( const float& x, const float& y, const float& z ) const
// {
//   if ( x < y )
//     return - bin( y, x, z );

//   unsigned clbin = 0;

//   float mindist = std::numeric_limits< float >::max();
//   float dist;

//   // Brute-force search over all the reference bins.
//   for ( auto xybin : _xybin )
//     if ( ( dist = distance( xybin.first, first( x, y, z ), second( x, y, z ) ) ) < mindist )
//     {
//       clbin = xybin.second;
//       mindist = dist;
//     }

//   return _from0 + clbin;
// }

