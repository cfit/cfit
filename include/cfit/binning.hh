#ifndef __BINNING_HH__
#define __BINNING_HH__


#include <vector>
#include <string>

class Binning
{
  typedef std::pair< std::pair< float, float >, unsigned > Datum;
  typedef std::vector< Datum >                             Data;

  class ByX
  {
  public:
    bool operator()( const Datum& p1, const Datum& p2 )
    {
      if ( p1.first.first == p2.first.first )
        return p1.first.second < p2.first.second;

      return p1.first.first < p2.first.first;
    }
  };

  class ByY
  {
  public:
    bool operator()( const Datum& p1, const Datum& p2 )
    {
      if ( p1.first.second == p2.first.second )
        return p1.first.first < p2.first.first;

      return p1.first.second < p2.first.second;
    }
  };

  class PosX
  {
  public:
    bool operator()( const Datum& p1, const float& p2 )
    {
      return p1.first.first < p2;
      // return std::fabs( p1.first.first  - p2.first.first );
    }

    bool operator()( const float& p1, const Datum& p2 )
    {
      return p1 < p2.first.first;
      // return std::fabs( p1.first.second - p2.first.second );
    }
  };

  class PosY
  {
  public:
    bool operator()( const Datum& p1, const float& p2 )
    {
      return p1.first.second < p2;
      // return std::fabs( p1.first.second - p2.first.second );
    }

    bool operator()( const float& p1, const Datum& p2 )
    {
      return p1 < p2.first.second;
      // return std::fabs( p1.first.second - p2.first.second );
    }
  };

  class Dist
  {
  private:
    float _x;
    float _y;
  public:
    Dist( const float& x, const float& y )
      : _x( x ), _y( y )
      {}

    bool operator()( const Datum& p1, const Datum& p2 )
    {
      // return distance( p1.first, _x, _y ) < distance( p2.first, _x, _y );

      return std::pow( p1.first.first  - _x, 2 ) + std::pow( p1.first.second - _y, 2 ) <
             std::pow( p2.first.first  - _x, 2 ) + std::pow( p2.first.second - _y, 2 );

      // return std::sqrt( std::pow( p.first.first  - _x, 2 ) +
      //                   std::pow( p.first.second - _y, 2 ) );
    }
  };

private:
  Data _xybinx;
  Data _xybiny;

  static float distance( const std::pair< float, float >& a, const std::pair< float, float >& b );
  static float distance( const std::pair< float, float >& a, const float& bx, const float& by );

public:
  Binning() = default;
  // Binning( const std::string& binsfile );
  Binning( const Data& xybin );

  std::size_t size() const { return _xybinx.size(); }

  // Find the bin corresponding to position (x,y) over the phase space.
  int bin( const float& x, const float& y ) const;
};

#endif

