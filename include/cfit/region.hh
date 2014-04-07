#ifndef __REGION_HH__
#define __REGION_HH__


#include <map>
#include <string>
#include <utility>

class Region
{
private:
  typedef std::map< const std::string, std::pair< double, double > > spMap;

  spMap _limits;

public:
  Region() {}

  // Setters.
  void setLimit( const std::string& name, const double& lower, const double& upper )
  {
    _limits[ name ] = std::make_pair( lower, upper );
  }

  // Getters.
  const spMap&                       limits     ()                          const { return _limits;               }
  const bool                         hasLimit   ( const std::string& name ) const { return _limits.count( name ); }
  const std::vector< std::string >   limitedPars()                          const
  {
    std::vector< std::string > names;

    typedef std::map< const std::string, std::pair< double, double > >::const_iterator mIter;
    for ( mIter limit = _limits.begin(); limit != _limits.end(); ++limit )
      names.push_back( limit->first );

    return names;
  }

  const std::pair< double, double >& limit      ( const std::string& name ) const
  {
    return _limits.at( name );
  }
};

#endif

