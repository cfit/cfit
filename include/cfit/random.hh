#ifndef __RANDOM_HH__
#define __RANDOM_HH__

#include <random>

class Random
{
private:
  static std::uniform_real_distribution< double > _uniform;
  static std::default_random_engine _e;

public:
  static std::default_random_engine& engine()
  {
    return _e;
  }

  static void setSeed( const unsigned& seed ) { _e.seed( seed ); }

  static const double flat   ( const double& min = 0.0, const double& max = 1.0 )
  {
    return min + ( max - min ) * _uniform( engine() );
  }
  static const double uniform( const double& min = 0.0, const double& max = 1.0 )
  {
    return min + ( max - min ) * _uniform( engine() );
  }
};


#endif

