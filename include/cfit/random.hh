#ifndef __RANDOM_HH__
#define __RANDOM_HH__

#include <random>

class Random
{
private:
  static std::default_random_engine               _engine;

  static std::uniform_real_distribution< double > _uniform;
  static std::normal_distribution      < double > _normal;

public:
  static std::default_random_engine& engine()
  {
    return _engine;
  }

  static void setSeed( const unsigned& seed ) { _engine.seed( seed ); }

  static const double flat   ( const double& min = 0.0, const double& max = 1.0 )
  {
    return min + ( max - min ) * _uniform( engine() );
  }

  static const double uniform( const double& min = 0.0, const double& max = 1.0 )
  {
    return min + ( max - min ) * _uniform( engine() );
  }

  static const double normal( const double& mu = 0.0, const double& sigma = 1.0 )
  {
    return mu + sigma * _normal( _engine );
  }
};


#endif

