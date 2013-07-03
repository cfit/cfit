
#include <cfit/random.hh>

std::uniform_real_distribution< double > Random::_uniform = std::uniform_real_distribution< double >();

std::default_random_engine Random::_e = std::default_random_engine();

