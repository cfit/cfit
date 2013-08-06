
#include <cfit/random.hh>

std::default_random_engine Random::_engine = std::default_random_engine();

std::uniform_real_distribution< double > Random::_uniform = std::uniform_real_distribution< double >();
std::normal_distribution      < double > Random::_normal  = std::normal_distribution      < double >();

