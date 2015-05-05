
#include <cfit/random.hh>

std::mt19937_64                          Random::_engine  = std::mt19937_64();

std::uniform_real_distribution< double > Random::_uniform = std::uniform_real_distribution< double >();
std::normal_distribution      < double > Random::_normal  = std::normal_distribution      < double >();

