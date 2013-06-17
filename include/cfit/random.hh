
#include <random>

class Random
{
public:
  static std::default_random_engine& engine()
  {
    static std::default_random_engine e{};
    return e;
  }
};

