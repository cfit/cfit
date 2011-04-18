#ifndef __FUNCTORS_HH__
#define __FUNCTORS_HH__

struct Select1st
{
  template <typename T>
  typename T::first_type operator()( T container ) const
  {
    return container.first;
  }
};


struct Select2nd
{
  template <typename T>
  typename T::second_type operator()( T container ) const
  {
    return container.second;
  }
};

#endif

