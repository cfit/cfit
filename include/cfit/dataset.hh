#ifndef __DATASET_HH__
#define __DATASET_HH__

#include <map>
#include <vector>
#include <utility>

#include <cfit/variable.hh>

class SelectValue
{
private:
  std::string _name;

public:
  SelectValue( const std::string& name )
    : _name( name )
  {}

  double operator()( const std::map< std::string, std::pair< double, double > >& container ) const
  {
    return container.find( _name )->second.first;
  }
};


class Dataset
{
private:
  std::vector< std::map< std::string, std::pair< double, double > > > _data;
//   std::map< std::string, std::vector< std::pair< double, double > > > _data;

public:
  Dataset()  {};
  ~Dataset() {};

  void push( const std::string& field, const double value, const double error = 0. );

  // Getters.
  int                        size  ()                                      const;
  double                     value ( const std::string& field, int entry ) const;
  double                     error ( const std::string& field, int entry ) const;
  std::vector< double >      values( const std::string& field )            const;
  std::vector< double >      errors( const std::string& field )            const;
  std::vector< std::string > fields()                                      const;
//void                       dump  ()                                      const;

#ifdef MPI_ON
  // Scatter the data through all the processes in an MPI communicator.
  void scatter();
#endif
};

#endif

