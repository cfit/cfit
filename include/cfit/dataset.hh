#ifndef __DATASET_HH__
#define __DATASET_HH__

#include <string>
#include <map>
#include <vector>
#include <utility>

class Dataset
{
private:
  std::map< std::string, std::vector< std::pair< double, double > > > _data;

public:
  Dataset()  {};
  ~Dataset() {};

  void push( const std::string& field, const double& value, const double& error = 0. );

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

