#ifndef __DATASET_HH__
#define __DATASET_HH__

#include <string>
#include <map>
#include <vector>
#include <utility>

#include <cfit/exceptions.hh>

class Dataset
{
private:
  std::map< std::string, std::vector< std::pair< double, double > > > _data;

public:
  Dataset()  {};
  ~Dataset() {};

  void push( const std::string& field, const double& value, const double& error = 0. );
  void push( const std::map< std::string, double >& event ); // Map of fields and values.

  // Getters.
  std::size_t                size  ()                                      const;
  double                     value ( const std::string& field, int entry ) const throw( DataException );
  double                     error ( const std::string& field, int entry ) const throw( DataException );
  std::vector< double >      values( const std::string& field )            const throw( DataException );
  std::vector< double >      errors( const std::string& field )            const throw( DataException );
  std::vector< std::string > fields()                                      const;
//void                       dump  ()                                      const;

#ifdef MPI_ON
  // Scatter the data through all the processes in an MPI communicator.
  void scatter();
#endif
};

#endif

