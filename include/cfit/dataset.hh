#ifndef __DATASET_HH__
#define __DATASET_HH__

#include <string>
#include <map>
#include <vector>
#include <utility>

#include <cfit/exceptions.hh>
#include <cfit/region.hh>

class Dataset
{
private:
  typedef std::map< std::string, std::vector< std::pair< double, double > > > data_type;
  typedef std::map< std::string, std::pair< double, double > >                datum_type;
  data_type _data;

public:
  Dataset()  {};
  ~Dataset() {};

  void push( const std::string& field, const double& value, const double& error = 0. );
  void push( const std::map< std::string, double >& event ); // Map of fields and values.
  void push( const std::map< std::string, std::pair< double, double > >& event );

  // Getters.
  bool                       empty ()                                      const;
  std::size_t                size  ()                                      const;
  const datum_type           entry ( const std::size_t& index )            const;
  double                     value ( const std::string& field, int entry ) const throw( DataException );
  double                     error ( const std::string& field, int entry ) const throw( DataException );
  std::vector< double >      values( const std::string& field )            const throw( DataException );
  std::vector< double >      errors( const std::string& field )            const throw( DataException );
  std::vector< std::string > fields()                                      const;
//void                       dump  ()                                      const;
  const Dataset              slice( const Region& region )                 const;

#ifdef MPI_ON
  // Scatter the data through all the processes in an MPI communicator.
  void scatter();
#endif
};

#endif

