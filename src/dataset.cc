
#include <utility>
#include <string>
#include <vector>
#include <algorithm>

#include <cfit/functors.hh>
#include <cfit/dataset.hh>

#ifdef MPI_ON
#include <mpi.h>
#endif

// Add event from field, value and error.
void Dataset::push( const std::string& field, const double& value, const double& error )
{
  _data[ field ].push_back( std::make_pair( value, error ) );
}


// Getters.
int Dataset::size() const
{
  return _data.begin()->second.size();
}


double Dataset::value( const std::string& field, int entry ) const
{
  return _data.find( field )->second[ entry ].first;
}


double Dataset::error( const std::string& field, int entry ) const
{
  return _data.find( field )->second[ entry ].second;
}


std::vector< double > Dataset::values( const std::string& field ) const
{
  std::vector< double > vals;
  const std::vector< std::pair< double, double > >& entries = _data.find( field )->second;

  std::transform( entries.begin(), entries.end(), std::back_inserter( vals ), Select1st() );

  return vals;
}


std::vector< double > Dataset::errors( const std::string& field ) const
{
  std::vector< double > errs;
  const std::vector< std::pair< double, double > >& entries = _data.find( field )->second;

  std::transform( entries.begin(), entries.end(), std::back_inserter( errs ), Select2nd() );

  return errs;
}


std::vector< std::string > Dataset::fields() const
{
  std::vector< std::string > fieldVect;

  std::transform( _data.begin(), _data.end(), std::back_inserter( fieldVect ), Select1st() );

  return fieldVect;
}



#ifdef MPI_ON
void Dataset::scatter()
{
  const MPI::Comm& world = MPI::COMM_WORLD;
  const int size = world.Get_size();
  const int rank = world.Get_rank();

  int root = 0;

  int nFields;
  int nData;    // Number of data this process is assigned.
  int nAllData; // Total number of data in all processes.

  char* fieldName;
  int   fieldLength;

  // The root process is the one that has read the data, so it knows about it.
  if ( rank == root )
    {
      nFields  = _data.size();
      nAllData = this->size();
    }

  // Broadcast the number of fields and data to be scattered.
  world.Bcast( &nFields , 1, MPI::INT, root );
  world.Bcast( &nAllData, 1, MPI::INT, root );

  // Compute how many events this process must receive from the root and allocate memory for them.
  nData = nAllData / size + ( rank < nAllData % size );

  // Arrays that will receive sent data.
  double* values = new double[ nData ];
  double* errors = new double[ nData ];

  if ( rank == root )
    {
      // Vectors of the size and offset of the dataset to be delivered to each process.
      int count [ size ];
      int offset[ size ];

      // Compute how many events to send each process, and their offsets.
      for ( int proc = 0; proc < size; proc++ )
	{
	  count [ proc ] =   nAllData / size + ( proc < nAllData % size );
	  offset[ proc ] = ( nAllData / size ) * rank + std::min( nAllData % size, rank );
	}

      double* allValues = new double[ nAllData ];
      double* allErrors = new double[ nAllData ];

      typedef std::map< std::string, std::vector< std::pair< double, double > > >::const_iterator dataIter;
      for ( dataIter field = _data.begin(); field != _data.end(); field++ )
	{
	  // Determine the name of the field and its length.
	  fieldName   = const_cast< char* >( field->first.c_str() );
	  fieldLength = field->first.size() + 1;

	  // Convert the C++ container objects to C arrays.
	  std::transform( field->second.begin(), field->second.end(), allValues, Select1st() );
	  std::transform( field->second.begin(), field->second.end(), allErrors, Select2nd() );

	  // Broadcast the name of the field.
	  world.Bcast( &fieldLength, 1          , MPI::INT , root );
	  world.Bcast( fieldName   , fieldLength, MPI::CHAR, root );

	  // Scatter the data corresponding to this field.
	  world.Scatterv( allValues, count, offset, MPI::DOUBLE, values, nData, MPI::DOUBLE, root );
	  world.Scatterv( allErrors, count, offset, MPI::DOUBLE, errors, nData, MPI::DOUBLE, root );

	  _data[ field->first ].clear();
	  for ( int datum = 0; datum < nData; datum++ )
	    this->push( field->first, values[ datum ], errors[ datum ] );
	}
    }
  else
    {
      for ( int field = 0; field < nFields; field++ )
	{
	  world.Bcast( &fieldLength, 1          , MPI::INT , root );
	  fieldName = new char[ fieldLength ];
	  world.Bcast( fieldName   , fieldLength, MPI::CHAR, root );

	  world.Scatterv( 0, 0, 0, MPI::DOUBLE, values, nData, MPI::DOUBLE, root );
	  world.Scatterv( 0, 0, 0, MPI::DOUBLE, errors, nData, MPI::DOUBLE, root );

	  _data[ fieldName ].clear();
	  for ( int datum = 0; datum < nData; datum++ )
	    this->push( fieldName, values[ datum ], errors[ datum ] );

	  delete[] fieldName;
	}
    }
}
#endif

