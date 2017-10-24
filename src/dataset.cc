
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


// Add event from map of fields and values.
void Dataset::push( const std::map< std::string, double >& event )
{
  typedef std::map< std::string, double >::const_iterator fIter;
  for ( fIter entry = event.begin(); entry != event.end(); ++entry )
    _data[ entry->first ].push_back( std::make_pair( entry->second, 0.0 ) );
}


// Add event from map of fields and pair of values and errors.
void Dataset::push( const datum_type& event )
{
  for ( datum_type::const_iterator entry = event.begin(); entry != event.end(); ++entry )
    _data[ entry->first ].push_back( entry->second );
}


// Getters.
bool Dataset::empty() const
{
  return _data.empty();
}

std::size_t Dataset::size() const
{
  if ( _data.empty() )
    return 0;

  return _data.begin()->second.size();
}


const Dataset::datum_type Dataset::entry( const std::size_t& index ) const
{
  Dataset::datum_type ret;

  for ( std::map< std::string, std::vector< std::pair< double, double > > >::const_iterator b = _data.begin();
        b != _data.end(); ++b )
    ret.emplace( b->first, b->second.at( index ) );

  return ret;
}



double Dataset::value( const std::string& field, int entry ) const throw( DataException )
{
  if ( ! _data.count( field ) )
    throw DataException( "Dataset: requested variable " + field + " does not exist in dataset" );

  return _data.find( field )->second[ entry ].first;
}


double Dataset::error( const std::string& field, int entry ) const throw( DataException )
{
  if ( ! _data.count( field ) )
    throw DataException( "Dataset: requested variable " + field + " does not exist in dataset" );

  return _data.find( field )->second[ entry ].second;
}


std::vector< double > Dataset::values( const std::string& field ) const throw( DataException )
{
  std::vector< double > vals;

  if ( ! _data.count( field ) )
    throw DataException( "Dataset: requested variable " + field + " does not exist in dataset" );

  const std::vector< std::pair< double, double > >& entries = _data.find( field )->second;

  std::transform( entries.begin(), entries.end(), std::back_inserter( vals ), Select1st() );

  return vals;
}


std::vector< double > Dataset::errors( const std::string& field ) const throw( DataException )
{
  std::vector< double > errs;

  if ( ! _data.count( field ) )
    throw DataException( "Dataset: requested variable " + field + " does not exist in dataset" );

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



const Dataset Dataset::slice( const Region& region ) const
{
  std::map< const std::string, std::pair< double, double > > limits = region.limits();

  Dataset ret;

  bool accept; // To decide whether a given entry passes the region cuts.
  std::size_t nentries = this->size();
  for ( std::size_t e = 0; e < nentries; ++e )
  {
    // Check if this entry passes all the region cuts.
    const datum_type& entry = this->entry( e );
    accept = true;
    for ( std::map< const std::string, std::pair< double, double > >::const_iterator limit = limits.begin();
          accept && ( limit != limits.end() ); ++limit )
    {
      accept &= entry.at( limit->first ).first > limit->second.first;
      accept &= entry.at( limit->first ).first < limit->second.second;
    }

    // If this entry passes all the region cuts, include it in the result dataset.
    if ( accept )
      ret.push( entry );
  }

  return ret;
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

