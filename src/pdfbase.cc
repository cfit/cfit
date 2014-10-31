
#include <vector>
#include <string>
#include <algorithm>

#include <cfit/pdfbase.hh>


unsigned PdfBase::_cacheIdxReal    = 0;
unsigned PdfBase::_cacheIdxComplex = 0;


void PdfBase::fix( const std::string& name ) throw( PdfException )
{
  if ( ! _parMap.count( name ) )
    throw PdfException( "Cannot fix unexisting parameter " + name + "." );

  _parMap[ name ].fix();
}


void PdfBase::release( const std::string& name ) throw( PdfException )
{
  if ( ! _parMap.count( name ) )
    throw PdfException( "Cannot release unexisting parameter " + name + "." );

  _parMap[ name ].release();
}


const std::vector< std::string > PdfBase::varNames() const
{
  std::vector< std::string > varNames;
  std::transform( _varMap.begin(), _varMap.end(), std::back_inserter( varNames ), Select1st() );
  return varNames;
}

