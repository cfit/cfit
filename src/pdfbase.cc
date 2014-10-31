
#include <vector>
#include <string>
#include <algorithm>

#include <cfit/pdfbase.hh>


unsigned PdfBase::_cacheIdxReal    = 0;
unsigned PdfBase::_cacheIdxComplex = 0;


const std::vector< std::string > PdfBase::varNames() const
{
  std::vector< std::string > varNames;
  std::transform( _varMap.begin(), _varMap.end(), std::back_inserter( varNames ), Select1st() );
  return varNames;
}

