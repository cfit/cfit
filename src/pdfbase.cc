
#include <vector>
#include <string>
#include <algorithm>

#include <cfit/pdfbase.hh>

const std::vector< std::string > PdfBase::varNames() const
{
  std::vector< std::string > varNames;
  std::transform( _varMap.begin(), _varMap.end(), std::back_inserter( varNames ), Select1st() );
  return varNames;
}

