#ifndef __PDFBASE_HH__
#define __PDFBASE_HH__

#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/exceptions.hh>
#include <cfit/functors.hh>


class PdfBase
{
protected:
  std::map< std::string, Variable  > _varMap;
  std::map< std::string, Parameter > _parMap;

public:
  virtual ~PdfBase() {};

  // Setters.
  virtual void setVars( const std::vector< double >& vars ) throw( PdfException ) = 0;
  virtual void setPars( const std::vector< double >& pars ) throw( PdfException ) = 0;

  virtual void setVar ( const std::string& name, const double& val, const double& err = -1. ) throw( PdfException ) = 0;
  virtual void setPar ( const std::string& name, const double& val, const double& err = -1. ) throw( PdfException ) = 0;

  // Getters.
  const std::map< std::string, Variable  >& getVars()  const { return _varMap;        }
  const std::map< std::string, Parameter >& getPars()  const { return _parMap;        }
  const unsigned                            nVars()    const { return _varMap.size(); }
  const unsigned                            nPars()    const { return _parMap.size(); }
  const std::vector< std::string >          varNames() const;

  // Before evaluating the pdf at all data points, cache anything common to
  //    all points (usually compute the norm).
  virtual void cache() = 0;

  virtual const double evaluate()                                    const throw( PdfException ) = 0;
  virtual const double evaluate( const std::vector< double >& vars ) const throw( PdfException ) = 0;
};

#endif
