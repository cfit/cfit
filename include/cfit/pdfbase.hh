#ifndef __PDFBASE_HH__
#define __PDFBASE_HH__

#include <string>
#include <vector>
#include <map>
#include <complex>
#include <algorithm>

#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/exceptions.hh>
#include <cfit/functors.hh>


class Dataset;
class Region;
class Minimizer;

class PdfBase
{
  friend Minimizer;

private:
  static unsigned _cacheIdxReal;
  static unsigned _cacheIdxComplex;

protected:
  // Make a dataset available to a pdf such that it can compute values to be cached.
  virtual const std::map< unsigned, std::vector< double > > cacheReal( const Dataset& data )
  {
    return std::map< unsigned, std::vector< double > >();
  }

  virtual const std::map< unsigned, std::vector< std::complex< double > > > cacheComplex( const Dataset& data )
  {
    return std::map< unsigned, std::vector< std::complex< double > > >();
  }

  std::map< std::string, Variable  > _varMap;
  std::map< std::string, Parameter > _parMap;

public:
  virtual ~PdfBase() {}

  // Setters.
  virtual void setVars( const std::vector< double >& vars ) throw( PdfException ) = 0;
  virtual void setPars( const std::vector< double >& pars ) throw( PdfException ) = 0;

  virtual void setVar ( const std::string& name, const double& val, const double& err = -1. ) throw( PdfException ) = 0;
  virtual void setPar ( const std::string& name, const double& val, const double& err = -1. ) throw( PdfException ) = 0;

  void fix    ( const std::string& name ) throw( PdfException );
  void release( const std::string& name ) throw( PdfException );

  // Getters.
  const std::map< std::string, Variable  >& getVars()  const { return _varMap;        }
  const std::map< std::string, Parameter >& getPars()  const { return _parMap;        }
  const unsigned                            nVars()    const { return _varMap.size(); }
  const unsigned                            nPars()    const { return _parMap.size(); }
  const std::vector< std::string >          varNames() const;
  const bool                                isFixed()  const
  {
    return std::all_of( _parMap.begin(), _parMap.end(),
                        []( const std::pair< std::string, Parameter >& par ){ return par.second.isFixed(); } );
  }

  const unsigned& nCachedReal()    const { return _cacheIdxReal;    }
  const unsigned& nCachedComplex() const { return _cacheIdxComplex; }

  virtual const unsigned& assignCacheIdxReal()    { return _cacheIdxReal++;    }
  virtual const unsigned& assignCacheIdxComplex() { return _cacheIdxComplex++; }

  // Before evaluating the pdf at all data points, cache anything common to
  //    all points (usually compute the norm).
  virtual void cache() = 0;

  // Evaluate functions.
  virtual const double evaluate()                                    const throw( PdfException ) = 0; // For variables already set. To be obsolete.
  virtual const double evaluate( const std::vector< double >& vars ) const throw( PdfException ) = 0; // For any pdf.
  virtual const double evaluate( const double& value )               const throw( PdfException )      // For pdfs of a single variable.
  {
    throw PdfException( "PdfBase::evaluate: evaluate( value ) has been called on a pdf with more than one variable." );
  }

  virtual const double evaluate( const std::vector< double                 >& vars,
                                 const std::vector< double                 >&     ,
                                 const std::vector< std::complex< double > >&       ) const throw( PdfException ) = 0;

  virtual const std::map< std::string, double > generate()           const throw( PdfException ) = 0;

  virtual const double project( const std::string& varName,
                                const double&      value    ) const throw( PdfException ) = 0;
  virtual const double project( const std::string& var1,
                                const std::string& var2,
                                const double&      val1,
                                const double&      val2     ) const throw( PdfException ) = 0;

  virtual const double project( const std::string& varName,
                                const double&      value  ,
                                const Region&      region   ) const throw( PdfException ) = 0;
  virtual const double project( const std::string& var1  ,
                                const std::string& var2  ,
                                const double&      val1  ,
                                const double&      val2  ,
                                const Region&      region   ) const throw( PdfException ) = 0;

  const bool dependsOn( const std::string& var ) const
  {
    const std::vector< std::string >& vars = varNames();
    return std::find( vars.begin(), vars.end(), var ) != vars.end();
  }
};

#endif
