#ifndef __PDFMODEL_HH__
#define __PDFMODEL_HH__

#include <string>
#include <vector>

#include <cfit/exceptions.hh>
#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/pdfbase.hh>
#include <cfit/pdfexpr.hh>
#include <cfit/region.hh>


class FunctionMinimum;

class PdfModel : public PdfBase
{
  friend class PdfExpr;

protected:
  std::vector< std::string > _varOrder;
  std::vector< std::string > _parOrder;

  void push( const Variable&  var ) throw( PdfException );
  void push( const Parameter& par ) throw( PdfException );

  const Variable&  getVar( const Variable&  var ) const;
  const Variable&  getVar( const int&       idx ) const;
  const Parameter& getPar( const Parameter& par ) const;
  const Parameter& getPar( const int&       idx ) const;

  virtual const double area( const double& min, const double& max ) const throw( PdfException );

public:
  virtual PdfModel* copy() const = 0;

  void setVar ( const std::string& name, const double& val, const double& err = -1. ) throw( PdfException );
  void setPar ( const std::string& name, const double& val, const double& err = -1. ) throw( PdfException );

  virtual void setVars( const std::vector< double >&              vars ) throw( PdfException );
  virtual void setPars( const std::vector< double >&              pars ) throw( PdfException );
  virtual void setVars( const std::map< std::string, Variable  >& vars ) throw( PdfException );
  virtual void setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException );
  virtual void setPars( const FunctionMinimum&                    min  ) throw( PdfException );

  virtual       void   cache() {};
  virtual const double evaluate()                                    const throw( PdfException ) = 0;
  virtual const double evaluate( const std::vector< double >& vars ) const throw( PdfException ) = 0;
  virtual const double evaluate( const double& value )               const throw( PdfException ) // For pdfs of a single variable.
  {
    throw PdfException( "PdfModel::evaluate: evaluate( value ) has been called on a pdf with more than one variable." );
  }

  virtual const std::map< std::string, double > generate()           const throw( PdfException )
  {
    throw PdfException( "Generate error: attempting to generate with a model without generate() implementation" );
  }

  // Projection function for pdfs with one single variable.
  virtual const double project( const std::string& varName, const double& value ) const throw( PdfException )
  {
    if ( this->dependsOn( varName ) )
      return this->evaluate( value );

    return 1.0;
  }


  // Projection function for pdfs with one single variable.
  virtual const double project( const std::string& varName, const double& value, const Region& region ) const throw( PdfException )
  {
    if ( this->dependsOn( varName ) )
      return this->evaluate( value );

    const   std::map< const std::string, std::pair< double, double > > limits = region.limits();
    typedef std::map< const std::string, std::pair< double, double > > spMap;
    typedef spMap::const_iterator                                      spmIter;

    for ( spmIter limit = limits.begin(); limit != limits.end(); ++limit )
      if ( this->dependsOn( limit->first ) )
        return this->area( limit->second.first, limit->second.second );

    return 1.0;
  }


  // Projection function for pdfs with two variables.
  virtual const double project( const std::string& var1, const std::string& var2,
                                const double&      val1, const double&      val2 ) const throw( PdfException )
  {
    if ( this->dependsOn( var1 ) && this->dependsOn( var2 ) )
    {
      std::vector< double > vals;
      vals.push_back( val1 );
      vals.push_back( val2 );

      return this->evaluate( vals );
    }

    return 1.0;
  }


  // Projection function for pdfs with two variables.
  virtual const double project( const std::string& var1, const std::string& var2,
                                const double&      val1, const double&      val2,
                                const Region&      region                         ) const throw( PdfException )
  {
    if ( this->dependsOn( var1 ) && this->dependsOn( var2 ) )
    {
      std::vector< double > vals;
      vals.push_back( val1 );
      vals.push_back( val2 );

      return this->evaluate( vals );
    }

    const   std::map< const std::string, std::pair< double, double > > limits = region.limits();
    typedef std::map< const std::string, std::pair< double, double > > spMap;
    typedef spMap::const_iterator                                      spmIter;

    for ( spmIter limit = limits.begin(); limit != limits.end(); ++limit )
      if ( this->dependsOn( limit->first ) )
        return this->area( limit->second.first, limit->second.second );

    return 1.0;
  }

  friend const PdfExpr operator+( const PdfModel&      left, const PdfModel&      right );
  friend const PdfExpr operator*( const PdfModel&      left, const PdfModel&      right );

  friend const PdfExpr operator+( const PdfModel&      left, const PdfExpr&       right );
  friend const PdfExpr operator*( const PdfModel&      left, const PdfExpr&       right );

  friend const PdfExpr operator+( const PdfExpr&       left, const PdfModel&      right );
  friend const PdfExpr operator*( const PdfExpr&       left, const PdfModel&      right );

  friend const PdfExpr operator*( const Parameter&     left, const PdfModel&      right );
  friend const PdfExpr operator*( const PdfModel&      left, const Parameter&     right );
  friend const PdfExpr operator/( const PdfModel&      left, const Parameter&     right );

  friend const PdfExpr operator*( const ParameterExpr& left, const PdfModel&      right );
  friend const PdfExpr operator*( const PdfModel&      left, const ParameterExpr& right );
  friend const PdfExpr operator/( const PdfModel&      left, const ParameterExpr& right );

  friend const PdfExpr operator*( const double&        left, const PdfModel&      right );
  friend const PdfExpr operator*( const PdfModel&      left, const double&        right );
  friend const PdfExpr operator/( const PdfModel&      left, const double&        right );
};

#endif
