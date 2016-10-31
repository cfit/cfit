#ifndef __PDFEXPR_HH__
#define __PDFEXPR_HH__

#include <sstream>

#include <string>
#include <vector>
#include <map>

#include <cfit/exceptions.hh>
#include <cfit/parameter.hh>
#include <cfit/parameterexpr.hh>
#include <cfit/pdfbase.hh>
#include <cfit/pdfmodel.hh>
#include <cfit/operation.hh>


class FunctionMinimum;
class Region;

class PdfExpr : public PdfBase
{
private:
  std::string                  _expression;
  std::vector< Operation::Op > _opers;
  std::vector< double        > _ctnts;
  std::vector< Parameter     > _parms;
  std::vector< PdfModel*     > _pdfs;

  std::map< std::string, std::pair< double, double > > _limits;
  double _scale;

  // Clean up the content of all the PdfExpr containers.
  void clear();

  void append( const PdfModel&      model );
  void append( const PdfExpr&       pdf   );
  void append( const Parameter&     par   );
  void append( const ParameterExpr& expr  );
  void append( const double&        ctnt  );
  void append( const Operation::Op& oper  );

  template< class L, class R >
  PdfExpr( const L& left, const R& right, const Operation::Op& oper )
    : _scale( 1.0 )
  {
    append( left  );
    append( right );
    append( oper  );
  }

  const std::map< std::string, double > generateFull() const throw( PdfException );


  template < class T >
  static void insert(       std::map< unsigned, std::vector< T > >& cache ,
                      const std::map< unsigned, std::vector< T > >& cached )
  {
    cache.insert( cached.begin(), cached.end() );
  }

  const std::map< unsigned, std::vector<               double   > > cacheReal   ( const Dataset& data );
  const std::map< unsigned, std::vector< std::complex< double > > > cacheComplex( const Dataset& data );

public:
  PdfExpr() : _scale( 1.0 ) {};
  PdfExpr( const PdfModel& model )
    : _scale( 1.0 )
  {
    append( model );
  }

  PdfExpr( const ParameterExpr& expr )
    : _scale( 1.0 )
  {
    append( expr );
  }

  PdfExpr( const PdfExpr& right );

  PdfExpr* copy() const { return new PdfExpr( *this ); }

  ~PdfExpr();

  void setPar( const std::string& name, const double& val, const double& err = -1. ) throw( PdfException );

  void setPars( const std::vector< double >&              pars ) throw( PdfException );
  void setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException );
  void setPars( const FunctionMinimum&                    min  ) throw( PdfException );

  void         cache();
  const double evaluate( const std::vector< double >& vars ) const throw( PdfException );

  const std::map< std::string, double > generate() const throw( PdfException );
  const double evaluate( const std::vector< double                 >& vars  ,
                         const std::vector< double                 >& cacheR,
                         const std::vector< std::complex< double > >& cacheC  ) const throw( PdfException );

  const double project( const std::string& varName,
                        const double&      value    ) const throw( PdfException );
  const double project( const std::string& var1,
                        const std::string& var2,
                        const double&      val1,
                        const double&      val2     ) const throw( PdfException );

  const double project( const std::string& varName,
                        const double&      value  ,
                        const Region&      region   ) const throw( PdfException );
  const double project( const std::string& var1,
                        const std::string& var2,
                        const double&      val1,
                        const double&      val2,
                        const Region&      region   ) const throw( PdfException );


  // Not sure if this is a dirty hack.
  void setLimits( const Variable& var, const double& min, const double& max );

  const double yield() const;

  double area()                                                               const throw( PdfException );
  double area( const std::string& var, const double& min, const double& max ) const throw( PdfException );

  // Assignment operator.
  const PdfExpr& operator= ( const PdfModel&      right );
  const PdfExpr& operator= ( const PdfExpr&       right );
  const PdfExpr& operator= ( const ParameterExpr& right );

  // Assignment operators with pdf objects.
  const PdfExpr& operator+=( const PdfModel&      right ) throw( PdfException );
  const PdfExpr& operator*=( const PdfModel&      right ) throw( PdfException );
  const PdfExpr& operator+=( const PdfExpr&       right ) throw( PdfException );
  const PdfExpr& operator*=( const PdfExpr&       right ) throw( PdfException );

  // Assignment operators with parameter objects and constants.
  const PdfExpr& operator*=( const Parameter&     right );
  const PdfExpr& operator/=( const Parameter&     right );
  const PdfExpr& operator*=( const ParameterExpr& right );
  const PdfExpr& operator/=( const ParameterExpr& right );
  const PdfExpr& operator*=( const double&        right );
  const PdfExpr& operator/=( const double&        right );

  std::vector< std::string > commonVars() const throw( PdfException );

  // Binary operators that need access to this class.
  friend const PdfExpr operator+( const PdfModel&      left, const PdfModel&      right );
  friend const PdfExpr operator*( const PdfModel&      left, const PdfModel&      right );

  friend const PdfExpr operator+( const PdfModel&      left, const PdfExpr&       right );
  friend const PdfExpr operator*( const PdfModel&      left, const PdfExpr&       right );

  friend const PdfExpr operator+( const PdfExpr&       left, const PdfModel&      right );
  friend const PdfExpr operator*( const PdfExpr&       left, const PdfModel&      right );

  friend const PdfExpr operator+( const PdfExpr&       left, const PdfExpr&       right );
  friend const PdfExpr operator*( const PdfExpr&       left, const PdfExpr&       right );

  friend const PdfExpr operator^( const PdfExpr&       left, const PdfExpr&       right ) throw( PdfException );

  friend const PdfExpr operator*( const Parameter&     left, const PdfModel&      right );
  friend const PdfExpr operator*( const PdfModel&      left, const Parameter&     right );
  friend const PdfExpr operator/( const PdfModel&      left, const Parameter&     right );

  friend const PdfExpr operator*( const ParameterExpr& left, const PdfModel&      right );
  friend const PdfExpr operator*( const PdfModel&      left, const ParameterExpr& right );
  friend const PdfExpr operator/( const PdfModel&      left, const ParameterExpr& right );

  friend const PdfExpr operator*( const double&        left, const PdfModel&      right );
  friend const PdfExpr operator*( const PdfModel&      left, const double&        right );
  friend const PdfExpr operator/( const PdfModel&      left, const double&        right );

  friend const PdfExpr operator*( const Parameter&     left, const PdfExpr&       right );
  friend const PdfExpr operator*( const PdfExpr&       left, const Parameter&     right );
  friend const PdfExpr operator/( const PdfExpr&       left, const Parameter&     right );

  friend const PdfExpr operator*( const ParameterExpr& left, const PdfExpr&       right );
  friend const PdfExpr operator*( const PdfExpr&       left, const ParameterExpr& right );
  friend const PdfExpr operator/( const PdfExpr&       left, const ParameterExpr& right );

  friend const PdfExpr operator*( const double&        left, const PdfExpr&       right );
  friend const PdfExpr operator*( const PdfExpr&       left, const double&        right );
  friend const PdfExpr operator/( const PdfExpr&       left, const double&        right );
};

#endif
