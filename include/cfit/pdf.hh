#ifndef __PDF_HH__
#define __PDF_HH__

#include <sstream>

#include <string>
#include <vector>
#include <map>

#include <cfit/pdfexception.hh>
#include <cfit/parameter.hh>
#include <cfit/parameterexpr.hh>
#include <cfit/pdfbase.hh>
#include <cfit/pdfmodel.hh>
#include <cfit/operation.hh>


class Pdf : public PdfBase
{
private:
  std::string                  _expression;
  std::vector< Operation::Op > _opers;
  std::vector< double        > _consts;
  std::vector< Parameter     > _parVect;
  std::vector< PdfModel*     > _pdfVect;

  void append( const PdfModel&      model );
  void append( const Pdf&           pdf   );
  void append( const Parameter&     par   );
  void append( const ParameterExpr& expr  );
  void append( const double&        ctnt  );

  template< class L, class R >
  Pdf( const L& left, const R& right, const Operation::Op& oper )
  {
    append( left  );
    append( right );

    _expression += "b"; // b = binary operation.
    _opers.push_back( oper );
  }

  // Execute binary operations.
  static double operate( const double& x, const double& y, const Operation::Op& oper ) throw( PdfException );

  // Execute unary operations.
  static double operate( const double& x,                  const Operation::Op& oper ) throw( PdfException );

public:
  Pdf() {};

  void setVars( const std::vector< double >& vars ) throw( PdfException );
  void setPars( const std::vector< double >& pars ) throw( PdfException );

  void   cache();
  double evaluate()                                    const throw( PdfException );
  double evaluate( const std::vector< double >& vars ) const throw( PdfException );

  std::vector< std::string > commonVars() const throw( PdfException );

  // ALERTA AMB ELS OPERADORS. NO S'HAURIA DE PODER FER gauss1 += gauss2, O HAURIA D'ESTAR BEN CONTROLAT.
  friend const Pdf operator+( const PdfModel&      left, const PdfModel&      right );
  friend const Pdf operator*( const PdfModel&      left, const PdfModel&      right );

  friend const Pdf operator+( const PdfModel&      left, const Pdf&           right );
  friend const Pdf operator*( const PdfModel&      left, const Pdf&           right );

  friend const Pdf operator+( const Pdf&           left, const PdfModel&      right );
  friend const Pdf operator*( const Pdf&           left, const PdfModel&      right );

  friend const Pdf operator+( const Pdf&           left, const Pdf&           right );
  friend const Pdf operator*( const Pdf&           left, const Pdf&           right );

  friend const Pdf operator^( const Pdf&           left, const Pdf&           right ) throw( PdfException );

  friend const Pdf operator*( const Parameter&     left, const PdfModel&      right );
  friend const Pdf operator*( const PdfModel&      left, const Parameter&     right );
  friend const Pdf operator/( const PdfModel&      left, const Parameter&     right );

  friend const Pdf operator*( const ParameterExpr& left, const PdfModel&      right );
  friend const Pdf operator*( const PdfModel&      left, const ParameterExpr& right );
  friend const Pdf operator/( const PdfModel&      left, const ParameterExpr& right );

  friend const Pdf operator*( const double&        left, const PdfModel&      right );
  friend const Pdf operator*( const PdfModel&      left, const double&        right );
  friend const Pdf operator/( const PdfModel&      left, const double&        right );

  friend const Pdf operator*( const Parameter&     left, const Pdf&           right );
  friend const Pdf operator*( const Pdf&           left, const Parameter&     right );
  friend const Pdf operator/( const Pdf&           left, const Parameter&     right );

  friend const Pdf operator*( const ParameterExpr& left, const Pdf&           right );
  friend const Pdf operator*( const Pdf&           left, const ParameterExpr& right );
  friend const Pdf operator/( const Pdf&           left, const ParameterExpr& right );

  friend const Pdf operator*( const double&        left, const Pdf&           right );
  friend const Pdf operator*( const Pdf&           left, const double&        right );
  friend const Pdf operator/( const Pdf&           left, const double&        right );
};

#endif
