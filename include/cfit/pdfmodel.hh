#ifndef __PDFMODEL_HH__
#define __PDFMODEL_HH__

#include <string>
#include <vector>

#include <cfit/pdfexception.hh>
#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/pdfbase.hh>
#include <cfit/pdf.hh>


class FunctionMinimum;

class PdfModel : public PdfBase
{
  friend class Pdf;

protected:
  std::vector< std::string > _varOrder;
  std::vector< std::string > _parOrder;

  void push( const Variable&  var ) throw( PdfException );
  void push( const Parameter& par ) throw( PdfException );

  const Variable&  getVar( int index ) const;
  const Parameter& getPar( int index ) const;

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

  // ALERTA AMB ELS OPERADORS. NO S'HAURIA DE PODER FER gauss1 += gauss2, O HAURIA D'ESTAR BEN CONTROLAT.
  friend const Pdf operator+( const PdfModel&      left, const PdfModel&      right );
  friend const Pdf operator*( const PdfModel&      left, const PdfModel&      right );

  friend const Pdf operator+( const PdfModel&      left, const Pdf&           right );
  friend const Pdf operator*( const PdfModel&      left, const Pdf&           right );

  friend const Pdf operator+( const Pdf&           left, const PdfModel&      right );
  friend const Pdf operator*( const Pdf&           left, const PdfModel&      right );

  friend const Pdf operator*( const Parameter&     left, const PdfModel&      right );
  friend const Pdf operator*( const PdfModel&      left, const Parameter&     right );
  friend const Pdf operator/( const PdfModel&      left, const Parameter&     right );

  friend const Pdf operator*( const ParameterExpr& left, const PdfModel&      right );
  friend const Pdf operator*( const PdfModel&      left, const ParameterExpr& right );
  friend const Pdf operator/( const PdfModel&      left, const ParameterExpr& right );

  friend const Pdf operator*( const double&        left, const PdfModel&      right );
  friend const Pdf operator*( const PdfModel&      left, const double&        right );
  friend const Pdf operator/( const PdfModel&      left, const double&        right );
};

#endif
