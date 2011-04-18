#ifndef __PDFMODEL_HH__
#define __PDFMODEL_HH__

#include <string>
#include <vector>

#include <cfit/pdfexception.hh>
#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/pdfbase.hh>
#include <cfit/pdf.hh>


class Pdf;


class PdfModel : public PdfBase
{
protected:
  std::vector< std::string > _varOrder;
  std::vector< std::string > _parOrder;

  void push( const Variable&  var ) throw( PdfException );
  void push( const Parameter& par ) throw( PdfException );

  const Variable&  getVar( int index ) const { return _vars.find( _varOrder[ index ] )->second; }
  const Parameter& getPar( int index ) const { return _pars.find( _parOrder[ index ] )->second; }

public:
  void setVars( const std::vector< double >& vars ) throw( PdfException );
  void setPars( const std::vector< double >& pars ) throw( PdfException );

  virtual void   cache() {};
  virtual double evaluate()                                    const throw( PdfException ) = 0;
  virtual double evaluate( const std::vector< double >& vars ) const throw( PdfException ) = 0;

  // ALERTA AMB ELS OPERADORS. NO S'HAURIA DE PODER FER gauss1 += gauss2, O HAURIA D'ESTAR BEN CONTROLAT.
  friend Pdf& operator+( const PdfModel& left, const PdfModel& right );
  friend Pdf& operator*( const PdfModel& left, const PdfModel& right );

  friend Pdf& operator+( const Pdf& left, const PdfModel& right );
  friend Pdf& operator*( const Pdf& left, const PdfModel& right );

  friend Pdf& operator+( const Pdf& left, const Pdf& right );
  friend Pdf& operator*( const Pdf& left, const Pdf& right );
  friend Pdf& operator^( const Pdf& left, const Pdf& right ) throw( PdfException );

  friend Pdf& operator+( const Parameter& left, const PdfModel& right );
  friend Pdf& operator*( const Parameter& left, const PdfModel& right );
};

#endif
