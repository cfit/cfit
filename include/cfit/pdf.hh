#ifndef __PDF_HH__
#define __PDF_HH__

#include <sstream>

#include <string>
#include <vector>
#include <map>

#include <cfit/pdfexception.hh>
#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/pdfbase.hh>
#include <cfit/pdfmodel.hh>


class PdfModel;

class Pdf : public PdfBase
{
protected:
  std::vector< PdfModel*   > _pdfVect;
  std::vector< std::string > _parVect;

  std::string _expression;

public:
  void setVars( const std::vector< double >& vars ) throw( PdfException );
  void setPars( const std::vector< double >& pars ) throw( PdfException );

  void   cache();
  double evaluate()                                    const throw( PdfException );
  double evaluate( const std::vector< double >& vars ) const throw( PdfException );

  std::vector< std::string > commonVars() const throw( PdfException );

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
